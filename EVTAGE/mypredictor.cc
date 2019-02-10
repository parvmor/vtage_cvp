#include <inttypes.h>
#include <iostream>
#include <map>
#include <utility>
#include <cassert>

#include "../cvp_kit/cvp.h"
#include "mypredictor.h"

int seq_commit;

#define NOTLLCMISS (actual_latency < 150)
#define NOTL2MISS (actual_latency < 60)
#define NOTL1MISS (actual_latency < 12)
#define FASTINST (actual_latency == 1)
#define MFASTINST (actual_latency < 3)

void getPredVtage(ForUpdate *U, uint64_t &predicted_value) {
  // whether to predict?
  bool predvtage = false;
  uint64_t pc = U->pc;
  uint64_t PCindex = ((pc) ^ (pc >> 2) ^ (pc >> 5)) % PREDSIZE;
  uint64_t PCbank = (PCindex >> LOGBANK) << LOGBANK;
  for (int i = 1; i <= NHIST; i++) {
    U->GI[i] = (gi(i, pc) + (PCbank + (i << LOGBANK))) % PREDSIZE;
    U->GTAG[i] = gtag(i, pc);
  }
  U->GTAG[0] = (pc ^ (pc >> 4) ^ (pc >> TAGWIDTH)) & ((1 << TAGWIDTH) - 1);
  U->GI[0] = PCindex;
  U->HitBank = -1;

  // Find the match with longest path history
  for (int i = NHIST; i >= 0; i--) {
    // Compare the tag in each bank.
    if (Vtage[U->GI[i]].tag == U->GTAG[i]) {
      U->HitBank = i;
      break;
    }
  }

  // when a misprediction is encountered on VTAGE, we do not predict with
  // VTAGE for 128 instructions; does not bring significant speed-up, but
  // reduces #misprediction significantly: mispreds tend to be clustered
  if (LastMispVT >= 128)
    if (U->HitBank >= 0) {
      // Found a match here.
      int index = Vtage[U->GI[U->HitBank]].hashpt;
      if (index < 3 * BANKDATA) {
        // the hash and the data are both present
        predicted_value = LDATA[index].data;
        // predict only when we have enough confidence.
        predvtage = ((Vtage[U->GI[U->HitBank]].conf >= MAXCONFID));
      }
    }
  U->predvtage = predvtage;
  return;
}

bool getPrediction(uint64_t seq_no, uint64_t pc, uint8_t piece,
                   uint64_t &predicted_value) {

  ForUpdate *U;
  // find the struct for current inflight instruction.
  U = &Update[seq_no & (MAXINFLIGHT - 1)];
  // pc + piece seems to be unique here.
  U->pc = pc + piece;
  getPredVtage(U, predicted_value);
  return U->predvtage;
}

// whether to update confidence probabilistically.
bool vtageupdateconf(ForUpdate *U, uint64_t actual_value, int actual_latency) {

#define LOWVAL ((abs(2 * ((int64_t)actual_value) + 1) < (1 << 16)) + (actual_value == 0))

#define UPDATECONF                                                             \
  ((random() &                                                                 \
    (((1 << (LOWVAL + NOTLLCMISS + 2 * FASTINST + NOTL2MISS + NOTL1MISS +      \
             ((U->INSTTYPE != loadInstClass) || NOTL1MISS))) -                 \
      1))) == 0)

  switch (U->INSTTYPE) {
  case aluInstClass:
  case fpInstClass:
  case slowAluInstClass:
  case undefInstClass:
  case loadInstClass:
  case storeInstClass:
    return (UPDATECONF);
    break;
  case uncondIndirectBranchInstClass:
    return (true);
    break;
  default:
    return (false);
  };
}

// whether to update usefulness counter probabilistically.
bool VtageUpdateU(ForUpdate *U, uint64_t actual_value, int actual_latency) {

#define UPDATEU                                                                \
  ((!U->prediction_result) &&                                                  \
   ((random() &                                                                \
     ((1 << (LOWVAL + 2 * NOTL1MISS + (U->INSTTYPE != loadInstClass) +         \
             FASTINST +                                                        \
             2 * (U->INSTTYPE == aluInstClass) * (U->NbOperand < 2))) -        \
      1)) == 0))

  switch (U->INSTTYPE) {
  case aluInstClass:
  case fpInstClass:
  case slowAluInstClass:
  case undefInstClass:
  case loadInstClass:
  case storeInstClass:
    return (UPDATEU);
    break;
  case uncondIndirectBranchInstClass:
    return (true);
    break;
  default:
    return (false);
  };
}

bool VtageAllocateOrNot(ForUpdate *U, uint64_t actual_value, int actual_latency,
                        bool MedConf) {
  bool X = false;

  // some black magic here!
  switch (U->INSTTYPE) {
  case undefInstClass:
  case aluInstClass:
  case storeInstClass:
    if (((U->NbOperand >= 2) & ((random() & 15) == 0)) ||
        ((U->NbOperand < 2) & ((random() & 63) == 0)))
  case fpInstClass:
  case slowAluInstClass:
  case loadInstClass:
        X = (((
          random() & (
              (2 << ((U->INSTTYPE != loadInstClass) + LOWVAL + NOTLLCMISS + NOTL2MISS + 1 * NOTL1MISS + 2 * FASTINST)) - 1)
          ) == 0) ||
          MedConf);
    break;
  case uncondIndirectBranchInstClass:
    X = true;
    break;
  default:
    X = false;
  };

  return X;
}

void UpdateVtagePred(ForUpdate *U, uint64_t actual_value, int actual_latency) {

  bool MedConf = false;
  // Keep hash pointer out of range as of now.
  uint64_t HashData =
      ((actual_value ^ (actual_value >> 7) ^ (actual_value >> 13) ^
        (actual_value >> 21) ^ (actual_value >> 29) ^ (actual_value >> 34) ^
        (actual_value >> 43) ^ (actual_value >> 52) ^ (actual_value >> 57)) &
       (BANKDATA - 1)) +
      3 * BANKDATA;

  bool ShouldWeAllocate = true;
  if (U->HitBank != -1) {
    // there was  an  hitting entry in VTAGE
    uint64_t index = U->GI[U->HitBank];

    // the entry might disappear, so check tag
    if (Vtage[index].tag == U->GTAG[U->HitBank]) {
      uint64_t indindex = Vtage[index].hashpt;
      ShouldWeAllocate =
          ((indindex >= 3 * BANKDATA) && (indindex != HashData)) ||
          ((indindex < 3 * BANKDATA) && (LDATA[indindex].data != actual_value));
      if (!ShouldWeAllocate) {
        // the predicted result is satisfactory: either a good hash without
        // data, or a pointer on the correct data
        if (Vtage[index].conf < MAXCONFID)
          if (vtageupdateconf(U, actual_value, actual_latency))
            Vtage[index].conf++;

        if (Vtage[index].u < MAXU)
          if ((VtageUpdateU(U, actual_value, actual_latency)) ||
              (Vtage[index].conf == MAXCONFID))
            Vtage[index].u++;

        if (indindex < 3 * BANKDATA)
          if (LDATA[indindex].u < 3)
            if (Vtage[index].conf == MAXCONFID)
              LDATA[indindex].u++;

        if (indindex >= 3 * BANKDATA) {
          // allocate a data entry when confidence is reasonable
          if (Vtage[index].conf >= MAXCONFID - 1) {
            int X[3];
            for (int i = 0; i < 3; i++)
              X[i] = (((actual_value) ^ (actual_value >> (LOGLDATA + (i + 1))) ^
                       (actual_value >> (3 * (LOGLDATA + (i + 1)))) ^
                       (actual_value >> (4 * (LOGLDATA + (i + 1)))) ^
                       (actual_value >> (5 * (LOGLDATA + (i + 1)))) ^
                       (actual_value >> (6 * (LOGLDATA + (i + 1)))) ^
                       (actual_value >> (2 * (LOGLDATA + (i + 1))))) &
                      ((1 << LOGLDATA) - 1)) +
                     i * (1 << LOGLDATA);
            bool done = false;
            for (int i = 0; i < 3; i++) {
              if (LDATA[X[i]].data == actual_value) {
                // the data is already present
                Vtage[index].hashpt = X[i];
                done = true;
                break;
              }
            }
            if (!done)
              if ((random() & 3) == 0) {
                // data absent: let us try try to steal an entry
                int i = (((uint64_t)random()) % 3);
                bool done = false;
                for (int j = 0; j < 3; j++) {
                  if ((LDATA[X[i]].u == 0)) {
                    LDATA[X[i]].data = actual_value;
                    LDATA[X[i]].u = 1;
                    Vtage[index].hashpt = X[i];
                    done = true;
                    break;
                  }
                  i++;
                  i = i % 3;
                }
                if (U->INSTTYPE == loadInstClass)
                  if (!done) {
                    if ((LDATA[X[i]].u == 0)) {
                      LDATA[X[i]].data = actual_value;
                      LDATA[X[i]].u = 1;
                      Vtage[index].hashpt = X[i];
                    } else if ((random() & 31) == 0)
                      LDATA[X[i]].u--;
                  }
              }
          }
        }

      } else {
        Vtage[index].hashpt = HashData;
        if ((Vtage[index].conf > MAXCONFID / 2) ||
            ((Vtage[index].conf == MAXCONFID / 2) & (Vtage[index].u == 3)) ||
            ((Vtage[index].conf > 0) & (Vtage[index].conf < MAXCONFID / 2)))
          MedConf = true;

        if (Vtage[index].conf == MAXCONFID) {

          Vtage[index].u = (Vtage[index].conf == MAXCONFID);
          Vtage[index].conf -= (MAXCONFID + 1) / 4;
        } else {
          Vtage[index].conf = 0;
          Vtage[index].u = 0;
        }
      }
    }
  }

  if (!U->prediction_result)
    // Don't waste your time allocating if it is predicted by the other
    // component
    if (ShouldWeAllocate) {
      // avoid allocating too often
      if (VtageAllocateOrNot(U, actual_value, actual_latency, MedConf)) {
        int ALL = 0;
        int NA = 0;
        int DEP = (U->HitBank + 1) + ((random() & 7) == 0);
        if (U->HitBank == 0)
          DEP++;

        if (U->HitBank == -1) {
          if (random() & 7)
            DEP = random() & 1;
          else
            DEP = 2 + ((random() & 7) == 0);
        }

        if (DEP > 1) {

          for (int i = DEP; i <= NHIST; i++) {
            int index = U->GI[i];
            if ((Vtage[index].u == 0) &&
                ((Vtage[index].conf == MAXCONFID / 2) ||
                 (Vtage[index].conf <= (random() & MAXCONFID))))
            // slightly favors the entries with real confidence
            {
              Vtage[index].hashpt = HashData;
              Vtage[index].conf =
                  MAXCONFID /
                  2; // set to 3  for faster warming to  high confidence
              Vtage[index].tag = U->GTAG[i];
              ALL++;

              break;

            } else {
              NA++;
            }
          }
        } else {

          for (int j = 0; j <= 1; j++) {
            int i = (j + DEP) & 1;

            int index = U->GI[i];
            if ((Vtage[index].u == 0) &&
                ((Vtage[index].conf == MAXCONFID / 2) ||
                 (Vtage[index].conf <= (random() & MAXCONFID)))) {
              Vtage[index].hashpt = HashData;
              Vtage[index].conf = MAXCONFID / 2;
              if (U->NbOperand == 0)
                if (U->INSTTYPE == aluInstClass)
                  Vtage[index].conf = MAXCONFID;
              Vtage[index].tag = U->GTAG[i];
              ALL++;
              break;
            } else {
              NA++;
            }
          }
        }

        TICK += NA - (5 * ALL);
        if (TICK < 0)
          TICK = 0;
        if (TICK >= MAXTICK) {

          for (int i = 0; i < PREDSIZE; i++)
            if (Vtage[i].u > 0)
              Vtage[i].u--;
          TICK = 0;
        }
      }
    }
}

void updatePredictor(uint64_t seq_no, uint64_t actual_addr,
                     uint64_t actual_value, uint64_t actual_latency) {
  ForUpdate *U;
  // Find inflight struct.
  U = &Update[seq_no & (MAXINFLIGHT - 1)];
  // This parameter is set by speculative update if instruction is eligible for prediction.
  if (U->todo == 1) {
    UpdateVtagePred(U, actual_value, (int)actual_latency);
    U->todo = 0;
  }
  seq_commit = seq_no;
}

void speculativeUpdate(uint64_t seq_no, bool eligible, uint8_t prediction_result,
                       uint64_t pc, uint64_t next_pc, InstClass insn, uint8_t piece,
                       uint64_t src1, uint64_t src2, uint64_t src3, uint64_t dst) {

  // Why??
  // the framework does not allow to filter the predictions, so we predict every instruction

  ForUpdate *U;
  // Current inflight instruction struct
  U = &Update[seq_no & (MAXINFLIGHT - 1)];
  // Assume we predicted correctly and increment
  // Will be set to zero in case of misprediction
  LastMispVT++;

  if (eligible) {
    // Number of operands.
    U->NbOperand = (src1 != 0xdeadbeef) + (src2 != 0xdeadbeef) + (src3 != 0xdeadbeef);
    U->todo = 1;
    U->INSTTYPE = insn;
    U->prediction_result = (prediction_result == 1);
    // If we predicted incorrectly
    if (prediction_result == 0) {
      assert (U->predvtage);
      LastMispVT = 0;
    }
  }

  bool isCondBr = insn == condBranchInstClass;
  bool isUnCondBr = insn == uncondIndirectBranchInstClass ||
                    insn == uncondDirectBranchInstClass;
  // update path history
  if (isCondBr || isUnCondBr) {
    if (pc != next_pc - 4) {
      for (int i = 7; i > 0; i--) {
        gpath[i] = (gpath[i] << 1) ^ ((gpath[i - 1] >> 63) & 1);
      }
      gpath[0] = (gpath[0] << 1) ^ (pc >> 2);
      gtargeth = (gtargeth << 1) ^ (next_pc >> 2);
    }
  }
}

void beginPredictor(int argc_other, char **argv_other) { }

void endPredictor() { }
