#include <inttypes.h>
#include <iostream>
#include <map>
#include <utility>

#include "../cvp_kit/cvp.h"
#include "mypredictor.h"
int seq_commit;

#define NOTLLCMISS (actual_latency < 150)
#define NOTL2MISS (actual_latency < 60)
#define NOTL1MISS (actual_latency < 12)
#define FASTINST (actual_latency == 1)
#define MFASTINST (actual_latency < 3)

void getPredVtage(ForUpdate *U, uint64_t &predicted_value) {
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

  for (int i = NHIST; i >= 0; i--) {
    if (Vtage[U->GI[i]].tag == U->GTAG[i]) {
      U->HitBank = i;
      break;
    }
  }

  if (LastMispVT >= 128)
    // when a misprediction is encountered on VTAGE, we do not predict with
    // VTAGE for 128 instructions; does not bring significant speed-up, but
    // reduces the misprediction number significantly: mispredictions tend to be
    // clustered
    if (U->HitBank >= 0) {
      int index = Vtage[U->GI[U->HitBank]].hashpt;
      if (index < 3 * BANKDATA) {
        // the hash and the data are both present
        predicted_value = LDATA[index].data;
        predvtage = ((Vtage[U->GI[U->HitBank]].conf >= MAXCONFID));
      }
    }
  U->predvtage = predvtage;
}

void getPredStride(ForUpdate *U, uint64_t &predicted_value, uint64_t seq_no) {
  bool predstride = false;

  int B[NBWAYSTR];
  int TAG[NBWAYSTR];
  uint64_t pc = U->pc;
  // use a 3-way skewed-associative structure  for the stride predictor
  for (int i = 0; i < NBWAYSTR; i++) {
    // B[i] index in way i ; TAG[i] tag in way i;
    B[i] = ((((pc) ^ (pc >> (2 * LOGSTR - i)) ^ (pc >> (LOGSTR - i)) ^
              (pc >> (3 * LOGSTR - i))) *
             NBWAYSTR) +
            i) %
           (NBWAYSTR * (1 << LOGSTR));
    int j = (NBWAYSTR - i);
    if (j < 0)
      j = 0;

    TAG[i] = ((pc >> (LOGSTR - j)) ^ (pc >> (2 * LOGSTR - j)) ^
              (pc >> (3 * LOGSTR - j)) ^ (pc >> (4 * LOGSTR - j))) &
             ((1 << TAGWIDTHSTR) - 1);
    U->B[i] = B[i];
    U->TAGSTR[i] = TAG[i];
  }

  int STHIT = -1;
  for (int i = 0; i < NBWAYSTR; i++) {
    if (STR[B[i]].tag == TAG[i]) {
      STHIT = B[i];
      break;
    }
  }
  U->STHIT = STHIT;
  if (STHIT >= 0)
    if (SafeStride >= 0) { // hit
      uint64_t LastCommitedValue = STR[STHIT].LastValue;

      if (STR[STHIT].conf >= MAXCONFIDSTR / 4) {
        int inflight = 0;
        // compute the number of inflight instances of the instruction
        for (uint64_t i = seq_commit + 1; i < seq_no; i++) {
          inflight += (Update[i & (MAXINFLIGHT - 1)].pc == pc);
        }
        predicted_value =
            (uint64_t)((int64_t)LastCommitedValue +
                       ((inflight + 1) * ((int64_t)STR[STHIT].Stride)));
        predstride = true;
      }
    }
  U->predstride = predstride;
}

bool getPrediction(uint64_t seq_no, uint64_t pc, uint8_t piece,
                   uint64_t &predicted_value) {

  ForUpdate *U;
  U = &Update[seq_no & (MAXINFLIGHT - 1)];
  U->pc = pc + piece;
  U->predvtage = false;
  U->predstride = false;
#ifdef PREDSTRIDE
  getPredStride(U, predicted_value, seq_no);
#endif
#ifdef PREDVTAGE
  getPredVtage(U, predicted_value);
#endif
  return (U->predstride || U->predvtage);
}

// Update of the Stride predictor
// function determining whether to  update or not confidence on a correct
// prediction
bool strideupdateconf(ForUpdate *U, uint64_t actual_value, int actual_latency,
                      int stride) {
#define UPDATECONFSTR                                                          \
  (((!U->prediction_result) || (U->predstride)) &&                             \
   ((random() & ((1 << (NOTLLCMISS + NOTL2MISS + NOTL1MISS + 2 * MFASTINST +   \
                        2 * (U->INSTTYPE != loadInstClass))) -                 \
                 1)) == 0))
  return (UPDATECONFSTR &
          ((abs(stride) > 1) || (U->INSTTYPE != loadInstClass) ||
           ((stride == -1) & ((random() & 1) == 0)) ||
           ((stride == 1) & ((random() & 3) == 0))));
}

// Allocate or not if instruction absent from the predictor
bool StrideAllocateOrNot(ForUpdate *U, uint64_t actual_value,
                         int actual_latency) {
  bool X = false;
#define LOGPBINVSTR 4
  switch (U->INSTTYPE) {
  case aluInstClass:
  case storeInstClass:
    X = ((random() & ((1 << (LOGPBINVSTR + 2)) - 1)) == 0);
    break;
  case fpInstClass:
    X = ((random() & ((1 << (LOGPBINVSTR)) - 1)) == 0);
    break;
  case slowAluInstClass:
    X = ((random() & ((1 << (LOGPBINVSTR)) - 1)) == 0);
    break;
  case loadInstClass:
    X = ((random() &
          ((1 << (NOTLLCMISS + NOTL2MISS + NOTL1MISS + MFASTINST)) - 1)) == 0);
    break;
  };
  return (X);
}

void UpdateStridePred(ForUpdate *U, uint64_t actual_value, int actual_latency) {

  int B[NBWAYSTR];
  int TAG[NBWAYSTR];
  for (int i = 0; i < NBWAYSTR; i++) {
    B[i] = U->B[i];
    TAG[i] = U->TAGSTR[i];
  }
  int STHIT = -1;

  for (int i = 0; i < NBWAYSTR; i++) {
    if (STR[B[i]].tag == TAG[i]) {
      STHIT = B[i];
      break;
    }
  }

  if (STHIT >= 0) {
    uint64_t LastValue = STR[STHIT].LastValue;
    uint64_t Value =
        (uint64_t)((int64_t)LastValue + (int64_t)STR[STHIT].Stride);
    int64_t INTER = abs(2 * ((int64_t)actual_value - (int64_t)LastValue) - 1);

    uint64_t stridetoalloc =
        (INTER < (1 << LOGSTRIDE))
            ? (uint64_t)((int64_t)actual_value - (int64_t)LastValue)
            : 0;

    STR[STHIT].LastValue = actual_value;

    // special case when the stride is not determined
    if (STR[STHIT].NotFirstOcc > 0) {
      if (Value == actual_value) {

        if (STR[STHIT].conf < MAXCONFIDSTR) {
          if (strideupdateconf(U, actual_value, actual_latency,
                               (int)stridetoalloc))
            STR[STHIT].conf++;
        }

        if (STR[STHIT].u < 3)
          if (strideupdateconf(U, actual_value, actual_latency,
                               (int)stridetoalloc))
            STR[STHIT].u++;
        if (STR[STHIT].conf >= MAXCONFIDSTR / 4)
          STR[STHIT].u = 3;
      } else {
        // misprediction

        if (STR[STHIT].conf > (1 << (WIDTHCONFIDSTR - 3))) {
          STR[STHIT].conf -= (1 << (WIDTHCONFIDSTR - 3));
        } else {
          STR[STHIT].conf = 0;
           STR[STHIT].u = 0;
        }
        // this allows to restart a new sequence with a different stride
        STR[STHIT].NotFirstOcc = 0;
      }
    } else {
      if (stridetoalloc != 0 || true) {
        STR[STHIT].Stride = stridetoalloc;
      } else {
        STR[STHIT].Stride = 0xffff;
        STR[STHIT].conf = 0;
        STR[STHIT].u = 0;
      }
      STR[STHIT].NotFirstOcc++;
    }
  } else if (!U->prediction_result) {
    if (StrideAllocateOrNot(U, actual_value, actual_latency)) {
      int X = random() % NBWAYSTR;
      bool done = false;
      // the target entry is not a stride candidate
      for (int i = 0; i < NBWAYSTR; i++) {
        STHIT = B[X];
        if (STR[STHIT].conf == 0) {
          STR[STHIT].conf = 1; // just to allow not to ejected before testing if
                               // possible stride candidate
          STR[STHIT].u = 0;
          STR[STHIT].tag = TAG[X];
          STR[STHIT].Stride = 0;
          STR[STHIT].NotFirstOcc = 0;
          STR[STHIT].LastValue = actual_value;
          done = true;
          break;
        }
        X = (X + 1) % NBWAYSTR;
      }
      // the target entry has not been useful recently
      if (!done)
        for (int i = 0; i < NBWAYSTR; i++) {
          STHIT = B[X];
          if (STR[STHIT].u == 0) {
            STR[STHIT].conf = 1;
            STR[STHIT].u = 0;
            STR[STHIT].tag = TAG[X];
            STR[STHIT].Stride = 0;
            STR[STHIT].NotFirstOcc = 0;
            STR[STHIT].LastValue = actual_value;
            done = true;
            break;
          }
          X = (X + 1) % NBWAYSTR;
        }
      // if unable to allocate: age some target entry
      if (!done) {
        if ((random() & ((1 << (2 + 2 * (STR[STHIT].conf > (MAXCONFIDSTR) / 8) +
                                2 * (STR[STHIT].conf >= MAXCONFIDSTR / 4))) -
                         1)) == 0)
          STR[STHIT].u--;
      }
    }
  }
  //////////////////////////////
}

/////////Update of  VTAGE
// function determining whether to  update or not confidence on a correct
// prediction

bool vtageupdateconf(ForUpdate *U, uint64_t actual_value, int actual_latency) {
#define LOWVAL                                                                 \
  ((abs(2 * ((int64_t)actual_value) + 1) < (1 << 16)) + (actual_value == 0))

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

// Update of the U counter or not
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

  switch (U->INSTTYPE) {
  case undefInstClass:
  case aluInstClass:
  case storeInstClass:
    if (((U->NbOperand >= 2) & ((random() & 15) == 0)) ||
        ((U->NbOperand < 2) & ((random() & 63) == 0)))
    case fpInstClass:
    case slowAluInstClass:
    case loadInstClass:
        X = ((random() & ((2 << (
            ((U->INSTTYPE != loadInstClass) || (NOTL1MISS)) +
              LOWVAL + NOTLLCMISS + NOTL2MISS + NOTL1MISS + 2 * FASTINST)) - 1)) == 0) ||
           MedConf;

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
    // the entry may have dissappeared in the interval between the prediction
    // and  the commit

    if (Vtage[index].tag == U->GTAG[U->HitBank]) {
      //  update the prediction
      uint64_t indindex = Vtage[index].hashpt;
      ShouldWeAllocate =
          ((indindex >= 3 * BANKDATA) & (indindex != HashData)) ||
          ((indindex < 3 * BANKDATA) & (LDATA[indindex].data != actual_value));
      if (!ShouldWeAllocate) {
        // the predicted result is satisfactory: either a good hash without
        // data, or a pointer on the correct data
        ShouldWeAllocate = false;

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

          // try to allocate an effective data entry when the confidence has
          // reached a reasonable level
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
                    } else if ((random() & 3) == 0)
                      LDATA[X[i]].u--;
                  }
              }
          }
        }

      }

      else {
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

        TICK += (NA - (5 * ALL));
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
  U = &Update[seq_no & (MAXINFLIGHT - 1)];
  if (U->todo == 1) {
#ifdef PREDVTAGE
    UpdateVtagePred(U, actual_value, (int)actual_latency);
#endif
#ifdef PREDSTRIDE
    UpdateStridePred(U, actual_value, (int)actual_latency);
#endif
    U->todo = 0;
  }
  seq_commit = seq_no;
}

void speculativeUpdate(
    uint64_t seq_no, bool eligible, uint8_t prediction_result,
    uint64_t pc, uint64_t next_pc, InstClass insn, uint8_t piece,
    uint64_t src1, uint64_t src2, uint64_t src3, uint64_t dst) {

  ForUpdate *U;
  U = &Update[seq_no & (MAXINFLIGHT - 1)];

  LastMispVT++;

  if (eligible) {

    U->NbOperand =
        (src1 != 0xdeadbeef) + (src2 != 0xdeadbeef) + (src3 != 0xdeadbeef);
    U->todo = 1;
    U->INSTTYPE = insn;
    U->prediction_result = (prediction_result == 1);
    if (SafeStride < (1 << 15) - 1)
      SafeStride++;
    if (prediction_result != 2) {
      if (prediction_result) {
        if (U->predstride)
          if (SafeStride < (1 << 15) - 1)
            SafeStride += 4 * (1 + (insn == loadInstClass));
      } else {
        if (U->predvtage)
          LastMispVT = 0;
        if (U->predstride)
          SafeStride -= 1024;
      }
    }
  }

  bool isCondBr = insn == condBranchInstClass;
  bool isUnCondBr = insn == uncondIndirectBranchInstClass ||
                    insn == uncondDirectBranchInstClass;
  // path history
  // just to have a longer history without (software) management
  if ((isCondBr) || (isUnCondBr))
    if (pc != next_pc - 4) {
      for (int i = 7; i > 0; i--)
        gpath[i] = (gpath[i] << 1) ^ ((gpath[i - 1] >> 63) & 1);
      gpath[0] = (gpath[0] << 1) ^ (pc >> 2);
      gtargeth = (gtargeth << 1) ^ (next_pc >> 2);
    }
}

void beginPredictor(int argc_other, char **argv_other) { }

void endPredictor() { }
