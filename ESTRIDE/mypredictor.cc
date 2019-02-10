#include <cassert>
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

void getPredStride(ForUpdate *U, uint64_t &predicted_value, uint64_t seq_no) {
  bool predstride = false;
  int B[NBWAYSTR];
  int TAG[NBWAYSTR];
  uint64_t pc = U->pc;

  // use a 3-way skewed-associative structure for the stride predictor.
  // find the B's and TAG's for the current pc.
  for (int i = 0; i < NBWAYSTR; i++) {
    // B[i] index in way i ; TAG[i] tag in way i;
    // Hashing
    B[i] = ((((pc) ^ (pc >> (2 * LOGSTR - i)) ^ (pc >> (LOGSTR - i)) ^
              (pc >> (3 * LOGSTR - i))) *
             NBWAYSTR) +
            i) %
           (NBWAYSTR * (1 << LOGSTR));
    int j = (NBWAYSTR - i);
    assert(j >= 0);

    // Orthogonal hashing
    TAG[i] = ((pc >> (LOGSTR - j)) ^ (pc >> (2 * LOGSTR - j)) ^
              (pc >> (3 * LOGSTR - j)) ^ (pc >> (4 * LOGSTR - j))) &
             ((1 << TAGWIDTHSTR) - 1);

    U->B[i] = B[i];
    U->TAGSTR[i] = TAG[i];
  }

  int STHIT = -1;
  for (int i = 0; i < NBWAYSTR; i++) {
    if (STR[B[i]].tag == TAG[i]) {
      // a hit in the stride table.
      STHIT = B[i];
      break;
    }
  }
  U->STHIT = STHIT;

  if (STHIT >= 0) {
    if (SafeStride >= 0) {
      uint64_t LastCommitedValue = STR[STHIT].LastValue;
      if (STR[STHIT].conf >= MAXCONFIDSTR / 4) {
        // Stride predictor has enough confidence.
        int inflight = 0;
        // compute the number of inflight instances of the instruction
        for (uint64_t i = seq_commit + 1; i < seq_no; i++) {
          inflight += (Update[i & (MAXINFLIGHT - 1)].pc == pc);
        }
        // Do signed addition as stride can be negative too!
        predicted_value =
            (uint64_t)((int64_t)LastCommitedValue +
                       ((inflight + 1) * ((int64_t)STR[STHIT].Stride)));
        predstride = true;
      }
    }
  }
  U->predstride = predstride;
}

bool getPrediction(uint64_t seq_no, uint64_t pc, uint8_t piece,
                   uint64_t &predicted_value) {
  ForUpdate *U;
  U = &Update[seq_no & (MAXINFLIGHT - 1)];
  U->pc = pc + piece;
  getPredStride(U, predicted_value, seq_no);
  return U->predstride;
}

// Function determining wheter to update confidence
bool strideupdateconf(ForUpdate *U, uint64_t actual_value, int actual_latency,
                      int stride) {
  // Some black magic in numbers
#define UPDATECONFSTR                                                          \
  (true &&                             \
   ((random() & ((1 << (NOTLLCMISS + NOTL2MISS + NOTL1MISS + 2 * MFASTINST +   \
                        2 * (U->INSTTYPE != loadInstClass))) -                 \
                 1)) == 0))
  return (UPDATECONFSTR &
          ((abs(stride) > 1) || (U->INSTTYPE != loadInstClass) ||
           ((stride == -1) & ((random() & 1) == 0)) ||
           ((stride == 1) & ((random() & 3) == 0))));
}

// Function to determine whether to allocate a missing pc
// into the predictor probabilisitically
bool StrideAllocateOrNot(ForUpdate *U, uint64_t actual_value,
                         int actual_latency) {
  #define LOGPBINVSTR 4
  bool X = false;
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
  // Fill in the boiler plate variables from the inflight instruction
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
    // If we had a hit in the table
    uint64_t LastValue = STR[STHIT].LastValue;
    // No need to consider inflight instructions here because we will change the last value.
    uint64_t Value =
        (uint64_t)((int64_t)LastValue + (int64_t)STR[STHIT].Stride);
    int64_t INTER = abs(2 * ((int64_t)actual_value - (int64_t)LastValue) - 1);

    uint64_t stridetoalloc =
        (INTER < (1 << LOGSTRIDE))
            ? (uint64_t)((int64_t)actual_value - (int64_t)LastValue)
            : 0;
    STR[STHIT].LastValue = actual_value;

    if (STR[STHIT].NotFirstOcc > 0) {
      // Stride is known
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
      // Stride is unknown
      if (stridetoalloc != 0 or true) {
        // Allocate only non-zero strides
        STR[STHIT].Stride = stridetoalloc;
      } else {
        // Stride zero so do not allocate
        STR[STHIT].Stride = 0xffff;
        STR[STHIT].conf = 0;
        STR[STHIT].u = 0;
      }
      // Now stride is known
      STR[STHIT].NotFirstOcc++;
    }
  } else {
    assert(!U->prediction_result);
    if (StrideAllocateOrNot(U, actual_value, actual_latency)) {
      // Now you can allocate
      int X = random() % NBWAYSTR;
      bool done = false;

      // Try to evict a non-confident entry if possible else
      // a non-useful entry.

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
}

void updatePredictor(uint64_t seq_no, uint64_t actual_addr,
                     uint64_t actual_value, uint64_t actual_latency) {
  ForUpdate *U;
  U = &Update[seq_no & (MAXINFLIGHT - 1)];
  if (U->todo == 1) {
    UpdateStridePred(U, actual_value, (int)actual_latency);
    U->todo = 0;
  }
  seq_commit = seq_no;
}

void speculativeUpdate(uint64_t seq_no, bool eligible, uint8_t prediction_result,
                       uint64_t pc, uint64_t next_pc, InstClass insn, uint8_t piece,
                       uint64_t src1, uint64_t src2, uint64_t src3, uint64_t dst) {
  ForUpdate *U;
  U = &Update[seq_no & (MAXINFLIGHT - 1)];
  if (eligible) {
    U->todo = 1;
    U->INSTTYPE = insn;
    U->NbOperand = (src1 != 0xdeadbeef) + (src2 != 0xdeadbeef) + (src3 != 0xdeadbeef);
    U->prediction_result = (prediction_result == 1);
    if (SafeStride < (1 << 15) - 1) {
        SafeStride++;
    }
    if (prediction_result != 2) {
      // We predicted something.
      assert(U->predstride);
      if (prediction_result) {
        // Prediction was correct
        if (SafeStride < (1 << 15) - 1) {
            SafeStride += 4 * (1 + (insn == loadInstClass));
        }
      } else {
        // Penalize on wrong prediction.
        SafeStride -= 1024;
      }
    }
  }
}

void beginPredictor(int argc_other, char **argv_other) {  }

void endPredictor() {  }
