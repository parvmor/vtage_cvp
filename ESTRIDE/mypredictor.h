#include <deque>
#include <vector>
#include <cassert>

// 8KB
// 4.026 //3.729 Stride only // 3.437 for TAGE  only
#define LOGSTR 4
#define NBWAYSTR 3
#define TAGWIDTHSTR 14
#define LOGSTRIDE 20

#define WIDTHCONFIDSTR 5
#define MAXCONFIDSTR ((1 << WIDTHCONFIDSTR) - 1)
#define UWIDTH 2
#define MAXU ((1 << UWIDTH) - 1)

#define MINSTRIDE -(1 << (LOGSTRIDE - 1))
#define MAXSTRIDE (-MINSTRIDE - 1)

// The E-Stride predictor
struct strdata {
  uint64_t LastValue;   // 64 bits
  uint64_t Stride;      // LOGSTRIDE bits
  uint8_t conf;         // WIDTHCONFIDSTR bits
  uint16_t tag;         // TAGWIDTHSTR bits
  uint16_t NotFirstOcc; // 1 bits
  int u;                // UWIDTH bits
  // 65 + UWIDTH + LOGSTRIDE + WIDTHCONFIDSTR + TAGWIDTHSTR bits
};

static strdata STR[NBWAYSTR * (1 << LOGSTR)]; // 48 entries

static int SafeStride = 0; // 16 bits

struct ForUpdate {
  bool predstride;
  bool prediction_result;
  uint8_t todo;
  uint64_t pc;
  int B[NBWAYSTR];
  int TAGSTR[NBWAYSTR];
  int STHIT;
  int8_t INSTTYPE;
  int8_t NbOperand;
};

#define MAXINFLIGHT 256
static ForUpdate Update[MAXINFLIGHT]; // there may be 256 instructions inflight
