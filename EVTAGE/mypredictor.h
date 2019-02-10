#include <deque>
#include <vector>

// 8KB
// 65378 bits
#define UWIDTH 2
#define LOGLDATA 7
#define LOGBANK 5
#define TAGWIDTH 11
#define BANKS 47
#define BANKDATA (1 << LOGLDATA)
#define BANKSIZE (1 << LOGBANK)
#define PREDSIZE (BANKS * BANKSIZE)

#define NHIST 7
int HL[NHIST + 1] = {0, 0, 1, 3, 6, 12, 18, 30};

#define WIDTHCONFID 3
#define MAXCONFID ((1 << WIDTHCONFID) - 1)
#define MAXU ((1 << UWIDTH) - 1)

// Global path history: 8 * 64 = 512 bit path history
static uint64_t gpath[8] = {0, 0, 0, 0, 0, 0, 0, 0};
// Target history
static uint64_t gtargeth = 0;

// Struct for value table
struct longdata {
  uint64_t data; // 64 - 7 = 57 bit data
  uint8_t u; // 2 bit useful counter
};
static longdata LDATA[3 * BANKDATA];

// struct for bank tables
struct vtentry {
  uint64_t hashpt; // LOGLDATA + 2 bits
  uint8_t conf;    // WIDTHCONFID bits
  uint16_t tag;    // TAGWIDTH bits
  uint8_t u;       // 2 bits
};

static vtentry Vtage[PREDSIZE];

#define MAXTICK 1024
// for managing replacement on the VTAGE entries
static int TICK; // 10 bits
// for tracking the last misprediction on VTAGE
static int LastMispVT = 0; // 8 bits

// index function for VTAGE (use the global path history): just a complex hash function
uint32_t gi(int i, uint64_t pc) {
  // resulting hash value
  uint64_t res = 0;

  // The fuck is this ??? It is same as min(64, HL[i])
  int hl = (HL[i] < 64) ? (HL[i] % 64) : 64;

  // Take lower hl bits of global path history
  uint64_t inter = (hl < 64) ? (((1 << hl) - 1) & gpath[0]) : gpath[0];
  // some hash function
  inter ^= (pc >> (i)) ^ (pc);
  for (int t = 0; t < 8; t++) {
    res ^= inter;
    inter ^= ((inter & 15) << 16);
    inter >>= (LOGBANK - ((NHIST - i + LOGBANK - 1) % (LOGBANK - 1)));
  }

  hl = (hl < (HL[NHIST] + 1) / 2) ? hl : ((HL[NHIST] + 1) / 2);

  // Take lower hl bits of global target history
  inter ^= (hl < 64) ? (((1 << hl) - 1) & gtargeth) : gtargeth;
  for (int t = 0; t <= hl / LOGBANK; t++) {
    res ^= inter;
    inter ^= ((inter & 15) << 16);
    inter >>= LOGBANK;
  }

  // 8 KB does not has history length >= 64
  // if (HL[i] >= 64) {
  //   int REMAIN = HL[i] - 64;
  //   hl = REMAIN;
  //   int PT = 1;

  //   while (REMAIN > 0) {

  //     inter ^= ((hl < 64) ? (((1 << hl) - 1) & gpath[PT]) : gpath[PT]);
  //     for (int t = 0; t < 8; t++) {
  //       res ^= inter;
  //       inter ^= ((inter & 15) << 16);

  //       inter >>= (LOGBANK - ((NHIST - i + LOGBANK - 1) % (LOGBANK - 1)));
  //     }
  //     REMAIN = REMAIN - 64;
  //     PT++;
  //   }
  // }
  // Mask it to the bank size.
  return ((uint32_t)res & (BANKSIZE - 1));
}

// tags for VTAGE: just another complex hash function "orthogonal" to the index function
uint32_t gtag(int i, uint64_t pc) {
  // The end result variable.
  uint64_t res = 0;

  // The fuck is this again ???
  int hl = (HL[i] < 64) ? (HL[i] % 64) : 64;

  // Take lower hl bits of global path history
  uint64_t inter = (hl < 64) ? (((1 << hl) - 1) & gpath[0]) : gpath[0];
  // some hash function
  inter ^= ((pc >> (i)) ^ (pc >> (5 + i)) ^ (pc));
  for (int t = 0; t < 8; t++) {
    res ^= inter;
    inter ^= ((inter & 31) << 14);
    inter >>= (LOGBANK - ((NHIST - i + LOGBANK - 2) % (LOGBANK - 1)));
  }

  // hash with global target history
  hl = (hl < (HL[NHIST] + 1) / 2) ? hl : ((HL[NHIST] + 1) / 2);
  // some hashing
  inter ^= ((hl < 64) ? (((1 << hl) - 1) & gtargeth) : gtargeth);
  for (int t = 0; t <= hl / TAGWIDTH; t++) {
    res ^= inter;
    inter ^= ((inter & 15) << 16);
    inter >>= TAGWIDTH;
  }

  // 8KB has history length less than 64 always.
  // if (HL[i] >= 64) {
  //   int REMAIN = HL[i] - 64;
  //   hl = REMAIN;
  //   int PT = 1;

  //   while (REMAIN > 0) {

  //     inter ^= ((hl < 64) ? (((1 << hl) - 1) & gpath[PT]) : gpath[PT]);
  //     for (int t = 0; t < 8; t++) {
  //       res ^= inter;
  //       inter ^= ((inter & 31) << 14);
  //       inter >>= (TAGWIDTH - (NHIST - i - 1));
  //     }
  //     REMAIN = REMAIN - 64;
  //     PT++;
  //   }
  // }

  // Mask it to tag size.
  return ((uint32_t)res & ((1 << TAGWIDTH) - 1));
}

////// for managing speculative state and forwarding information to the back-end
struct ForUpdate {
  bool predvtage;
  bool prediction_result;
  uint8_t todo;
  uint64_t pc;
  uint32_t GI[NHIST + 1];
  uint32_t GTAG[NHIST + 1];
  int STHIT;
  int HitBank;
  int8_t INSTTYPE;
  int8_t NbOperand;
};

#define MAXINFLIGHT 256
static ForUpdate Update[MAXINFLIGHT]; // there may be 256 instructions inflight
