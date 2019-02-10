#ifndef __PREDICTOR_H__
#define __PREDICTOR_H__

#include <inttypes.h>

#define USE_WIDTH (2)
#define USE_MAX ((1 << USE_WIDTH) - 1)

#define CONF_WIDTH (4)
#define CONF_MAX ((1 << CONF_WIDTH) - 1)
#define CONF_THRESHOLD (CONF_MAX)

#define STRIDE_LOG (25)
#define STRIDE_MIN (-(1 << (STRIDE_LOG - 1)))
#define STRIDE_MAX ((1 << (STRIDE_LOG - 1)) - 1)

#define TAGE_HIST_LEN (10)
static uint64_t tage_hist[TAGE_HIST_LEN] = { 0, 1, 3, 6, 13, 20, 30, 36, 45, 57 };

#define TAGE_TAG_LOG (16)
#define TAGE_BANK_LOG (6)
#define TAGE_BANK_SIZE (1 << TAGE_BANK_LOG)
#define TAGE_BANK_CNT (60)
#define TAGE_SIZE (TAGE_BANK_CNT * TAGE_BANK_SIZE)

static uint64_t commited_seq;

static uint64_t global_val_hist = 0;
static uint64_t global_pth_hist = 0;
static uint64_t global_tgt_hist = 0;

struct tage_entry_t {
    uint64_t last_value;
    uint64_t stride; // Sign is maintained at prediction level.
    uint64_t tag;
    uint64_t confidence;
    uint64_t useful;
    uint64_t first_occurrence;
};
static tage_entry_t tage_entries[TAGE_SIZE];

// TODO: Run experiments to check the impact.
//
#define MISPRED_THRESHOLD (128)
static int64_t last_mispred = 0;

// TODO: Run experiments to check the impact.
//
// A counter to age the entries of VTAGE by
// decreasing the useful entries by 1.
#define AGE_THRESHOLD (1024)
static int64_t age_counter = 0;

// TODO: Experiment with rolling hash by converting
// hashee into strings. As of now stealing the hash
// functions from EVES implementation.
//
// tage_index computes hash of `pc` with the
// global history to index into `i`th vtage bank.
uint64_t tage_index(int i, uint64_t pc)
{
    uint64_t hash_val = 0;
    uint64_t hist_len, aux;
    hist_len = tage_hist[i];
    aux = (pc >> i) ^ pc;

    // hash with global path history
    aux ^= ((1 << hist_len) - 1) & global_pth_hist;
    for (int j = 0; j < TAGE_HIST_LEN; j++) {
        hash_val ^= aux;
        aux ^= (aux & 15) << 16;
        aux >>= (TAGE_BANK_LOG - (TAGE_HIST_LEN - j) % (TAGE_BANK_LOG - 1));
    }

    // hash with global target history
    aux ^= ((1 << hist_len) - 1) & global_tgt_hist;
    for (int j = 0; j < TAGE_HIST_LEN; j++) {
        hash_val ^= aux;
        aux ^= (aux & 15) << 16;
        aux >>= (TAGE_BANK_LOG - (TAGE_HIST_LEN - j) % (TAGE_BANK_LOG - 1));
    }

    // hash with global value history
    aux ^= ((1 << hist_len) - 1) & global_val_hist;
    for (int j = 0; j < TAGE_HIST_LEN; j++) {
        hash_val ^= aux;
        aux ^= (aux & 15) << 16;
        aux >>= (TAGE_BANK_LOG - (TAGE_HIST_LEN - j) % (TAGE_BANK_LOG - 1));
    }

    return hash_val & (TAGE_BANK_SIZE - 1);
}

// TODO: Experiment with rolling hash by converting
// hashee into strings. As of now stealing the hash
// functions from EVES implementation.
//
// tage_tag computes hash of `pc` with the
// global history to get tag of `i`th vtage bank.
uint64_t tage_tag(int i, uint64_t pc)
{
    uint64_t hash_val = 0;
    uint64_t hist_len, aux;
    hist_len = tage_hist[i];
    aux = (pc >> (i + 5)) ^ (pc >> i) ^ pc;

    // hash with global path history
    aux ^= ((1 << hist_len) - 1) & global_pth_hist;
    for (int j = 0; j < TAGE_HIST_LEN; j++) {
        hash_val ^= aux;
        aux ^= (aux & 31) << 14;
        aux >>= (TAGE_TAG_LOG - (TAGE_HIST_LEN - j) % (TAGE_TAG_LOG - 1));
    }

    // hash with global target history
    aux ^= ((1 << hist_len) - 1) & global_tgt_hist;
    for (int j = 0; j < TAGE_HIST_LEN; j++) {
        hash_val ^= aux;
        aux ^= (aux & 31) << 14;
        aux >>= (TAGE_TAG_LOG - (TAGE_HIST_LEN - j) % (TAGE_TAG_LOG - 1));
    }

    // hash with global value history
    aux ^= ((1 << hist_len) - 1) & global_val_hist;
    for (int j = 0; j < TAGE_HIST_LEN; j++) {
        hash_val ^= aux;
        aux ^= (aux & 31) << 14;
        aux >>= (TAGE_TAG_LOG - (TAGE_HIST_LEN - j) % (TAGE_TAG_LOG - 1));
    }

    return hash_val & ((1 << TAGE_TAG_LOG) - 1);
}

#define INFLIGHT_MAX (1 << 9)
struct inflight_entry_t {
    bool predicted;
    bool prediction_result;
    bool to_update;
    uint64_t global_index[TAGE_HIST_LEN + 1];
    uint64_t global_tag[TAGE_HIST_LEN + 1];
    uint64_t hit_bank;
    uint64_t pc;
    uint64_t predicted_val;
    uint8_t instruction_type;
    uint8_t operands_cnt;
};
static inflight_entry_t inflight_entries[INFLIGHT_MAX];

#endif
