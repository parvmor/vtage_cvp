#include <cassert>
#include <inttypes.h>
#include <iostream>
#include <utility>

#include "mypredictor.h"
#include "../cvp_kit/cvp.h"

#define NOT_L1_MISS  (actual_latency < 12)
#define NOT_L2_MISS  (actual_latency < 60)
#define NOT_LLC_MISS (actual_latency < 150)
#define FAST_INST (actual_latency <= 3)

// beginPredictor does nothing as of now.
void beginPredictor(int argc_other, char **argv_other)
{
    return;
}

// endPredictor does nothing as of now.
void endPredictor()
{
    return;
}

// TODO: Check the hashing done here.
// Something seems to be wrong.
//
// fill_hashed_values fills the inflight_entry with
// indices and tags of the pc in all tage banks.
void fill_hashed_values(inflight_entry_t *inflight_entry)
{
    uint64_t pc = inflight_entry->pc;
    uint64_t pc_index = ((pc >> 0) ^ (pc >> 2) ^ (pc >> 5)) % TAGE_SIZE;
    uint64_t pc_bank = (pc_index >> TAGE_BANK_LOG) << TAGE_BANK_LOG;

    inflight_entry->global_index[0] = pc_index;
    inflight_entry->global_tag[0] =
        ((pc >> 0) ^ (pc >> 4) ^ (pc >> TAGE_TAG_LOG)) & ((1 << TAGE_TAG_LOG) - 1);

    for (int i = 0; i < TAGE_HIST_LEN; i++) {
        inflight_entry->global_index[i + 1] =
            (tage_index(i, pc) + (pc_bank + (i << TAGE_BANK_LOG))) % TAGE_SIZE;
        inflight_entry->global_tag[i] = tage_tag(i, pc);
    }
}

bool getPrediction(uint64_t seq_no, uint64_t pc, uint8_t piece, uint64_t& predicted_value)
{
    inflight_entry_t *inflight_entry;
    inflight_entry = inflight_entries + (seq_no & (INFLIGHT_MAX - 1));
    inflight_entry->pc = pc + piece;
    inflight_entry->predicted_val = 0xdeadbeef;
    inflight_entry->predicted = false;
    inflight_entry->hit_bank = -1;

    fill_hashed_values(inflight_entry);

    for (int i = TAGE_HIST_LEN; i >= 0; i--) {
        uint64_t tag_1 = tage_entries[inflight_entry->global_index[i]].tag;
        uint64_t tag_2 = inflight_entry->global_tag[i];
        if (tag_1 == tag_2) {
            inflight_entry->hit_bank = i;
            break;
        }
    }

    if (last_mispred >= MISPRED_THRESHOLD && inflight_entry->hit_bank >= 0) {
        uint64_t tage_idx = inflight_entry->global_index[inflight_entry->hit_bank];
        tage_entry_t *tage_entry = tage_entries + tage_idx;
        if (tage_entry->confidence >= CONF_THRESHOLD) {
            int64_t last_value = int64_t(tage_entry->last_value);
            int64_t stride = int64_t(tage_entry->stride);
            int64_t inflight_cnt = 0;
            for (uint64_t i = commited_seq + 1; i < seq_no; i++) {
                inflight_cnt += (inflight_entries[i & (INFLIGHT_MAX - 1)].pc == pc);
            }
            predicted_value = uint64_t(last_value + (inflight_cnt + 1) * stride);
            inflight_entry->predicted_val = predicted_value;
            inflight_entry->predicted = true;
        }
    }

    return inflight_entry->predicted;
}

void speculativeUpdate(uint64_t seq_no, bool eligible, uint8_t prediction_result,
                       uint64_t pc, uint64_t next_pc, InstClass insn, uint8_t piece,
                       uint64_t src1, uint64_t src2, uint64_t src3, uint64_t dst)
{
    return;
}

void updatePredictor(uint64_t seq_no, uint64_t actual_addr,
                     uint64_t actual_value, uint64_t actual_latency)
{
    return;
}
