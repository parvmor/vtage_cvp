#include <cassert>
#include <inttypes.h>
#include <iostream>
#include <cstring>
#include <utility>

#include "mypredictor.h"
#include "../cvp_kit/cvp.h"

// beginPredictor does nothing as of now.
void beginPredictor(int argc_other, char **argv_other)
{
    memset(tage_entries, 0, sizeof(tage_entry_t) * TAGE_SIZE);
    memset(inflight_entries, 0, sizeof(inflight_entry_t) * INFLIGHT_MAX);
    return;
}

// endPredictor does nothing as of now.
void endPredictor()
{
    return;
}

// TODO: Check the hashing done here.
// Something feels to be wrong.
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

// getPrediction application can be seen in `cvp.h`
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
            for (uint64_t i = committed_seq + 1; i < seq_no; i++) {
                inflight_cnt += (inflight_entries[i & (INFLIGHT_MAX - 1)].pc == pc);
            }
            predicted_value = uint64_t(last_value + (inflight_cnt + 1) * stride);
            inflight_entry->predicted_val = predicted_value;
            inflight_entry->predicted = true;
        }
    }

    return inflight_entry->predicted;
}

// speculativeUpdate application can be seen in `cvp.h`
void speculativeUpdate(uint64_t seq_no, bool eligible, uint8_t prediction_result,
                       uint64_t pc, uint64_t next_pc, InstClass insn, uint8_t piece,
                       uint64_t src1, uint64_t src2, uint64_t src3, uint64_t dst)
{
    inflight_entry_t *inflight_entry;
    inflight_entry = inflight_entries + (seq_no & (INFLIGHT_MAX - 1));
    if (eligible) {
        last_mispred += 1;
        inflight_entry->operands_cnt =
            int8_t((src1 != 0xdeadbeef) + (src2 != 0xdeadbeef) + (src3 != 0xdeadbeef));
        inflight_entry->to_update = true;
        inflight_entry->instruction_type = insn;
        inflight_entry->prediction_result = prediction_result;
        if (prediction_result == 0) {
            last_mispred = 0;
        }
    }

    bool is_branch = false;
    is_branch |= insn == condBranchInstClass;
    is_branch |= insn == uncondDirectBranchInstClass;
    is_branch |= insn == uncondIndirectBranchInstClass;
    if (is_branch) {
        // TODO: Check why this?
        // `srv_0.gz` does not shows any changes.
        // if (pc != next_pc - 4) {
        //    global_pth_hist = (global_pth_hist << 1) ^ (pc >> 2);
        //    global_tgt_hist = (global_tgt_hist << 1) ^ (next_pc >> 2);
        // }
        global_pth_hist = (global_pth_hist << 1) ^ (pc >> 2);
        global_tgt_hist = (global_tgt_hist << 1) ^ (next_pc >> 2);
    }

    return;
}

// is_instruction_critical returns *probabilistically* whether
// an instruction is critical depending on the latency and
// value of instruction. Heuristics are:
// 1) High valued instructions are more critical.
// 2) High latency instructions are more critical.
bool is_instruction_critical(inflight_entry_t *inflight_entry,
                             uint64_t actual_value, uint64_t actual_latency)
{
    bool critical = false;
    switch (inflight_entry->instruction_type) {
    case aluInstClass:
    case fpInstClass:
    case slowAluInstClass:
    case undefInstClass:
    case loadInstClass:
    case storeInstClass:
        critical = IS_INST_CRITICAL;
        break;
    case uncondIndirectBranchInstClass:
        critical = true;
        break;
    case uncondDirectBranchInstClass:
    case condBranchInstClass:
    default:
        critical = false;
        break;
    }
    return critical;
}

void update_tage_predictor(inflight_entry_t *inflight_entry,
                           uint64_t actual_value, uint64_t actual_latency)
{
    if (inflight_entry->hit_bank != -1) {
        // we had a bank hit in tage table
        int64_t hit_bank = inflight_entry->hit_bank;
        uint64_t tage_idx = inflight_entry->global_index[hit_bank];
        uint64_t tage_tag = inflight_entry->global_tag[hit_bank];
        tage_entry_t *tage_entry = tage_entries + tage_idx;

        // Ideally, the instruction should be inflight.
        // But it may happen that it was replaced.
        // As of now raise an error in such a case.
        assert(tage_entry->tag == tage_tag);

        int64_t last_value = int64_t(tage_entry->last_value);
        int64_t stride = int64_t(tage_entry->stride);
        int64_t actual_stride = int64_t(actual_value) - last_value;
        bool stride_in_limit = STRIDE_MIN <= actual_stride && actual_stride <= STRIDE_MAX;
        // Assuming updatePrediction sequence and getPrediction sequence matches
        // Cannot use inflight_entry->predicted_val since prediction is not necessary
        uint64_t predicted_val = uint64_t(last_value + stride);

        if (tage_entry->first_occurrence == 1) {
            tage_entry->first_occurrence = 0;
            if (stride_in_limit) {
                tage_entry->stride = uint64_t(actual_stride);
            } else {
                tage_entry->stride = uint64_t(STRIDE_MAX + 1);
                tage_entry->confidence = 0;
                tage_entry->useful = 0;
            }
        } else {
            if (predicted_val == actual_value) {
                // could have been a correct prediction
                if (tage_entry->confidence < CONF_MAX) {
                    if (is_instruction_critical(inflight_entry, actual_value, actual_latency)) {
                        tage_entry->confidence += 1;
                    }
                }

                if (tage_entry->useful < USE_MAX) {
                    bool flag = is_instruction_critical(inflight_entry,
                                                        actual_value,
                                                        actual_latency);
                    flag |= tage_entry->confidence == CONF_MAX;
                    if (flag) {
                        tage_entry->useful += 1;
                    }
                }

                if (tage_entry->useful < USE_MAX) {
                    if (tage_entry->confidence == CONF_MAX) {
                        tage_entry->useful += 1;
                    }
                }
            } else {
                // incorrect prediction
                if (tage_entry->confidence >= CONF_THRESHOLD) {
                    tage_entry->confidence -= CONF_THRESHOLD;
                    if (tage_entry->useful >= 1) {
                        tage_entry->useful -= 1;
                    }
                } else {
                    tage_entry->confidence = 0;
                    tage_entry->useful = 0;
                }
                tage_entry->first_occurrence = 1;
            }
        }

    }

    if (inflight_entry->prediction_result != 1) {
        if (is_instruction_critical(inflight_entry, actual_value, actual_latency)) {
            int64_t hit_bank = inflight_entry->hit_bank + 1;
            for (int64_t i = hit_bank; i <= TAGE_HIST_LEN; i++) {
                uint64_t tage_idx = inflight_entry->global_index[i];
                tage_entry_t *tage_entry = tage_entries + tage_idx;
                bool flag = tage_entry->confidence < CONF_THRESHOLD;
                flag |= tage_entry->confidence < (rng() & CONF_MAX);
                flag &= tage_entry->useful == 0;
                if (flag) {
                    tage_entry->last_value = actual_value;
                    tage_entry->stride = 0;
                    tage_entry->tag = inflight_entry->global_tag[i];
                    tage_entry->first_occurrence = 1;
                    tage_entry->useful = 0;
                    tage_entry->confidence = 1;
                    age_counter -= 3;
                    break;
                }
                age_counter += 1;
            }
        }
        if (age_counter >= AGE_THRESHOLD) {
            for (int i = 0; i < TAGE_SIZE; i++) {
                if (tage_entries[i].useful > 0) {
                    tage_entries[i].useful -= 1;
                }
            }
            age_counter = 0;
        }
    }
    return;
}

// updatePredictor application can be seen in `cvp.h`
void updatePredictor(uint64_t seq_no, uint64_t actual_addr,
                     uint64_t actual_value, uint64_t actual_latency)
{
    inflight_entry_t *inflight_entry;
    inflight_entry = inflight_entries + (seq_no & (INFLIGHT_MAX - 1));

    if (inflight_entry->to_update) {
        update_tage_predictor(inflight_entry, actual_value, actual_latency);
        inflight_entry->to_update = false;
    }

    global_val_hist = (global_val_hist << 1) ^ (actual_value >> 2);

    // Change the committed sequence number.
    committed_seq = seq_no;

    return;
}
