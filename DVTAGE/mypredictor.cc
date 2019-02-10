#include <cassert>
#include <inttypes.h>
#include <iostream>
#include <utility>

#include "mypredictor.h"
#include "../cvp_kit/cvp.h"

void beginPredictor(int argc_other, char **argv_other)
{
    return;
}

void endPredictor()
{
    return;
}

bool getPrediction(uint64_t seq_no, uint64_t pc, uint8_t piece, uint64_t& predicted_value)
{
    return false;
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
