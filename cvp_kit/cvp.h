
////////////////////////////////////////////////////////////////////////////////
//
// Interface for the 1st Championship on Value Prediction (CVP-1).
// Conceived by Arthur Perais, Rami Al Sheikh, and Eric Rotenberg.
// Authored by Eric Rotenberg.
// File first created: Feb. 7, 2018.
//
// (C) 2018 Eric Rotenberg
//
////////////////////////////////////////////////////////////////////////////////

// Instruction type.
// Example use:
//
// InstClass insn;
// ...
// if (insn == InstClass::aluInstClass) ...

enum InstClass : uint8_t
{
    aluInstClass = 0,
    loadInstClass = 1,
    storeInstClass = 2,
    condBranchInstClass = 3,
    uncondDirectBranchInstClass = 4,
    uncondIndirectBranchInstClass = 5,
    fpInstClass = 6,
    slowAluInstClass = 7,
    undefInstClass = 8
};


//
// getPrediction()
//
// Return value:
// "true" if microarch. simulator should speculate based on the prediction for this instruction.
// "false" if it should not speculate based on the prediction for this instruction.
// This allows contestants to decide between the potential speedup of speculation vs.
// the potential penalty of a squash from ROB-head due to misspeculation.
//
// Input arguments:
// 1. sequence number: the dynamic micro-instruction number
// 2. program counter (pc) of instruction
// 3. piece: Some instructions operate on values that are wider than 64 bits.
//           These are split into multiple 64-bit pieces and getPrediction() is called
//           for pieces 0, 1, ..., n, consecutively, for the same instruction.
//           Note that sequence number is incremented whether for an instruction or piece
//           of an instruction (hence sequence number is dynamic micro-instruction number).
//
// Output argument:
// 1. predicted value, whether or not the simulator is directed to speculate
//
extern
bool getPrediction(uint64_t seq_no, uint64_t pc, uint8_t piece, uint64_t& predicted_value);

//
// speculativeUpdate()
//
// This function is called immediately after getPrediction() for the just-predicted instruction.
//
// A key argument to this function is "prediction_result".
// * If getPrediction() instructed the simulator to speculate, then "prediction_result" will reveal whether or not
//   the predicted value is correct, immediately after getPrediction().
//   In other words, if the contestant took the risk of speculating, he/she has the privilege of knowing the outcome of
//   speculation immediately.  This privilege is justified whether or not the prediction was correct.  If the prediction
//   was incorrect, all instructions prior to the mispredicted one will be retired before the next call to getPrediction(),
//   including draining the window of all pending calls to updatePredictor().  Effectively, architectural state becomes
//   visible after a value misprediction and before the next value prediction.
//   On the other hand, if the prediction was correct, exposing it as correct supports a speculative-update policy
//   without the need for value predictor fix-ups, for all contestants.
// * On the other hand, if getPrediction() instructed the simulator to NOT speculate, then the contestant forfeits the
//   privilege described above.  "prediction_result" will signify that the outcome of the prediction will not be
//   revealed immediately.  The contestant must wait until the non-speculative updatePredictor() function to see
//   whether or not the prediction was correct.
//
extern
void speculativeUpdate(uint64_t seq_no,    		// dynamic micro-instruction # (starts at 0 and increments indefinitely)
                       bool eligible,			// true: instruction is eligible for value prediction. false: not eligible.
        		       uint8_t prediction_result,	// 0: incorrect, 1: correct, 2: unknown (not revealed)
        		       // Note: can assemble local and global branch history using pc, next_pc, and insn.
        		       uint64_t pc,
        		       uint64_t next_pc,
        		       InstClass insn,
        		       uint8_t piece,
        		       // Note: up to 3 logical source register specifiers, up to 1 logical destination register specifier.
        		       // 0xdeadbeef means that logical register does not exist.
        		       // May use this information to reconstruct architectural register file state (using log. reg. and value at updatePredictor()).
        		       uint64_t src1,
        		       uint64_t src2,
        		       uint64_t src3,
        		       uint64_t dst);

//
// updatePredictor()
//
// This is called for a micro-instruction when it is retired.
//
// Generally there is a delay between the getPrediction()/speculativeUpdate() calls and the corresponding updatePredictor() call,
// which is the delay between fetch and retire.  This delay manifests to contestants as multiple unrelated
// getPrediction(y,z,...)/speculativeUpdate(y,z,...) calls between an instruction x's getPrediction(x)/speculativeUpdate(x) calls and
// its updatePredictor(x) call.  The delay goes away when there is a value misprediction because of the complete-squash recovery model.
// After a value misprediction, the window is drained of all pending updatePredictor() calls before the next call to getPrediction()/speculativeUpdate().
//
extern
void updatePredictor(uint64_t seq_no,		// dynamic micro-instruction #
        		     uint64_t actual_addr,	// load or store address (0xdeadbeef if not a load or store instruction)
        		     uint64_t actual_value,	// value of destination register (0xdeadbeef if instr. is not eligible for value prediction)
        		     uint64_t actual_latency);	// actual execution latency of instruction

//
// beginPredictor()
//
// This function is called by the simulator before the start of simulation.
// It can be used for arbitrary initialization steps for the contestant's code.
// Note that the contestant may also parse additional contestant-specific command-line arguments that come
// after the final required argument, i.e., after the .gz trace file.
//
extern
void beginPredictor(int argc_other, char **argv_other);

//
// endPredictor()
//
// This function is called by the simulator at the end of simulation.
// It can be used by the contestant to print out other contestant-specific measurements.
//
extern
void endPredictor();
