/*
* This File chain neural network compute
*
* author: xiangxyq
* E-mail: xiangxyq@foxmail.com
*/

#include <vector>
#include <cstdio>
#include <cassert>
#include <cstring>
#include <iostream>
#include "chain-compute.h"
#include "../configs/nnet3-conf.h"
#include "compents.h"
#include "../models/final-mdl.h"
#include "../utils/asr-math.h"
#include "../configs/feat-conf.h"

typedef struct
{
	int32 dim;
	int32 data[INDEXES_MAX_DIM];
} CuArray;

static int32 program_counter_ = 0;
static std::vector<int32> pending_commands_;
static const int32 id2pdf_id_[ID2PDF_ID_LEN] = ID2PDF_ID_VAL;
static BaseFloat *weights_params_, *bias_params_;
static int32 weights_rows_, weights_cols_, weights_stride_, bias_dim_;
static Command cmds[NNET_COMMANDS_LEN] = NNET_COMMANDS;
static SubMatrixInfo submatrices[SUBMATRICES_LEN] = SUBMATRICES;
static MatrixInfo matrices[MATRICES_LEN] = MATRICES;
static Matrix matrices_[MATRICES_LEN];
static std::string node_names_[NODE_NUMS] = NODE_NAMES;
static CuArray indexes_[INDEXES_DIM] = INDEXES;
static CompentType commpent_[COMPENTS_LEN] = COMPENTS;
static Matrix output, current_log_post_;
static int32 current_log_post_subsampled_offset_ = -1;
static int32 num_chunks_computed_ = 0;

static void CheckNoPendingIo();
static void ExecuteCommand();
static bool IsWholeMatrix(int32 submatrix_index);
static int32 GetNodeIndex(const std::string &node_name);
static int32 GetIoMatrixIndex(const std::string &node_name, bool is_output);
static void GetSubMatrix(int32 submatrix_index, Matrix *m);
static void CopyRows(Matrix *dst, Matrix *src, const int32 *indices);
static void ComputerRun();
static void GetOutputDestructive(const std::string &node_name, Matrix *output);
static void AcceptInput(const std::string &node_name, Matrix *input);
static void AdvanceChunk(const std::vector<BaseFloat* > &feats);
static void EnsureFrameIsComputed(const std::vector<BaseFloat* > &feats, int32 subsampled_frame);
static int32 TransitionIdToPdfFast(int32 trans_id);


void CheckNoPendingIo()
{
	while (program_counter_ < static_cast<int32>(NNET_COMMANDS_LEN) &&
		(cmds[program_counter_].command_type == kAcceptInput ||
			cmds[program_counter_].command_type == kProvideOutput))
	{
		pending_commands_.push_back(program_counter_);
		program_counter_++;
	}

	for (int32 i = 0; i < pending_commands_.size(); ++i)
	{
		// the order here doesn't really matter; we go from back to front
		// as it's more efficient, not that efficiency really matters here.
		int32 command = pending_commands_[i];
		if (-1 == command)
		{
			std::cout << "Error pending commands id, need to check!" << std::endl;
			exit(1);
		}
		if (cmds[command].command_type == kAcceptInput) 
		{
			// we can't ignore if we needed input from the user that hasn't been
			// provided.
			int32 node = cmds[command].arg2;
			std::cout << "Cannot run computation -- we did not get input for node " << node << std::endl;
			exit(1);
		}
	}
	pending_commands_.clear();
}


void ExecuteCommand()
{
	const Command c = cmds[program_counter_];
	int32 m1, m2;

	switch (c.command_type)
	{
		case kAllocMatrix:
			m1 = submatrices[c.arg1].matrix_index;
			ResizeMatrix(matrices[m1].num_rows, matrices[m1].num_cols, kDefaultStride, &(matrices_[m1]));
			break;

		case kDeallocMatrix:
			m1 = submatrices[c.arg1].matrix_index;
			ResizeMatrix(0, 0, kDefaultStride, &(matrices_[m1]));
			break;

		case kSwapMatrix:
			m1 = submatrices[c.arg1].matrix_index;
			m2 = submatrices[c.arg2].matrix_index;
			SwapMatrix(&(matrices_[m1]), &(matrices_[m2]));
			break;

		case kPropagate: 
		{
			Matrix input, output;
			GetSubMatrix(c.arg3, &input);
			GetSubMatrix(c.arg4, &output);

			if (AffineComponent == commpent_[c.arg1] || NaturalGradientAffineComponent == commpent_[c.arg1])
			{
				Matrix tmp_weights;
				weights_params_ = weights[c.arg1];
				weights_rows_ = weights_rows[c.arg1];
				weights_cols_ = weights_cols[c.arg1];
				weights_stride_ = weights_stride[c.arg1];
				bias_params_ = bias[c.arg1];
				bias_dim_ = bias_len[c.arg1];

				tmp_weights.data = weights_params_;
				tmp_weights.mInfo.num_rows = weights_rows_;
				tmp_weights.mInfo.num_cols = weights_cols_;
				tmp_weights.mInfo.stride = weights_stride_;
				assert(NULL != weights_params_ || weights_rows_ != 0 || weights_cols_ != 0 || bias_params_ != NULL || bias_dim_ != 0);

				AffinePropagate(&input, &tmp_weights, bias_params_, bias_dim_, &output);
			}
			else if (RectifiedLinearComponent == commpent_[c.arg1])
				ReluPropagate(&input, &output);
			else
				std::cout << "Other net component, maybe this is error!" << std::endl;
			
			break;
		}

		case kMatrixCopy: 
		{
			Matrix dest, src;
			GetSubMatrix(c.arg1, &dest);
			GetSubMatrix(c.arg2, &src);
			CopyFromMat(&src, &dest);

			break;
		}

		case kCopyRows: 
		{
			Matrix dest, src;
			GetSubMatrix(c.arg1, &dest);
			GetSubMatrix(c.arg2, &src);
			const CuArray indexes = indexes_[c.arg3];
			CopyRows(&dest, &src, indexes.data);

			break;
		}

		case kSetConst:
		case kBackprop:
		case kBackpropNoModelUpdate:
		case kMatrixAdd:
		case kAddRows:
		case kCopyRowsMulti:
		case kCopyToRowsMulti:
		case kAddRowsMulti:
		case kAddToRowsMulti:
		case kAddRowRanges:
		case kCompressMatrix: case kDecompressMatrix:
		case kNoOperation: case kNoOperationPermanent: case kNoOperationMarker:
		case kNoOperationLabel:
			break;
		case kGotoLabel:
			assert(cmds[c.arg1].command_type == kNoOperationLabel);
			program_counter_ = c.arg1;
			break;
		default:
			std::cout << "Invalid command in computation" << std::endl;
			exit(1);
	}
}

void ComputerRun()
{
	if (program_counter_ >= NNET_COMMANDS_LEN)
	{
		std::cout << "Running computation that has finished: program-counter=" << program_counter_ << std::endl;
		exit(1);
	}
	CheckNoPendingIo();
	for (; program_counter_ < NNET_COMMANDS_LEN; ++program_counter_)
	{
		if (cmds[program_counter_].command_type == kAcceptInput ||
			cmds[program_counter_].command_type == kProvideOutput) {
			// We have hit a part of the computation that requires user
			// interaction, e.g. the end of the forward or backward phase.
			break;
		}
		ExecuteCommand();
	}
}

void AcceptInput(const std::string &node_name, Matrix *input)
{
	bool is_output = false;
	int32 matrix_index = GetIoMatrixIndex(node_name, is_output);
	const MatrixInfo matrix_info = matrices[matrix_index];

	if (input->mInfo.num_rows != matrix_info.num_rows) 
	{
		std::cout << "Num-rows mismatch for input " << node_name << " : " << matrix_info.num_rows \
			<< " in computation-request, " << input->mInfo.num_rows << " provided." << std::endl;
		exit(0);
	}
	if (input->mInfo.num_cols != matrix_info.num_cols) 
	{
		std::cout << "Num-cols mismatch for input " << node_name << " : " << matrix_info.num_cols \
				  << " in computation-request, " << input->mInfo.num_cols << " provided." << std::endl;
		exit(0);
	}
	if (input->mInfo.stride == input->mInfo.num_cols) 
	{
		SwapMatrix(&(matrices_[matrix_index]), input);
		FreeMatrix(input);
	}
	else 
	{
		ResizeMatrix(matrix_info.num_rows, matrix_info.num_cols, kStrideEqualNumCols, &(matrices_[matrix_index]));
		memcpy((BaseFloat*)(matrices_[matrix_index].data), (BaseFloat*)(input), sizeof(BaseFloat) * (matrix_info.num_rows * matrix_info.num_cols));
		ResizeMatrix(0, 0, kDefaultStride, input);
	}
}

int32 GetNodeIndex(const std::string &node_name)
{
	for (int32 i = 0; i < NODE_NUMS; ++i)
		if (node_names_[i] == node_name)
			return i;
	return -1;
}

int32 GetIoMatrixIndex(const std::string &node_name, bool is_output)
{
	int32 node_index = GetNodeIndex(node_name);
	if (node_index == -1)
	{
		std::cout << "network has no node named " << node_name.c_str() << std::endl;
		exit(1);
	}
	// first make sure all the I/O commands that we immediately expect, are listed
	// in 'pending_commands_'.
	while (program_counter_ < NNET_COMMANDS_LEN &&
		((cmds[program_counter_].command_type == kAcceptInput ||
			cmds[program_counter_].command_type == kProvideOutput ||
			cmds[program_counter_].command_type == kNoOperationMarker))) {
		if (cmds[program_counter_].command_type != kNoOperationMarker)
		{
			pending_commands_.push_back(program_counter_);
	    }
		program_counter_++;
	}
	
	for (int32 i = 0; i < pending_commands_.size(); ++i) 
	{
		int32 command = pending_commands_[i];
		bool this_command_is_output = (cmds[command].command_type == kProvideOutput);
		int32 this_submatrix_index = cmds[command].arg1,
			this_node_index = cmds[command].arg2;
		if (this_command_is_output == is_output && node_index == this_node_index) 
		{
			if (!is_output) 
			{
				pending_commands_.erase(pending_commands_.begin() + i);
				//pending_commands_.erase(pending_commands_.begin() + i);
				// don't erase the command for outputs, as that would prevent things
				// from being output twice, which is an unnecessary restriction.
			}
			if (!(IsWholeMatrix(this_submatrix_index)))
			{
				std::cout << "Getting input or output that is not a whole matrix (probably some optimization code needs to be changed)" << std::endl;
				exit(1);
			}

			return submatrices[this_submatrix_index].matrix_index;
		}
	}
	// if you get the following error it will likely be a bug in the calling code,
	// or possibly due to giving the wrong egs.
	assert("it is not expected at this point in the computation!");
	return 0;  // Suppress compiler warnings; this line will never be reached.
}

bool IsWholeMatrix(int32 submatrix_index)
{
	assert(submatrix_index > 0 && submatrix_index < SUBMATRICES_LEN);

	const SubMatrixInfo submat_info = submatrices[submatrix_index];
	const MatrixInfo mat_info = matrices[submat_info.matrix_index];

	return submat_info.row_offset == 0 && submat_info.col_offset == 0 &&
		submat_info.num_rows == mat_info.num_rows &&
		submat_info.num_cols == mat_info.num_cols;
}

void GetOutputDestructive(const std::string &node_name, Matrix *output)
{
	bool is_output = true;
	int32 matrix_index = GetIoMatrixIndex(node_name, is_output);
	assert(matrices_[matrix_index].mInfo.num_rows != 0);

	SwapMatrix(&(matrices_[matrix_index]), output);
	ResizeMatrix(0, 0, kDefaultStride, &(matrices_[matrix_index]));
}


void GetSubMatrix(int32 submatrix_index, Matrix *m)
{
	assert(submatrix_index < SUBMATRICES_LEN);

	const SubMatrixInfo info = submatrices[submatrix_index];
	const Matrix mat = matrices_[info.matrix_index];
	m->data = mat.data + info.col_offset + info.row_offset * (mat.mInfo.stride);
	m->mInfo.num_cols = info.num_cols;
	m->mInfo.num_rows = info.num_rows;
	m->mInfo.stride = mat.mInfo.stride;
}

void CopyRows(Matrix *dst, Matrix *src, const int32 *indices)
{
	assert(dst->mInfo.num_cols == src->mInfo.num_cols);

	int32 num_rows = dst->mInfo.num_rows, num_cols = dst->mInfo.num_cols,
		this_stride = dst->mInfo.stride;
	BaseFloat *this_data = dst->data;

	for (int32 r = 0; r < num_rows; ++r, this_data += this_stride)
	{
		int32 index = indices[r];
		if (index < 0)
			memset(this_data, 0, sizeof(BaseFloat) * dst->mInfo.num_cols);
		else
			memcpy(this_data, src->data + (index * src->mInfo.stride), num_cols * sizeof(BaseFloat));
	}
}

void AdvanceChunk(const std::vector<BaseFloat* > &feats)
{
	// Prepare the input data for the next chunk of features.
	// note: 'end' means one past the last.
	int32 begin_input_frame, end_input_frame;
	if (num_chunks_computed_ == 0) 
	{
		begin_input_frame = -FARMES_LEFT_CONTEXT;
		// note: end is last plus one.
		end_input_frame = FRAMES_PER_CHUNK + FRAMES_RIGHT_CONTEXT;
	}
	else 
	{
		// note: begin_input_frame will be the same as the previous end_input_frame.
		// you can verify this directly if num_chunks_computed_ == 0, and then by
		// induction.
		begin_input_frame = num_chunks_computed_ * FRAMES_PER_CHUNK + FRAMES_RIGHT_CONTEXT;
		end_input_frame = begin_input_frame + FRAMES_PER_CHUNK;
	}
	int32 num_feature_frames_ready = NumFramesReady(feats);
	int32 is_finished = IsLastFrame(feats, num_feature_frames_ready - 1);

	if (end_input_frame > num_feature_frames_ready && !is_finished)
	{
		// we shouldn't be attempting to read past the end of the available features
		// until we have reached the end of the input (i.e. the end-user called
		// InputFinished(), announcing that there is no more waveform; at this point
		// we pad as needed with copies of the last frame, to flush out the last of
		// the output.
		// If the following error happens, it likely indicates a bug in this
		// decodable code somewhere (although it could possibly indicate the
		// user asking for a frame that was not ready, which would be a misuse
		// of this class.. it can be figured out from gdb as in either case it
		// would be a bug in the code.
		assert("Attempt to access frame past the end of the available input");
	}

	Matrix feats_chunk;
	{ // this block sets 'feats_chunk'.
		Matrix this_feats;
		ResizeMatrix(end_input_frame - begin_input_frame, NUM_CEPS, kDefaultStride, &this_feats);
		for (int32 i = begin_input_frame; i < end_input_frame; ++i) 
		{
			int32 input_frame = i;
			if (input_frame < 0) input_frame = 0;
			if (input_frame >= num_feature_frames_ready)
				input_frame = num_feature_frames_ready - 1;
			BaseFloat *this_row = GetFrame(feats, input_frame);
			if (NULL != this_row)
				memcpy(this_feats.data + (i - begin_input_frame) * this_feats.mInfo.stride, this_row, sizeof(BaseFloat) * NUM_CEPS);
		}
		SwapMatrix(&feats_chunk, &this_feats);
	}
	AcceptInput("input", &feats_chunk);

	ComputerRun();

	{
		// Note: it's possible in theory that if you had weird recurrence that went
		// directly from the output, the call to GetOutputDestructive() would cause
		// a crash on the next chunk.  If that happens, GetOutput() should be used
		// instead of GetOutputDestructive().  But we don't anticipate this will
		// happen in practice.
		GetOutputDestructive("output", &output);
		// apply the acoustic scale
		//output.Scale(acoustic_scale);  //不需要sacle因为acoustic_scale = 1.0

		ResizeMatrix(0, 0, kDefaultStride, &current_log_post_);
		SwapMatrix(&output, &current_log_post_);
		FreeMatrix(&output);	
	}
	assert(current_log_post_.mInfo.num_rows == FRAMES_PER_CHUNK / FRAME_SUBSAMPLING_FACTOR &&
		current_log_post_.mInfo.num_cols == OUTPUT_DIM);
	num_chunks_computed_++;
	current_log_post_subsampled_offset_ = (num_chunks_computed_ - 1) * (FRAMES_PER_CHUNK / FRAME_SUBSAMPLING_FACTOR);
}

static void EnsureFrameIsComputed(const std::vector<BaseFloat* > &feats, int32 subsampled_frame)
{
	assert(subsampled_frame >= current_log_post_subsampled_offset_);

	while (subsampled_frame >= current_log_post_subsampled_offset_ +
		current_log_post_.mInfo.num_rows)
	{
		AdvanceChunk(feats);
	}

}

int32 TransitionIdToPdfFast(int32 trans_id)
{
	// Note: it's a little dangerous to assert this only in paranoid mode.
	// However, this function is called in the inner loop of decoders and
	// the assertion likely takes a significant amount of time.  We make
	// sure that past the end of thd id2pdf_id_ array there are big
	// numbers, which will make the calling code more likely to segfault
	// (rather than silently die) if this is called for out-of-range values.
	if (trans_id >= ID2PDF_ID_LEN)
	{
		std::cout << "Likely graph/model mismatch (graph built from wrong model?)" << std::endl;
		exit(1);
	}

	return id2pdf_id_[trans_id];
}

void InitChainCompute()
{
	program_counter_ = 0;
	current_log_post_subsampled_offset_ = -1;
	num_chunks_computed_ = 0;

	pending_commands_.clear();
	FreeChianCompute();
}

int32 NumFramesDecodeReady(const std::vector<BaseFloat* > &feats)
{
	// note: the ivector_features_ may have 2 or 3 fewer frames ready than
	// input_features_, but we don't wait for them; we just use the most recent
	// iVector we can.
	int32 features_ready = NumFramesReady(feats);
	if (features_ready == 0)
		return 0;
	bool input_finished = IsLastFrame(feats, features_ready - 1);

	int32 sf = FRAME_SUBSAMPLING_FACTOR;

	if (input_finished) 
	{
		// if the input has finished,... we'll pad with duplicates of the last frame
		// as needed to get the required right context.
		return (features_ready + sf - 1) / sf;
	}
	else 
	{
		// note: info_.right_context_ includes both the model context and any
		// extra_right_context_ (but this
		int32 non_subsampled_output_frames_ready = ((features_ready - FRAMES_RIGHT_CONTEXT) > 0) ? (features_ready - FRAMES_RIGHT_CONTEXT) : 0;
		int32 num_chunks_ready = non_subsampled_output_frames_ready / FRAMES_PER_CHUNK;
		// note: the division by the frame subsampling factor 'sf' below
		// doesn't need any attention to rounding because info_.frames_per_chunk
		// is always a multiple of 'sf' (see 'frames_per_chunk = GetChunksize..."
		// in decodable-simple-looped.cc).
		return num_chunks_ready * FRAMES_PER_CHUNK / sf;
	}
}

BaseFloat LogLikelihood(const std::vector<BaseFloat* > &feats, int32 subsampled_frame, int32 index)
{
	
	EnsureFrameIsComputed(feats, subsampled_frame);
	return GetMatPoVal(&current_log_post_,
		subsampled_frame - current_log_post_subsampled_offset_, TransitionIdToPdfFast(index));
}

void FreeChianCompute()
{
	FreeMatrix(&output);
	FreeMatrix(&current_log_post_);
	for (int32 i = 0; i < MATRICES_LEN; ++i)
		FreeMatrix(&matrices_[i]);
}
