/*
* This File chain neural network compute
*
* author: xiangxyq
* E-mail: xiangxyq@foxmail.com
*/

#ifndef CHAIN_COMPUTE_H
#define CHAIN_COMPUTE_H

#include "../utils/asr-types.h"
#include "../utils/matrix.h"
#include "../feat/compute-feat.h"

typedef enum 
{
	kAllocMatrix, kDeallocMatrix, kSwapMatrix, kSetConst,
	kPropagate, kBackprop, kBackpropNoModelUpdate,
	kMatrixCopy, kMatrixAdd, kCopyRows, kAddRows,
	kCopyRowsMulti, kCopyToRowsMulti, kAddRowsMulti, kAddToRowsMulti,
	kAddRowRanges, kCompressMatrix, kDecompressMatrix,
	kAcceptInput, kProvideOutput,
	kNoOperation, kNoOperationPermanent, kNoOperationMarker, kNoOperationLabel,
	kGotoLabel
} CommandType;

typedef enum
{
	RectifiedLinearComponent,
	NaturalGradientAffineComponent,
	AffineComponent,
	LogSoftmaxComponent
} CompentType;

typedef struct
{
	int32 matrix_index;  // index into "matrices": the underlying matrix.
	int32 num_rows;
	int32 num_cols;
	int32 row_offset;
	int32 col_offset;
} SubMatrixInfo;

typedef struct
{
	CommandType command_type;
	BaseFloat alpha;
	int32 arg1;
	int32 arg2;
	int32 arg3;
	int32 arg4;
	int32 arg5;
	int32 arg6;
	int32 arg7;
} Command;


extern void InitChainCompute();

extern int32 NumFramesDecodeReady(const std::vector<BaseFloat* > &feats);

extern BaseFloat LogLikelihood(const std::vector<BaseFloat* > &feats, int32 subsampled_frame, int32 index);

extern void FreeChianCompute();


#endif // CHAIN_COMPUTE_H
