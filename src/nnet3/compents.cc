/*
* This File support nnet3 compents calculation apis
*
* author: xiangxyq
* E-mail: xiangxyq@foxmail.com
*/

#include "compents.h"
#include "../utils/asr-math.h"

void AffinePropagate(Matrix *in, Matrix *weights, BaseFloat *bias, int32 bias_len, Matrix *out)
{

	// No need for asserts as they'll happen within the matrix operations.
	CopyRowsFromVec(bias, bias_len, out);// copies bias_params_ to each row of *out
	AddMatMat(1.0f, in, kNoTrans, weights, kTrans, 1.0f, out);
}

void ReluPropagate(Matrix *in, Matrix *out)
{
	// Apply rectified linear function (x >= 0 ? 1.0 : 0.0)
	CopyFromMat(in, out);
	MatApplyFloor(out, 0.0f);
}