/*
* This File support matrix apis
*
* author: xiangxyq
* E-mail: xiangxyq@foxmail.com
*/

#ifndef MATRIX_H
#define MATRIX_H

#include "asr-types.h"

typedef enum
{
	kDefaultStride,
	kStrideEqualNumCols,
} MatrixStrideType;

typedef struct
{
	int32 num_rows;
	int32 num_cols;
	int32 stride;
} MatrixInfo;

typedef struct {
	MatrixInfo mInfo;
	BaseFloat *data;
} Matrix;


extern void ResizeMatrix(int32 rows, int32 cols, MatrixStrideType stride_type, Matrix *matrix);

extern void SwapMatrix(Matrix *m1, Matrix *m2);

extern void CopyFromVec(BaseFloat *dest, const BaseFloat *src, int32 num_cols);

extern void CopyFromMat(const Matrix *src, Matrix *dst);

extern void CopyRowsFromVec(BaseFloat *rv, const int32 rv_len, Matrix *out);

extern void SetMatPosVal(Matrix *m, int32 r, int32 c, BaseFloat value);

extern BaseFloat GetMatPoVal(Matrix *m, int32 r, int32 c);

extern void FreeMatrix(Matrix *matrix);


#endif // MATRIX_H