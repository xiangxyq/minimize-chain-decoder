/*
* This File support matrix apis
*
* author: xiangxyq
* E-mail: xiangxyq@foxmail.com
*/

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <cstring>
#include "matrix.h"

void ResizeMatrix(int32 rows, int32 cols, MatrixStrideType stride_type, Matrix *matrix)
{
	if (rows * cols == 0)
	{
		assert(rows == 0 && cols == 0);
		FreeMatrix(matrix);
		return;
	}
	assert(rows > 0 && cols > 0);

	int32 skip, stride, size;
	// compute the size of skip and real cols
	skip = ((16 / sizeof(BaseFloat)) - cols % (16 / sizeof(BaseFloat))) % (16 / sizeof(BaseFloat));
	stride = cols + skip;
	size = (int32)(rows) * (int32)(stride) * sizeof(BaseFloat);

	matrix->mInfo.num_rows = rows;
	matrix->mInfo.num_cols = cols;
	matrix->mInfo.stride = (stride_type == kDefaultStride ? stride : cols);
	matrix->data = (BaseFloat *)malloc(size);
	memset((BaseFloat *)matrix->data, 0, size);
}

void SwapMatrix(Matrix *m1, Matrix *m2)
{
	Matrix tmp;
	tmp.mInfo.num_rows = m1->mInfo.num_rows;
	tmp.mInfo.num_cols = m1->mInfo.num_cols;
	tmp.mInfo.stride = m1->mInfo.stride;
	tmp.data = m1->data;

	m1->mInfo.num_rows = m2->mInfo.num_rows;
	m1->mInfo.num_cols = m2->mInfo.num_cols;
	m1->mInfo.stride = m2->mInfo.stride;
	m1->data = m2->data;

	m2->mInfo.num_rows = tmp.mInfo.num_rows;
	m2->mInfo.num_cols = tmp.mInfo.num_cols;
	m2->mInfo.stride = tmp.mInfo.stride;
	m2->data = tmp.data;
}

void CopyFromVec(BaseFloat *dest, const BaseFloat *src, int32 num_cols)
{
	if (dest != src)
		memcpy(dest, src, num_cols * sizeof(BaseFloat));
}

void CopyFromMat(const Matrix *src, Matrix *dst)
{
	assert(dst->mInfo.num_rows == src->mInfo.num_rows && dst->mInfo.num_cols == src->mInfo.num_cols);

	for (int32 i = 0; i < dst->mInfo.num_rows; ++i)
		CopyFromVec(dst->data + (i * dst->mInfo.stride), src->data + (i * src->mInfo.stride), dst->mInfo.num_cols);
}

void CopyRowsFromVec(BaseFloat *rv, const int32 rv_len, Matrix *out)
{
	if (rv_len == out->mInfo.num_cols)
	{
		for (int32 r = 0; r < out->mInfo.num_rows; ++r)
			memcpy(out->data + (r * out->mInfo.stride), rv, sizeof(BaseFloat) * out->mInfo.num_cols);
	}
	else if (rv_len == out->mInfo.num_rows * out->mInfo.num_cols)
	{
		if (out->mInfo.stride == out->mInfo.num_cols)
			memcpy(out->data, rv, sizeof(BaseFloat) * out->mInfo.num_rows * out->mInfo.num_cols);
		else
		{
			for (int32 r = 0; r < out->mInfo.num_rows; ++r) 
            {
				float *row_data = out->data + r * out->mInfo.stride;
				for (int32 c = 0; c < out->mInfo.num_cols; ++c) 
					row_data[c] = rv[c];
				rv += out->mInfo.num_cols;
			}
		}
	}
	else 
	{
		std::cout << "Wrong sized arguments" << std::endl;
		exit(1);
	}
}

void SetMatPosVal(Matrix *m, int32 r, int32 c, BaseFloat value)
{
	assert(r <= m->mInfo.num_rows && c <= m->mInfo.num_cols);

	*(m->data + r * m->mInfo.stride + c) = value;
}

BaseFloat GetMatPoVal(Matrix *m, int32 r, int32 c)
{
	assert(r <= m->mInfo.num_rows && c <= m->mInfo.num_cols);

	return *(m->data + r * m->mInfo.stride + c);
}

void FreeMatrix(Matrix *matrix)
{
	free(matrix->data);
	matrix->mInfo.num_rows = 0;
	matrix->mInfo.num_cols = 0;
	matrix->mInfo.stride = 0;
	matrix->data = NULL;
}
