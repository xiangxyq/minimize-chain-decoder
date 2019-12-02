/*
* math apis
*
* author: xiangxyq
* E-mail: xiangxyq@foxmail.com
*/

#include <cmath>
#include <cassert>
#include <cstring>
#include "asr-math.h"
#ifdef USE_OPENBLAS
#include "cblas.h"
#endif


// ###################### Below math use openblas interfaces Start ######################
BaseFloat VecVec(const BaseFloat *a, const BaseFloat *b, int32 len)
{
#ifdef USE_OPENBLAS
	return cblas_sdot(len, a, 1, b, 1);
#else
	BaseFloat sum = 0.0f;

	for (int i = 0; i < len; ++i)
		sum += a[i] * b[i];
	
	return sum;
#endif // USE_OPENBLAS
}

void AddMatVec(const BaseFloat alpha, const BaseFloat *M, int32 M_rows, int32 M_cols, int32 stride,
	MatrixTransposeType trans, const BaseFloat *v, int32 v_len, const BaseFloat beta, BaseFloat *result, int32 result_len)
{
	assert((trans == kNoTrans && M_cols == v_len && M_rows == result_len)
		   || (trans == kTrans && M_rows == v_len && M_cols == result_len));
#ifdef USE_OPENBLAS
	cblas_sgemv(CblasRowMajor, (CBLAS_TRANSPOSE)trans, M_rows,
		                              M_cols, alpha, M, stride, v, 1, beta, result, 1);
#else
	for(int i = 0; i < M_rows; ++i)
		result[i] = VecVec(M + i * M_cols, v, v_len);
#endif // USE_OPENBLAS
}

void AddMatMat(const BaseFloat alpha, const Matrix *A, MatrixTransposeType transA,
					const Matrix *B, MatrixTransposeType transB, const BaseFloat beta, Matrix *out)
{
	assert((transA == kNoTrans && transB == kNoTrans && A->mInfo.num_cols == B->mInfo.num_rows && A->mInfo.num_rows == out->mInfo.num_rows && B->mInfo.num_cols == out->mInfo.num_cols)
		|| (transA == kTrans && transB == kNoTrans && A->mInfo.num_rows == B->mInfo.num_rows && A->mInfo.num_cols == out->mInfo.num_rows && B->mInfo.num_cols == out->mInfo.num_cols)
		|| (transA == kNoTrans && transB == kTrans && A->mInfo.num_cols == B->mInfo.num_cols && A->mInfo.num_rows == out->mInfo.num_rows && B->mInfo.num_rows == out->mInfo.num_cols)
		|| (transA == kTrans && transB == kTrans && A->mInfo.num_rows == B->mInfo.num_cols && A->mInfo.num_cols == out->mInfo.num_rows && B->mInfo.num_rows == out->mInfo.num_cols));
	assert(A != out && B != out);

	if (out->mInfo.num_rows == 0) return;
#ifdef USE_OPENBLAS
	cblas_sgemm(CblasRowMajor, (CBLAS_TRANSPOSE)transA,
		(CBLAS_TRANSPOSE)transB,
		out->mInfo.num_rows, out->mInfo.num_cols, transA == kNoTrans ? A->mInfo.num_cols : A->mInfo.num_rows,
		alpha, A->data, A->mInfo.stride, B->data, B->mInfo.stride,
		beta, out->data, out->mInfo.stride);
#else
	BaseFloat result = 0.0f ;
	Matrix A1, B1, C;
	ResizeMatrix(A->mInfo.num_rows, A->mInfo.num_cols, kDefaultStride, &A1);
	for(int i = 0;  i < A1.mInfo.num_rows; i++)
		memcpy(A1.data + i * A1.mInfo.num_cols, A->data + i * A->mInfo.stride, A1.mInfo.num_cols * sizeof(BaseFloat));

	ResizeMatrix(B->mInfo.num_rows, B->mInfo.num_cols, kDefaultStride, &B1);
	for(int i = 0;  i < B1.mInfo.num_rows; i++)
		memcpy(B1.data + i * B1.mInfo.num_cols, B->data + i * B->mInfo.stride, B1.mInfo.num_cols * sizeof(BaseFloat));

	ResizeMatrix(A->mInfo.num_rows, B->mInfo.num_rows, kDefaultStride, &C);
	for(int i = 0 ; i < A->mInfo.num_rows ; i++)
	{
		for(int j = 0 ; j < B->mInfo.num_rows; j++)
		{
			result = 0.0f ;
			for(int k = 0 ; k < A->mInfo.num_cols ; ++k)  result += A1.data[i * A->mInfo.num_cols + k] * B1.data[ j* A->mInfo.num_cols + k] ;
			C.data[i * B->mInfo.num_rows + j] = result ;
		}
	}
	for (int i = 0;  i < C.mInfo.num_rows; i++)
		for(int j =0 ; j < C.mInfo.num_cols; j++)
			out->data[i * out->mInfo.stride + j] +=C.data[i * C.mInfo.num_cols + j];

	FreeMatrix(&A1);
	FreeMatrix(&B1);
	FreeMatrix(&C);
#endif // USE_OPENBLAS
}
// ###################### Above math use openblas interfaces End ######################

void VectorAdd(BaseFloat *vec, int32 len, BaseFloat c)
{
	for (int32 i = 0; i < len; ++i)
		vec[i] += c;
}


void ApplyFloor(BaseFloat *data, BaseFloat floor_val, int32 len)
{
	for (int32 i = 0; i < len; ++i)
		data[i] = (data[i] > floor_val) ? data[i] : floor_val;
}

void ApplyLog(BaseFloat *data, int32 len)
{
	for (int32 i = 0; i < len; ++i)
	{
		assert(data[i] >= 0.0);
		data[i] = logf(data[i]);
	}
}

BaseFloat VectorSum(BaseFloat *vec, int32 len)
{
	BaseFloat sum = 0;
	for (int32 i = 0; i < len; ++i)
		sum += vec[i];
	return sum;
}

void VectorMulElements(BaseFloat *vec1, BaseFloat *vec2, int32 len)
{
	for (int32 i = 0; i < len; ++i)
		vec1[i] *= vec2[i];
}

void MatApplyFloor(Matrix *M, const BaseFloat floor_val) 
{
	int32 num_rows = M->mInfo.num_rows, num_cols = M->mInfo.num_cols;
	for (int32 i = 0; i < num_rows; ++i) 
	{
		BaseFloat *data = M->data + i * M->mInfo.stride;
		for (int32 j = 0; j < num_cols; ++j)
			data[j] = (data[j] < floor_val ? floor_val : data[j]);
	}
}

void ComplexImExp(BaseFloat x, BaseFloat *a_re, BaseFloat *a_im)
{
	*a_re = cos(x);
	*a_im = sin(x);
}

void ComplexMul(const BaseFloat a_re, const BaseFloat a_im, BaseFloat *b_re, BaseFloat *b_im)
{
	BaseFloat tmp_re = (*b_re * a_re) - (*b_im * a_im);
	*b_im = *b_re * a_im + *b_im * a_re;
	*b_re = tmp_re;
}

void ComplexAddProduct(const BaseFloat a_re, const BaseFloat a_im,
	const BaseFloat b_re, const BaseFloat b_im,
	BaseFloat *c_re, BaseFloat *c_im)
{
	*c_re += b_re * a_re - b_im * a_im;
	*c_im += b_re * a_im + b_im * a_re;
}

bool ApproxEqual(BaseFloat a, BaseFloat b, BaseFloat relative_tolerance)
{
	// a==b handles infinities.
	if (a == b) return true;
	BaseFloat diff = fabs(a - b);
	if (diff == INFINITY || diff != diff) return false;  // diff is +inf or nan.
	return (diff <= relative_tolerance * (fabs(a) + fabs(b)));
}