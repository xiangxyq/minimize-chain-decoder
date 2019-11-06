/*
* math apis
*
* author: xiangxyq
* E-mail: xiangxyq@foxmail.com
*/

#ifndef ASR_MATH_H
#define ASR_MATH_H

#include "asr-types.h"

#ifndef M_PI
#define M_PI       (3.1415926535897932384626433832795)
#endif

#ifndef M_2PI
#define M_2PI      (6.283185307179586476925286766559005)
#endif

#ifndef M_SQRT1_2
#define M_SQRT1_2 0.7071067811865475244008443621048490
#endif

typedef enum
{
	kTrans = 112, // = CblasTrans
	kNoTrans = 111  // = CblasNoTrans
} MatrixTransposeType;

// ###################### Below math use openblas interfaces Start ######################
extern BaseFloat VecVec(const BaseFloat *a, const BaseFloat *b, int32 len);

extern void AddMatVec(const BaseFloat alpha, const BaseFloat *M, int32 M_rows, int32 M_cols, int32 stride,
	MatrixTransposeType trans, const BaseFloat *v, int32 v_len, const BaseFloat beta, BaseFloat *result, int32 result_len);
// ###################### Below math use openblas interfaces End ######################

extern void VectorAdd(BaseFloat *vec, int32 len, BaseFloat c);

extern void ApplyLog(BaseFloat *data, int32 len);

extern void ApplyFloor(BaseFloat *data, BaseFloat floor_val, int32 len);

extern BaseFloat VectorSum(BaseFloat *vec, int32 len);

extern void VectorMulElements(BaseFloat *vec1, BaseFloat *vec2, int32 len);

extern void ComplexImExp(BaseFloat x, BaseFloat *a_re, BaseFloat *a_im);

//! ComplexMul implements, inline, the complex multiplication b *= a.
extern void ComplexMul(const BaseFloat a_re, const BaseFloat a_im, BaseFloat *b_re, BaseFloat *b_im);

extern void ComplexAddProduct(const BaseFloat a_re, const BaseFloat a_im,
	const BaseFloat b_re, const BaseFloat b_im,
	BaseFloat *c_re, BaseFloat *c_im);


#endif //ASR_MATH_H
