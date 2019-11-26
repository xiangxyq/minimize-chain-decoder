/*
* This File support nnet3 compents calculation apis
*
* author: xiangxyq
* E-mail: xiangxyq@foxmail.com
*/

#ifndef COMPENTS_H
#define COMPENTS_H

#include "../utils/matrix.h"

extern void AffinePropagate(Matrix *in, Matrix *weights, BaseFloat *bias, int32 bias_len, Matrix *out);

extern void ReluPropagate(Matrix *in, Matrix *out);


#endif // COMPENTS_H