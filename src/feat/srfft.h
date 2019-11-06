/*
* fft compute
*
* author: xiangxyq
* E-mail: xiangxyq@foxmail.com
*/

#ifndef SRFFT_H
#define SRFFT_H

#include "feat-params.h"


extern void SplitRadixFftInit();

extern void SplitRadixFftCompute(BaseFloat *data);

extern void FreeSrfftMemory();


#endif // SRFFT_H
