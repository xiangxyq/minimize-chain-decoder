/*
* core decode file
*
* author: xiangxyq
* E-mail: xiangxyq@foxmail.com
*/

#ifndef DECODE_H
#define DECODE_H

#include <vector>
#include "../utils/asr-types.h"
#include "../feat/compute-feat.h"

//#define DEBUG

extern void InitDecoding();

extern void ResetDecoding();

extern void AdvanceDecoding(const std::vector<BaseFloat* > &feats);

extern void FinalizeDecoding();

extern bool GetLinearSymbolSequence(std::vector<int32> *isymbols_out, std::vector<int32> *osymbols_out, BaseFloat *tot_weight_out);

extern void DestroyDecode();


#endif // DECODE_H