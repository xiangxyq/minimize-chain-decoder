/*
* This File Packages various interfaces of ASR
*
* author: xiangxyq
* E-mail: xiangxyq@foxmail.com
*/

#ifndef ASR_APIS_H
#define ASR_APIS_H

#include <vector>
#include <iostream>
#include "../utils/asr-types.h"

extern void InitAsrModels();

extern void AsrDecodingBegin(std::vector<BaseFloat> &waveform);

extern bool AsrDecodingResult(std::vector<std::string> &result, BaseFloat *likelihood);

extern void AsrEnd();


#endif // ASR_APIS_H
