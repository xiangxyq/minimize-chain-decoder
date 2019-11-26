/*
* This File Packages various interfaces of ASR
*
* author: xiangxyq
* E-mail: xiangxyq@foxmail.com
*/

#include "asr-apis.h"
#include "../utils/hash-list.h"
#include "../feat/compute-feat.h"
#include "../models/fst.h"
#include "../models/word-symbols.h"
#include "decode/decode.h"
#include "../nnet3/chain-compute.h"

static std::vector<BaseFloat *> feats;

void InitAsrModels()
{
    ComputeFeatInit();
    InitDecoding();
}

void AsrDecodingBegin(std::vector<BaseFloat> &waveform)
{
    ComputeFeatures(waveform, feats);
    ResetDecoding();
    AdvanceDecoding(feats);
	FinalizeDecoding();
	FeatsFree(feats);
}

bool AsrDecodingResult(std::vector<std::string> &result, BaseFloat *likelihood)
{
    std::vector<int32> alignment, words_seq;

	if (!GetLinearSymbolSequence(&alignment, &words_seq, likelihood))
		return false;

	for (int i = 0; i < words_seq.size(); ++i)
		result.push_back(WordId2Str(words_seq[i]));

	BaseFloat per_likelihood = - (*likelihood) / alignment.size();  //compute average likelihood
	std::cout << "likelihood per frame is : " << per_likelihood << " over " << alignment.size() << " frames." << std::endl;

    return true;
}

void AsrEnd()
{
    FeatDestroy(feats);
    DestroyDecode();
    FreeChianCompute();
    FreeHashList();
}