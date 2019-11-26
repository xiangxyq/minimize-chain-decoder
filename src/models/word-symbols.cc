/*
* Word id to string
*
* author: xiangxyq
* E-mail: xiangxyq@foxmail.com
*/

#include "words.h"
#include "word-symbols.h"

static const int32 words_id[] = WORDS_IDS;
static const char *words_str[] = WORDS_STRS;

// search word
int FindWord(const int32 *sortedSeq, int32 seqLength, int32 keyData)
{
	for (int32 i = 0; i < seqLength; ++i)
	{
		if (sortedSeq[i] == keyData)
			return i;
	}

	return -1;
}

#if 0
// binary search, must be sorted
int32 Find(const int32 *sortedSeq, int32 seqLength, int32 keyData)
{
	int32 low = 0, mid, high = seqLength - 1;

	while (low <= high)
	{
		mid = (low + high) / 2;
		if (keyData < sortedSeq[mid])
			high = mid - 1;
		else if (keyData > sortedSeq[mid])
			low = mid + 1;
		else
			return mid;
	}
	return -1;
}
#endif

const std::string WordId2Str(int word_id)
{
	int32 idx = FindWord(words_id, WORDS_LEN, word_id);
	if (idx == -1)
	{
		std::cout << "Could not find word id : " << word_id << std::endl;
		exit(1);
	}

	return words_str[idx];
}