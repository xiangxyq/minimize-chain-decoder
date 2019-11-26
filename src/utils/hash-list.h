/*
* This File hash list operations
*
* author: xiangxyq
* E-mail: xiangxyq@foxmail.com
*/

#ifndef HASH_LIST_H
#define HASH_LIST_H

#include "tokens.h"
#include "asr-types.h"

struct Elem
{
	int key;
	Token *val;
	Elem *tail;
};

typedef struct
{
	int32 prev_bucket;  // index to next bucket (-1 if list tail).  Note:
	// list of buckets goes in opposite direction to list of Elems.
	Elem *last_elem;  // pointer to last element in this bucket (NULL if empty)
	//inline HashBucket(int i, struct Elem *e) : prev_bucket(i), last_elem(e) {}
} HashBucket;


extern void InitHashList();

extern void SetHashListSize(int32 size);

extern int32 GetHashListSize();

extern Elem *GetHashList();

extern Elem *HashListFindElem(int32 key);

extern Elem *NewHashListElem();

extern void HashListInsertToken(int32 key, Token *val);

extern Elem *ClearHashList();

extern void DeleteHashList(Elem *e);

extern void FreeHashList();


#endif // HASH_LIST_H