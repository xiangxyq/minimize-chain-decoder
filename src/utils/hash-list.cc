/*
* This File hash list operations
*
* author: xiangxyq
* E-mail: xiangxyq@foxmail.com
*/

#include <vector>
#include <cstdio>
#include <cassert>
#include <iostream>
#include "hash-list.h"

static Elem *list_head_; // head of currently stored list.
static Elem *freed_head_;  // head of list of currently freed elements. [ready for allocation]
static int32 bucket_list_tail_; // tail of list of active hash buckets.
static int32 hash_size_; // number of hash buckets.
static const int32 allocate_block_size_ = 1024;  // Number of Elements to
                                                 // allocate in one block.  Must be largish so storing allocated_ doesn't become a problem.
static std::vector<HashBucket> buckets_;
static std::vector<Elem *> allocated_;  // list of allocated blocks.

void InitHashList()
{
	list_head_ = NULL;
	bucket_list_tail_ = -1;  // invalid.
	hash_size_ = 0;
	freed_head_ = NULL;
}

void SetHashListSize(int32 size)
{
	hash_size_ = size;
	assert(list_head_ == NULL && bucket_list_tail_ == -1);  // make sure empty.

	if (size > buckets_.size())
		buckets_.resize(size);
}

int32 GetHashListSize()
{
	return hash_size_;
}

Elem *GetHashList()
{
	return list_head_;
}

Elem *HashListFindElem(int32 key)
{
	int32 index = static_cast<int32> (key) % hash_size_;

	HashBucket &bucket = buckets_[index];
	if (bucket.last_elem == NULL)
		return NULL;  // empty bucket.
	else 
	{
		Elem *head = (bucket.prev_bucket == static_cast<int32>(-1) ? list_head_ :
                buckets_[bucket.prev_bucket].last_elem->tail),*tail = bucket.last_elem->tail;
		for (Elem *e = head; e != tail; e = e->tail)
			if (e->key == key) return e;
		return NULL;  // Not found.
	}
}

Elem *NewHashListElem()
{
	if (freed_head_) {
		Elem *ans = freed_head_;
		freed_head_ = freed_head_->tail;
		return ans;
	}
	else {
		Elem *tmp = new Elem[allocate_block_size_];
		for (int i = 0; i + 1 < allocate_block_size_; i++)
			tmp[i].tail = tmp + i + 1;
		tmp[allocate_block_size_ - 1].tail = NULL;
		freed_head_ = tmp;
		allocated_.push_back(tmp);
		return NewHashListElem();
	}
}

void HashListInsertToken(int32 key, Token *val) 
{
	int32 index = (int32)key % hash_size_;
	HashBucket &bucket = buckets_[index];
	Elem *elem = NewHashListElem();
	elem->key = key;
	elem->val = val;

	if (bucket.last_elem == NULL)
	{  // Unoccupied bucket.  Insert at
	  // head of bucket list (which is tail of regular list, they go in
	  // opposite directions).
		if (bucket_list_tail_ == static_cast<int32>(-1))
		{
			// list was empty so this is the first elem.
			assert(list_head_ == NULL);
			list_head_ = elem;
		}
		else 
		{
			// link in to the chain of Elems
			buckets_[bucket_list_tail_].last_elem->tail = elem;
		}
		elem->tail = NULL;
		bucket.last_elem = elem;
		bucket.prev_bucket = bucket_list_tail_;
		bucket_list_tail_ = index;
	}
	else 
	{
		// Already-occupied bucket.  Insert at tail of list of elements within
		// the bucket.
		elem->tail = bucket.last_elem->tail;
		bucket.last_elem->tail = elem;
		bucket.last_elem = elem;
	}
}

Elem *ClearHashList()
{
	// Clears the hashtable and gives ownership of the currently contained list to the user.
	for (int32 cur_bucket = bucket_list_tail_; cur_bucket != static_cast<int32>(-1);
			cur_bucket = buckets_[cur_bucket].prev_bucket)
		buckets_[cur_bucket].last_elem = NULL;  // this is how we indicate "empty".

	bucket_list_tail_ = static_cast<int32>(-1);
	Elem *ans = list_head_;
	list_head_ = NULL;

	return ans;
}

void DeleteHashList(Elem *e) 
{
	e->tail = freed_head_;
	freed_head_ = e;
}

void FreeHashList() 
{
	// First test whether we had any memory leak within the
	// HashList, i.e. things for which the user did not call Delete().
	int32 num_in_list = 0, num_allocated = 0;
	for (Elem *e = freed_head_; e != NULL; e = e->tail)
		num_in_list++;
	for (int32 i = 0; i < allocated_.size(); ++i)
	{
		num_allocated += allocate_block_size_;
			delete[] allocated_[i];
	}
	buckets_.clear();
	allocated_.clear();

	if (num_in_list != num_allocated) 
	{
		std::cout << "Possible memory leak: " \
                  << num_in_list << " != " << num_allocated \
                  << " : you might have forgotten to call Delete on some Elems" << std::endl;
        exit(1);
	}

	InitHashList();
}
