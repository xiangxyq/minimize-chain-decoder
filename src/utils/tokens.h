/*
* This File token operations
*
* author: xiangxyq
* E-mail: xiangxyq@foxmail.com
*/

#ifndef TOKEN_H
#define TOKEN_H

#include "asr-types.h"

typedef struct StdToken Token;
typedef struct ForwardLink ForwardLink_T;

// Standard token type for LatticeFasterDecoder.  Each active HCLG
// (decoding-graph) state on each frame has one token.
struct StdToken
{
	// tot_cost is the total (LM + acoustic) cost from the beginning of the
	// utterance up to this point.  (but see cost_offset_, which is subtracted
	// to keep it in a good numerical range).
	BaseFloat tot_cost;

	// exta_cost is >= 0.  After calling PruneForwardLinks, this equals
	// the minimum difference between the cost of the best path, and the cost of
	// this is on, and the cost of the absolute best path, under the assumption
	// that any of the currently active states at the decoding front may
	// eventually succeed (e.g. if you were to take the currently active states
	// one by one and compute this difference, and then take the minimum).
	BaseFloat extra_cost;

	// 'links' is the head of singly-linked list of ForwardLinks, which is what we
	// use for lattice generation.
	ForwardLink_T *links;

	//'next' is the next in the singly-linked list of tokens for this frame.
	StdToken *next;
};

struct ForwardLink
{
	Token *next_tok;  // the next token [or NULL if represents final-state]
	int32 ilabel;  // ilabel on arc
	int32 olabel;  // olabel on arc
	BaseFloat graph_cost;  // graph cost of traversing arc (contains LM, etc.)
	BaseFloat acoustic_cost;  // acoustic cost (pre-scaled) of traversing arc
	ForwardLink_T *next;  // next in singly-linked list of forward arcs (arcs in the state-level lattice) from a token.
};

// head of per-frame list of Tokens (list is in topological order),
// and something saying whether we ever pruned it using PruneForwardLinks.
typedef struct 
{
	Token *toks;
	bool must_prune_forward_links;
	bool must_prune_tokens;
} TokenList;


extern Token *AddToken(BaseFloat tot_cost, BaseFloat extra_cost, ForwardLink *links, Token *next);

extern ForwardLink *AddForwardLink(Token *next_tok, int32 ilabel, int32 olabel,
					BaseFloat graph_cost, BaseFloat acoustic_cost, ForwardLink *next);


#endif // TOKEN_H