/*
* core decode file
*
* author: xiangxyq
* E-mail: xiangxyq@foxmail.com
*/

#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <iostream>
#include <cassert>
#include "decode.h"
#include "../configs/decode_conf.h"
#include "../nnet3/chain-compute.h"
#include "../utils/asr-math.h"
#include "../utils/tokens.h"
#include "../models/fst.h"
#include "../utils/hash-list.h"


static BaseFloat beam = BEAM;
static int32 max_active = MAX_ACTIVE;
static int32 min_active = MIN_ACTIVE;
static BaseFloat  lattice_beam = LATTICE_BEAM;
static int32 prune_interval = PRUNE_INTERVAL;
static BaseFloat beam_delta = BEAM_DELTA; // has nothing to do with beam_ratio
static BaseFloat hash_ratio = HASH_RATIO;
static BaseFloat prune_scale = PRUNE_SCALE;

static Fst_Strcut_Def *hclg_, *raw_fst_, *best_fst_;
static std::vector<BaseFloat> cost_offsets_;	 // This contains, for each
								  // frame, an offset that was added to the acoustic log-likelihoods on that
								  // frame in order to keep everything in a nice dynamic range i.e.  close to
								  // zero, to reduce roundoff errors.
static int32 num_toks_ = 0;
static bool decoding_finalized_;
static std::vector<TokenList> active_toks_;; // Lists of tokens, indexed by frame (members of TokenList are toks, must_prune_forward_links, must_prune_tokens).
static bool warned_;
/// For the meaning of the next 3 variables, see the comment for
/// decoding_finalized_ above., and ComputeFinalCosts().
static std::unordered_map<Token*, BaseFloat> final_costs_;
static BaseFloat final_relative_cost_;
static BaseFloat final_best_cost_;
static std::vector<StateId> queue_;  // temp variable used in ProcessNonemitting,
static std::vector<BaseFloat> tmp_array_; // used in GetCutoff.

static void DeleteElems(Elem *list);
static void DeleteForwardLinks(Token *tok);
static void ClearActiveTokens();
static int32 NumFramesDecoded();
static void ProcessNonemitting(BaseFloat cutoff);
static Token* FindOrAddToken(int32 state, int32 frame_plus_one, BaseFloat tot_cost, Token *backpointer, bool *changed);
static void PruneTokensForFrame(int32 frame_plus_one);
static void PruneForwardLinks(int32 frame_plus_one, bool *extra_costs_changed, bool *links_pruned, BaseFloat delta);
static void PruneActiveTokens(BaseFloat delta);
static BaseFloat ProcessEmitting(const std::vector<BaseFloat* > &feats);
static BaseFloat GetCutoff(Elem *list_head, int *tok_count, BaseFloat *adaptive_beam, Elem **best_elem);
static void PossiblyResizeHash(int num_toks);
static void ComputeFinalCosts(std::unordered_map<Token*, BaseFloat> *final_costs,
											BaseFloat *final_relative_cost, BaseFloat *final_best_cost);
static void TopSortTokens(Token *tok_list, std::vector<Token*> *topsorted_list);
static void PruneForwardLinksFinal();
static void GetBestPath();


void InitDecoding()
{
	assert(beam > 0.0 && max_active > 1 && lattice_beam > 0.0
			&& min_active <= max_active
			&& prune_interval > 0 && beam_delta > 0.0 && hash_ratio >= 1.0
			&& prune_scale > 0.0 && prune_scale < 1.0);
 
	InitHashList();
	InitChainCompute();

	InitFst();
	hclg_ = GetHclg();
	raw_fst_ = GetRawFst();
	best_fst_ = GetBestFst();
}

void ResetDecoding()
{
	InitChainCompute();
	// clean up from last time:;
	Elem *last_elem = ClearHashList();
	DeleteElems(last_elem);
	cost_offsets_.clear();
	queue_.clear();
	tmp_array_.clear();
	ClearActiveTokens();

	DeleteFstStates(raw_fst_);
	DeleteFstStates(best_fst_);
	warned_ = false;
	num_toks_ = 0;
	decoding_finalized_ = false;
	final_costs_.clear();	
	FreeHashList();

	StateId start_state = GetFstStartState(hclg_);
	assert(start_state != kNoStateId);
	active_toks_.resize(1);
	Token *start_tok = AddToken(0.0f, 0.0f, NULL, NULL);
	active_toks_[0].toks = start_tok;
	SetHashListSize(1000);
	HashListInsertToken(start_state, start_tok);
	num_toks_++;

	ProcessNonemitting(beam);
}



static void ProcessNonemitting(BaseFloat cutoff)
{
	assert(!active_toks_.empty());

	int32 frame = static_cast<int32>(active_toks_.size()) - 2;
	// Note: "frame" is the time-index we just processed, or -1 if
	// we are processing the nonemitting transitions before the
	// first frame (called from ResetDecoding()).

	// Processes nonemitting arcs for one frame.  Propagates within toks_.
	// Note-- this queue structure is is not very optimal as
	// it may cause us to process states unnecessarily (e.g. more than once),
	// but in the baseline code, turning this vector into a set to fix this
	// problem did not improve overall speed.
	assert(queue_.empty());

	if (GetHashList() == NULL) 
	{
		if (!warned_) {
			std::cout << "Error, no surviving tokens: frame is " << frame << std::endl;
			warned_ = true;
		}
	}

	for (const Elem *e = GetHashList(); e != NULL; e = e->tail)
	{
		StateId state = e->key;
		if (GetFstNumInputEpsilons(state, hclg_) != 0)
			queue_.push_back(state);
	}

	while (!queue_.empty())
	{
		StateId state = queue_.back();
		queue_.pop_back();

		Token *tok = HashListFindElem(state)->val;  // would segfault if state not in toks_ but this can't happen.
		BaseFloat cur_cost = tok->tot_cost;
		if (cur_cost > cutoff) // Don't bother processing successors.
			continue;
		// If "tok" has any existing forward links, delete them,
		// because we're about to regenerate them.  This is a kind
		// of non-optimality (remember, this is the simple decoder),
		// but since most states are emitting it's not a huge issue.
		DeleteForwardLinks(tok); // necessary when re-visiting
		tok->links = NULL;

		const State *state_node = GetFstState(state, hclg_);
		for (int i = 0; i < state_node->num_arc; ++i)
		{
			if (0 == state_node->arc[i].ilabel)
			{
				BaseFloat graph_cost = state_node->arc[i].arc_weight,
					tot_cost = cur_cost + graph_cost;
				if (tot_cost < cutoff) 
				{
					bool changed;
					Token *new_tok = FindOrAddToken(state_node->arc[i].next_state, frame + 1, tot_cost, tok, &changed);
					tok->links = AddForwardLink(new_tok, 0, state_node->arc[i].olabel, graph_cost, 0, tok->links);

					// "changed" tells us whether the new token has a different
					// cost from before, or is new [if so, add into queue].
					if (changed && GetFstNumInputEpsilons(state_node->arc[i].next_state, hclg_) != 0)
						queue_.push_back(state_node->arc[i].next_state);
				}
			}
		}
	} // while queue not empty
}

static BaseFloat ProcessEmitting(const std::vector<BaseFloat* > &feats)
{
	assert(active_toks_.size() > 0);

	int32 frame = active_toks_.size() - 1; // frame is the frame-index
										 // (zero-based) used to get likelihoods
										 // from the decodable object.
	active_toks_.resize(active_toks_.size() + 1);

	Elem *final_toks = ClearHashList(); // analogous to swapping prev_toks_ / cur_toks_
											   // in simple-decoder.h.   Removes the Elems from
									           // being indexed in the hash in toks_.
	Elem *best_elem = NULL;
	BaseFloat adaptive_beam;
	int tok_cnt;
	BaseFloat cur_cutoff = GetCutoff(final_toks, &tok_cnt, &adaptive_beam, &best_elem);

#ifdef DEBUG
	std::cout << "Adaptive beam on frame " << NumFramesDecoded() << " is " <<  adaptive_beam << std::endl;
#endif

	PossiblyResizeHash(tok_cnt);  // This makes sure the hash is always big enough.
	BaseFloat next_cutoff = std::numeric_limits<BaseFloat>::infinity(); // pruning "online" before having seen all tokens
	BaseFloat cost_offset = 0.0f; // Used to keep probabilities in a good dynamic range.

	// First process the best token to get a hopefully
	// reasonably tight bound on the next cutoff.  The only
	// products of the next block are "next_cutoff" and "cost_offset".
	if (best_elem) 
	{
		StateId state = best_elem->key;
		Token *tok = best_elem->val;
		cost_offset = - tok->tot_cost;

		const State *state_node = GetFstState(state, hclg_);
		for (int i = 0; i < state_node->num_arc; ++i)
		{
			if (0 != state_node->arc[i].ilabel)
			{
				BaseFloat new_weight = state_node->arc[i].arc_weight + cost_offset -
										LogLikelihood(feats, frame, state_node->arc[i].ilabel) + tok->tot_cost;
				if (new_weight + adaptive_beam < next_cutoff)
					next_cutoff = new_weight + adaptive_beam;
			}
		}
	}
	// Store the offset on the acoustic likelihoods that we're applying.
	// Could just do cost_offsets_.push_back(cost_offset), but we
	// do it this way as it's more robust to future code changes.
  	cost_offsets_.resize(frame + 1, 0.0);
  	cost_offsets_[frame] = cost_offset;

	// the tokens are now owned here, in final_toks, and the hash is empty.
	// 'owned' is a complex thing here; the point is we need to call DeleteElem
	// on each elem 'e' to let toks_ know we're done with them.
	for (Elem *e = final_toks, *e_tail; e != NULL; e = e_tail) 
	{
		// loop this way because we delete "e" as we go.
		StateId state = e->key;
		Token *tok = e->val;
		if (tok->tot_cost <= cur_cutoff) 
		{
			const State *state_node = GetFstState(state, hclg_);
			for (int i = 0; i < state_node->num_arc; ++i)
			{
				if (0 != state_node->arc[i].ilabel)
				{
					BaseFloat ac_cost = cost_offset -
						LogLikelihood(feats, frame, state_node->arc[i].ilabel),
						graph_cost = state_node->arc[i].arc_weight,
						cur_cost = tok->tot_cost,
						tot_cost = cur_cost + ac_cost + graph_cost;
					if (tot_cost > next_cutoff)
						continue;
					else if (tot_cost + adaptive_beam < next_cutoff)
						next_cutoff = tot_cost + adaptive_beam; // prune by best current token
																// Note: the frame indexes into active_toks_ are one-based, hence the + 1.
					Token *next_tok = FindOrAddToken(state_node->arc[i].next_state, frame + 1, tot_cost, tok, NULL);
					// NULL: no change indicator needed
					// Add ForwardLink from tok to next_tok (put on head of list tok->links)
					tok->links = AddForwardLink(next_tok, state_node->arc[i].ilabel, state_node->arc[i].olabel, graph_cost, ac_cost, tok->links);
				}
			}
		}
		e_tail = e->tail;
		DeleteHashList(e);
	}
	return next_cutoff;
}

static Token* FindOrAddToken(int state, int frame_plus_one, BaseFloat tot_cost,Token *backpointer, bool *changed) 
{
	// Returns the Token pointer.  Sets "changed" (if non-NULL) to true
	// if the token was newly created or the cost changed.
	assert(frame_plus_one < active_toks_.size());

	Token **toks = &(active_toks_[frame_plus_one].toks);
	Elem *e_found = HashListFindElem(state);
	if (e_found == NULL)
	{  // no such token presently.
		const BaseFloat extra_cost = 0.0;
		// tokens on the currently final frame have zero extra_cost
		// as any of them could end up
		// on the winning path.
		Token *new_tok = AddToken(tot_cost, extra_cost, NULL, *toks);

		// NULL: no forward links yet
		*toks = new_tok;
		num_toks_++;
		HashListInsertToken(state, new_tok);
		if (changed) *changed = true;
		return new_tok;
	}
	else 
	{
		Token *tok = e_found->val;  // There is an existing Token for this state.
		if (tok->tot_cost > tot_cost)
		{  // replace old token
			tok->tot_cost = tot_cost;
			if (changed) *changed = true;
		}
		else {
			if (changed) *changed = false;
		}
		return tok;
	}
}

static int NumFramesDecoded() 
{
	return active_toks_.size() - 1; 
}

void AdvanceDecoding(const std::vector<BaseFloat* > &feats)
{
	int max_num_frames = -1;  // 默认max_num_frames = -1
	assert(!active_toks_.empty() && !decoding_finalized_); //"You must call ResetDecoding() before AdvanceDecoding"

	int num_frames_ready = NumFramesDecodeReady(feats);// num_frames_ready must be >= num_frames_decoded, or else
													  // the number of frames ready must have decreased (which doesn't
													  // make sense) or the decodable object changed between calls (which isn't allowed).

	assert(num_frames_ready >= NumFramesDecoded());
	int target_frames_decoded = num_frames_ready;
	if (max_num_frames >= 0)
		target_frames_decoded = (target_frames_decoded < (NumFramesDecoded() + max_num_frames)) ? target_frames_decoded : (NumFramesDecoded() + max_num_frames) ;

	while (NumFramesDecoded() < target_frames_decoded) 
	{
		if (NumFramesDecoded() % prune_interval == 0)
			PruneActiveTokens(lattice_beam * prune_scale);
		BaseFloat cost_cutoff = ProcessEmitting(feats);
		ProcessNonemitting(cost_cutoff);
	}
}

void FinalizeDecoding()
{
	int final_frame_plus_one = NumFramesDecoded();
#ifdef DEBUG
	int32 num_toks_begin = num_toks_;
#endif
	// PruneForwardLinksFinal() prunes final frame (with final-probs), and
	// sets decoding_finalized_.
	PruneForwardLinksFinal();
	for (int32 f = final_frame_plus_one - 1; f >= 0; f--)
	{
		bool b1, b2; // values not used.
		BaseFloat dontcare = 0.0f; // delta of zero means we must always update
		PruneForwardLinks(f, &b1, &b2, dontcare);
		PruneTokensForFrame(f + 1);
	}
	PruneTokensForFrame(0);
#ifdef DEBUG
	std::cout << "pruned tokens from " << num_toks_begin << " to " << num_toks_ << std::endl;
#endif
}

void PruneActiveTokens(BaseFloat delta) 
{
	int32 cur_frame_plus_one = NumFramesDecoded();
#ifdef DEBUG
	int32 num_toks_begin = num_toks_;
#endif
	// The index "f" below represents a "frame plus one", i.e. you'd have to subtract
	// one to get the corresponding index for the decodable object.
	for (int32 f = cur_frame_plus_one - 1; f >= 0; f--) 
	{
		// Reason why we need to prune forward links in this situation:
		// (1) we have never pruned them (new TokenList)
		// (2) we have not yet pruned the forward links to the next f,
		// after any of those tokens have changed their extra_cost.
		if (active_toks_[f].must_prune_forward_links) 
		{
			bool extra_costs_changed = false, links_pruned = false;
			PruneForwardLinks(f, &extra_costs_changed, &links_pruned, delta);
			if (extra_costs_changed && f > 0) // any token has changed extra_cost
				active_toks_[f - 1].must_prune_forward_links = true;
			if (links_pruned) // any link was pruned
				active_toks_[f].must_prune_tokens = true;
			active_toks_[f].must_prune_forward_links = false; // job done
		}
		if (f + 1 < cur_frame_plus_one &&      // except for last f (no forward links)
			active_toks_[f + 1].must_prune_tokens) {
			PruneTokensForFrame(f + 1);
			active_toks_[f + 1].must_prune_tokens = false;
		}
	}
#ifdef DEBUG
	std::cout << "PruneActiveTokens: pruned tokens from " << num_toks_begin << " to " << num_toks_ << std::endl;
#endif
}

void PruneForwardLinks(int frame_plus_one, bool *extra_costs_changed,
	bool *links_pruned, BaseFloat delta) 
{
	// delta is the amount by which the extra_costs must change
	// If delta is larger,  we'll tend to go back less far
	//    toward the beginning of the file.
	// extra_costs_changed is set to true if extra_cost was changed for any token
	// links_pruned is set to true if any link in any token was pruned

	*extra_costs_changed = false;
	*links_pruned = false;

	assert(frame_plus_one >= 0 && frame_plus_one < active_toks_.size());

	if (active_toks_[frame_plus_one].toks == NULL) 
	{  // empty list; should not happen.
		if (!warned_) 
		{
			std::cout << "No tokens alive [doing pruning].. warning first time only for each utterance" << std::endl;
			warned_ = true;
		}
	}

	// We have to iterate until there is no more change, because the links
	// are not guaranteed to be in topological order.
	bool changed = true;  // difference new minus old extra cost >= delta ?
	while (changed) 
	{
		changed = false;
		for (Token *tok = active_toks_[frame_plus_one].toks; tok != NULL; tok = tok->next)
		{
			ForwardLink *link, *prev_link = NULL;
			// will recompute tok_extra_cost for tok.
			BaseFloat tok_extra_cost = std::numeric_limits<BaseFloat>::infinity();
			// tok_extra_cost is the best (min) of link_extra_cost of outgoing links
			for (link = tok->links; link != NULL; )
			{
				// See if we need to excise this link...
				Token *next_tok = link->next_tok;
				BaseFloat link_extra_cost = next_tok->extra_cost +
					((tok->tot_cost + link->acoustic_cost + link->graph_cost)
						- next_tok->tot_cost);  // difference in brackets is >= 0
				   // link_exta_cost is the difference in score between the best paths
				   // through link source state and through link destination state

				assert(link_extra_cost == link_extra_cost);  // check for NaN

				if (link_extra_cost > lattice_beam) 
				{  // excise link
					ForwardLink *next_link = link->next;
					if (prev_link != NULL) prev_link->next = next_link;
					else tok->links = next_link;
					delete link;
					link = next_link;  // advance link but leave prev_link the same.
					*links_pruned = true;
				}
				else
				{   // keep the link and update the tok_extra_cost if needed.
					if (link_extra_cost < 0.0) // this is just a precaution.
					{ 
						if (link_extra_cost < -0.01)
							std::cout << "Negative extra_cost: " << link_extra_cost << std::endl;
						link_extra_cost = 0.0;
					}
					if (link_extra_cost < tok_extra_cost)
						tok_extra_cost = link_extra_cost;
					prev_link = link;  // move to next link
					link = link->next;
				}
			}  // for all outgoing links
			if (fabs(tok_extra_cost - tok->extra_cost) > delta)
				changed = true;   // difference new minus old is bigger than delta
			tok->extra_cost = tok_extra_cost;	// will be +infinity or <= lattice_beam_.
												// infinity indicates, that no forward link survived pruning
		}  // for all Token on active_toks_[frame]
		if (changed) *extra_costs_changed = true;

		// Note: it's theoretically possible that aggressive compiler
		// optimizations could cause an infinite loop here for small delta and
		// high-dynamic-range scores.
	} // while changed
}

void PruneTokensForFrame(int frame_plus_one)
{
	assert(frame_plus_one >= 0 && frame_plus_one < active_toks_.size());

	Token **toks = &(active_toks_[frame_plus_one].toks);
	if (*toks == NULL)
		std::cout << "No tokens alive [doing pruning]" << std::endl;
	Token *tok, *next_tok, *prev_tok = NULL;
	for (tok = *toks; tok != NULL; tok = next_tok) 
	{
		next_tok = tok->next;
		if (tok->extra_cost == std::numeric_limits<BaseFloat>::infinity())
		{
			// token is unreachable from end of graph; (no forward links survived)
			// excise tok from list and delete tok.
			if (prev_tok != NULL) prev_tok->next = tok->next;
			else *toks = tok->next;
			delete tok;
			num_toks_--;
		}
		else 
		{  // fetch next Token
			prev_tok = tok;
		}
	}
}

BaseFloat GetCutoff(Elem *list_head, int *tok_count,
						BaseFloat *adaptive_beam, Elem **best_elem) 
{
	BaseFloat best_weight = std::numeric_limits<BaseFloat>::infinity();
	// positive == high cost == bad.
	int32 count = 0;
	if (max_active == std::numeric_limits<int32>::max() && min_active == 0)
	{
		for (Elem *e = list_head; e != NULL; e = e->tail, count++) 
		{
			BaseFloat w = (BaseFloat)(e->val->tot_cost);
			if (w < best_weight) 
			{
				best_weight = w;
				if (best_elem) *best_elem = e;
			}
		}
		if (tok_count != NULL) *tok_count = count;
		if (adaptive_beam != NULL) *adaptive_beam = beam;
		return best_weight + beam;
	}
	else 
	{
    	tmp_array_.clear();
		for (Elem *e = list_head; e != NULL; e = e->tail, count++)
		{
			BaseFloat w = e->val->tot_cost;
			tmp_array_.push_back(w);
			if (w < best_weight) {
				best_weight = w;
				if (best_elem) *best_elem = e;
			}
		}
		if (tok_count != NULL) *tok_count = count;

		BaseFloat beam_cutoff = best_weight + beam,
			min_active_cutoff = std::numeric_limits<BaseFloat>::infinity(),
			max_active_cutoff = std::numeric_limits<BaseFloat>::infinity();

		if (tmp_array_.size() > max_active) 
		{
			std::nth_element(tmp_array_.begin(),
                       tmp_array_.begin() + max_active,
                       tmp_array_.end());
     		max_active_cutoff = tmp_array_[max_active];
		}
		if (max_active_cutoff < beam_cutoff) // max_active is tighter than beam.
		{ 
			if (adaptive_beam)
				*adaptive_beam = max_active_cutoff - best_weight + beam_delta;
			return max_active_cutoff;
		}
		if (tmp_array_.size() > min_active) 
		{
			if (min_active == 0) min_active_cutoff = best_weight;
			else 
			{
				std::nth_element(tmp_array_.begin(),
                         tmp_array_.begin() + min_active,
                         tmp_array_.size() > max_active ?
                         tmp_array_.begin() + max_active :
                         tmp_array_.end());
        		min_active_cutoff = tmp_array_[min_active];
			}
		}
		if (min_active_cutoff > beam_cutoff) // min_active is looser than beam.
		{ 
			if (adaptive_beam)
				*adaptive_beam = min_active_cutoff - best_weight + beam_delta;
			return min_active_cutoff;
		}
		else 
		{
			*adaptive_beam = beam;
			return beam_cutoff;
		}
	}
}

void PossiblyResizeHash(int num_toks) 
{
	int32 new_sz = (int32)((BaseFloat)num_toks * hash_ratio);
	if (new_sz > GetHashListSize()) 
		SetHashListSize(new_sz);
}


void PruneForwardLinksFinal() 
{
	assert(!active_toks_.empty());
	int frame_plus_one = active_toks_.size() - 1;

	if (active_toks_[frame_plus_one].toks == NULL)  // empty list; should not happen.
	{
		std::cout << "No tokens alive at end of file" << std::endl;
		exit(1);
	}

	typedef typename std::unordered_map<Token*, BaseFloat>::const_iterator IterType;
	ComputeFinalCosts(&final_costs_, &final_relative_cost_, &final_best_cost_);
	decoding_finalized_ = true;
	// We call DeleteElems() as a nicety, not because it's really necessary;
	// otherwise there would be a time, after calling PruneTokensForFrame() on the
	// final frame, when toks_.GetList() or toks_.Clear() would contain pointers
	// to nonexistent tokens.
	DeleteElems(ClearHashList());

	// Now go through tokens on this frame, pruning forward links...  may have to
	// iterate a few times until there is no more change, because the list is not
	// in topological order.  This is a modified version of the code in
	// PruneForwardLinks, but here we also take account of the final-probs.
	bool changed = true;
	BaseFloat delta = 1.0e-05f;
	while (changed) 
	{
		changed = false;
		for (Token *tok = active_toks_[frame_plus_one].toks;
			tok != NULL; tok = tok->next) {
			ForwardLink *link, *prev_link = NULL;
			// will recompute tok_extra_cost.  It has a term in it that corresponds
			// to the "final-prob", so instead of initializing tok_extra_cost to infinity
			// below we set it to the difference between the (score+final_prob) of this token,
			// and the best such (score+final_prob).
			BaseFloat final_cost;
			if (final_costs_.empty())
				final_cost = 0.0f;
			else
			{
       	 		IterType iter = final_costs_.find(tok);
        		if (iter != final_costs_.end())
          			final_cost = iter->second;
        		else
          			final_cost = std::numeric_limits<BaseFloat>::infinity();
			}
			BaseFloat tok_extra_cost = tok->tot_cost + final_cost - final_best_cost_;
			// tok_extra_cost will be a "min" over either directly being final, or
			// being indirectly final through other links, and the loop below may
			// decrease its value:
			for (link = tok->links; link != NULL; )
			{
				// See if we need to excise this link...
				Token *next_tok = link->next_tok;
				BaseFloat link_extra_cost = next_tok->extra_cost +
					((tok->tot_cost + link->acoustic_cost + link->graph_cost)
						- next_tok->tot_cost);
				if (link_extra_cost > lattice_beam)
				{  // excise link
					ForwardLink *next_link = link->next;
					if (prev_link != NULL) prev_link->next = next_link;
					else tok->links = next_link;
					delete link;
					link = next_link; // advance link but leave prev_link the same.
				}
				else
				{ // keep the link and update the tok_extra_cost if needed.
					if (link_extra_cost < 0.0)
					{ // this is just a precaution.
						if (link_extra_cost < -0.01)
							std::cout << "WARN: Negative extra_cost: " << link_extra_cost << std::endl;
						link_extra_cost = 0.0;
					}
					if (link_extra_cost < tok_extra_cost)
						tok_extra_cost = link_extra_cost;
					prev_link = link;
					link = link->next;
				}
			}
			// prune away tokens worse than lattice_beam above best path.  This step
			// was not necessary in the non-final case because then, this case
			// showed up as having no forward links.  Here, the tok_extra_cost has
			// an extra component relating to the final-prob.
			if (tok_extra_cost > lattice_beam)
				tok_extra_cost = std::numeric_limits<BaseFloat>::infinity();
			// to be pruned in PruneTokensForFrame

			if (!ApproxEqual(tok->extra_cost, tok_extra_cost, delta))
				changed = true;
			tok->extra_cost = tok_extra_cost; // will be +infinity or <= lattice_beam_.
		}
	} // while changed
}


void ComputeFinalCosts(std::unordered_map<Token*, BaseFloat> *final_costs, 
									BaseFloat *final_relative_cost, BaseFloat *final_best_cost)
{
	assert(!decoding_finalized_);

	if (final_costs != NULL)
		final_costs->clear();
	const Elem *final_toks = GetHashList();
	BaseFloat infinity = std::numeric_limits<BaseFloat>::infinity();
	BaseFloat best_cost = infinity,
		best_cost_with_final = infinity;

	while (final_toks != NULL) 
	{
		StateId state = final_toks->key;
		Token *tok = final_toks->val;
		const Elem *next = final_toks->tail;
		BaseFloat final_cost = GetFstFinalWeight(state, hclg_);
		BaseFloat cost = tok->tot_cost,
			cost_with_final = cost + final_cost;
		best_cost = cost < best_cost ? cost : best_cost;
		best_cost_with_final = cost_with_final < best_cost_with_final ? cost_with_final : best_cost_with_final;
		if (final_costs != NULL && final_cost != infinity)
			(*final_costs)[tok] = final_cost;
		final_toks = next;
	}
	if (final_relative_cost != NULL) {
		if (best_cost == infinity && best_cost_with_final == infinity)
		{
			// Likely this will only happen if there are no tokens surviving.
			// This seems the least bad way to handle it.
			*final_relative_cost = infinity;
		}
		else {
			*final_relative_cost = best_cost_with_final - best_cost;
		}
	}
	if (final_best_cost != NULL) 
	{
		if (best_cost_with_final != infinity) // final-state exists.
		{ 
			*final_best_cost = best_cost_with_final;
		}
		else // no final-state exists.
		{ 
			*final_best_cost = best_cost;
		}
	}
}

void TopSortTokens(Token *tok_list, std::vector<Token*> *topsorted_list)
{
	std::unordered_map<Token*, int32> token2pos;
	typedef typename std::unordered_map<Token*, int32>::iterator IterType;
	int num_toks = 0;
	for (Token *tok = tok_list; tok != NULL; tok = tok->next)
		num_toks++;
	int cur_pos = 0;
	// We assign the tokens numbers num_toks - 1, ... , 2, 1, 0.
	// This is likely to be in closer to topological order than
	// if we had given them ascending order, because of the way
	// new tokens are put at the front of the list.
	for (Token *tok = tok_list; tok != NULL; tok = tok->next)
    	token2pos[tok] = num_toks - ++cur_pos;

  	std::unordered_set<Token*> reprocess;

	for (IterType iter = token2pos.begin(); iter != token2pos.end(); ++iter)
	{
    	Token *tok = iter->first;
   	 	int32 pos = iter->second;
		for (ForwardLink *link = tok->links; link != NULL; link = link->next)
		{
			if (link->ilabel == 0) 
			{
				// We only need to consider epsilon links, since non-epsilon links
				// transition between frames and this function only needs to sort a list
				// of tokens from a single frame.
				IterType following_iter = token2pos.find(link->next_tok);
				if (following_iter != token2pos.end()) // another token on this frame,so must consider it.
				{ 
					int next_pos = following_iter->second;
					if (next_pos < pos) // reassign the position of the next Token.
					{ 
            			following_iter->second = cur_pos++;
            			reprocess.insert(link->next_tok);
					}
				}
			}
		}
		// In case we had previously assigned this token to be reprocessed, we can
		// erase it from that set because it's "happy now" (we just processed it).
		reprocess.erase(tok);
	}

	int max_loop = 1000000, loop_count; // max_loop is to detect epsilon cycles.
	for (loop_count = 0; !reprocess.empty() && loop_count < max_loop; ++loop_count)
	{
    	std::vector<Token*> reprocess_vec;
    	for (typename std::unordered_set<Token*>::iterator iter = reprocess.begin(); iter != reprocess.end(); ++iter)
      		reprocess_vec.push_back(*iter);

		reprocess.clear();
		for (typename std::vector<Token*>::iterator iter = reprocess_vec.begin(); iter != reprocess_vec.end(); ++iter)
		{
			Token *tok = *iter;
			int32 pos = token2pos[tok];
			// Repeat the processing we did above (for comments, see above).
			for (ForwardLink *link = tok->links; link != NULL; link = link->next)
			{
				if (link->ilabel == 0) 
				{
					IterType following_iter = token2pos.find(link->next_tok);
          			if (following_iter != token2pos.end()) 
					{
            			int32 next_pos = following_iter->second;
           				if (next_pos < pos) 
						{
              				following_iter->second = cur_pos++;
              				reprocess.insert(link->next_tok);
            			}
          			}
				}
			}
		}
	}

	assert(loop_count < max_loop); //Epsilon loops exist in your decoding graph (this is not allowed!)

  	topsorted_list->clear();
  	topsorted_list->resize(cur_pos, NULL);  // create a list with NULLs in between.
  	for (IterType iter = token2pos.begin(); iter != token2pos.end(); ++iter)
    	(*topsorted_list)[iter->second] = iter->first;
}


void GetBestPath()
{
	if (GetFstStartState(raw_fst_) == kNoStateId) return;
	assert(GetFstStartState(raw_fst_) == 0);

	std::vector<std::pair<BaseFloat, StateId> > best_cost_and_pred(GetFstNumStates(raw_fst_) + 1);
	int32 superfinal = GetFstNumStates(raw_fst_);
	for (int s = 0; s <= GetFstNumStates(raw_fst_); ++s) 
	{
		best_cost_and_pred[s].first = std::numeric_limits<BaseFloat>::infinity();
		best_cost_and_pred[s].second = kNoStateId;
	}
	best_cost_and_pred[0].first = 0;
	for (int s = 0; s < GetFstNumStates(raw_fst_); ++s) 
	{
		BaseFloat my_cost = best_cost_and_pred[s].first;
		const State *state_node = GetFstState(s,raw_fst_);
		for (int j = 0; j < state_node->num_arc; ++j)
		{
			double arc_cost = state_node->arc[j].arc_weight,
				next_cost = my_cost + arc_cost;
			if (next_cost < best_cost_and_pred[state_node->arc[j].next_state].first) 
			{
				best_cost_and_pred[state_node->arc[j].next_state].first = next_cost;
				best_cost_and_pred[state_node->arc[j].next_state].second = s;
			}
		}
		BaseFloat final_cost = GetFstFinalWeight(s, raw_fst_),
			tot_final = my_cost + final_cost;
		if (tot_final < best_cost_and_pred[superfinal].first)
		{
			best_cost_and_pred[superfinal].first = tot_final;
			best_cost_and_pred[superfinal].second = s;
		}
	}

	std::vector<StateId> states; // states on best path.
	StateId cur_state = superfinal;
	while (cur_state != 0) 
	{
		StateId prev_state = best_cost_and_pred[cur_state].second;
		if (prev_state == kNoStateId) 
		{
			std::cout << "Warning: Failure in best-path algorithm for lattice (infinite costs?)" << std::endl;
			return; // return empty best-path.
		}
		states.push_back(prev_state);
		assert(cur_state != prev_state);

		cur_state = prev_state;
	}

	std::reverse(states.begin(), states.end());
	for (int i = 0; i < states.size(); ++i)
		AddFstState(best_fst_);
	for (int s = 0; s < states.size(); ++s) 
	{
		if (s == 0) SetFstStartState(s, best_fst_);
		if (static_cast<int32>(s + 1) < states.size()) // transition to next state.
		{ 
			bool have_arc = false;
			Arc cur_arc;

			const State *state_node = GetFstState(states[s], raw_fst_);
			for (int j = 0; j < state_node->num_arc; ++j)
			{
				if (state_node->arc[j].next_state == states[s + 1]) 
				{
					if (!have_arc ||
						state_node->arc[j].arc_weight < cur_arc.arc_weight) 
					{
						cur_arc = state_node->arc[j];
						have_arc = true;
					}
				
				}
			}

			assert(have_arc);

			Arc tmp;
			tmp.ilabel = cur_arc.ilabel;
			tmp.olabel = cur_arc.olabel;
			tmp.arc_weight = cur_arc.arc_weight;
			tmp.next_state = s + 1;
			AddFstArc(s, tmp, best_fst_);
		}
		else 
		{ // final-prob.
			SetFstFinalWeight(s, GetFstFinalWeight(states[s], raw_fst_), best_fst_);
		}
	}
}

bool GetRawFstTokens()
{
	std::unordered_map<Token*, BaseFloat> final_costs_local;
	const std::unordered_map<Token*, BaseFloat> &final_costs =
		(decoding_finalized_ ? final_costs_ : final_costs_local);
	if (!decoding_finalized_)
		ComputeFinalCosts(&final_costs_local, NULL, NULL);

	//DeleteRawFstStates();  //moved to ResetDecoding()
	// num-frames plus one (since frames are one-based, and we have an extra frame for the start-state).
	int32 num_frames = active_toks_.size() - 1;
	assert(num_frames > 0);

	const int32 bucket_count = num_toks_/2 + 3;
	std::unordered_map<Token*, StateId> tok_map(bucket_count);
	// First create all states.
 	std::vector<Token*> token_list;
	for (int f = 0; f <= num_frames; ++f) 
	{
		if (active_toks_[f].toks == NULL) 
		{
			std::cout << "GetRawLattice: no tokens active on frame " << f << ": not producing lattice." << std::endl;
			return false;
		}
		TopSortTokens(active_toks_[f].toks, &token_list);
		for (int i = 0; i < token_list.size(); ++i)
			if (token_list[i] != NULL)
				tok_map[token_list[i]] = AddFstState(raw_fst_);
	}

	// The next statement sets the start state of the output FST.  Because we
	// topologically sorted the tokens, state zero must be the start-state.
	SetFstStartState(0, raw_fst_);

	// Now create all arcs.
	for (int f = 0; f <= num_frames; ++f) 
	{
		for (Token *tok = active_toks_[f].toks; tok != NULL; tok = tok->next)
		{
			StateId cur_state = tok_map[tok];
			for (ForwardLink *l = tok->links;
				l != NULL;
				l = l->next) 
			{
        		typename std::unordered_map<Token*, StateId>::const_iterator
            										iter = tok_map.find(l->next_tok);
				assert(iter != tok_map.end());

				StateId nextstate = iter->second;
				BaseFloat cost_offset = 0.0f;
				if (l->ilabel != 0) 
				{  // emitting..
					assert(f >= 0 && f < cost_offsets_.size());
					cost_offset = cost_offsets_[f];
				}
				Arc arc;
				arc.ilabel = l->ilabel;
				arc.olabel = l->olabel;
				arc.arc_weight = Times(l->graph_cost, l->acoustic_cost - cost_offset);
				arc.next_state = nextstate;
				AddFstArc(cur_state, arc, raw_fst_);
			}
			if (f == num_frames)
			{
				if (!final_costs.empty())
				{
				    typename std::unordered_map<Token*, BaseFloat>::const_iterator
              												iter = final_costs.find(tok);				
					if (iter != final_costs.end())
						SetFstFinalWeight(cur_state, Times(iter->second, 0.0f), raw_fst_);
				}
				else 
				{
					SetFstFinalWeight(cur_state, Times(0.0f, 0.0f), raw_fst_);
				}
			}
		}
	}

	return (GetFstNumStates(raw_fst_) > 0);
}

bool GetLinearSymbolSequence(std::vector<int32> *isymbols_out, std::vector<int32> *osymbols_out, BaseFloat *tot_weight_out)
{
	BaseFloat tot_weight = 0.0f;
  	std::vector<int32> ilabel_seq;
  	std::vector<int32> olabel_seq;

	GetRawFstTokens();
	GetBestPath();
	StateId cur_state = GetFstStartState(best_fst_);
	if (cur_state == kNoStateId) 
	{  // empty sequence.
		if (tot_weight_out != NULL) *tot_weight_out = std::numeric_limits<BaseFloat>::infinity();
		return false;
	}
	while (1) 
	{
		Weight w = GetFstFinalWeight(cur_state, best_fst_);
		if (w != std::numeric_limits<BaseFloat>::infinity())
		{  // is final..
			tot_weight = Times(w, tot_weight);
			if (GetNumFstArcs(cur_state, best_fst_) != 0) return false;
			if (isymbols_out != NULL) *isymbols_out = ilabel_seq;
			if (osymbols_out != NULL) *osymbols_out = olabel_seq;
			if (tot_weight_out != NULL) *tot_weight_out = tot_weight;
			return true;
		}
		else 
		{
			if (GetNumFstArcs(cur_state, best_fst_) != 1) return false;
			const State *state_node = GetFstState(cur_state, best_fst_);
			const Arc arc = state_node->arc[0];
			tot_weight = Times(arc.arc_weight, tot_weight);
			if (arc.ilabel != 0) ilabel_seq.push_back(arc.ilabel);
			if (arc.olabel != 0) olabel_seq.push_back(arc.olabel);
			cur_state = arc.next_state;
		}
	}
}


void DeleteElems(Elem *list) 
{
	for (Elem *e = list, *e_tail; e != NULL; e = e_tail)
	{
		e_tail = e->tail;
		DeleteHashList(e);
	}
}

void ClearActiveTokens() 
{ // a cleanup routine, at utt end/begin
	for (int i = 0; i < active_toks_.size(); ++i) 
	{
		// Delete all tokens alive on this frame, and any forward
		// links they may have.
		for (Token *tok = active_toks_[i].toks; tok != NULL; ) 
		{
			DeleteForwardLinks(tok);
			Token *next_tok = tok->next;
			delete tok;
			num_toks_--;
			tok = next_tok;
		}
	}
	active_toks_.clear();

	assert(num_toks_ == 0);
}

void DeleteForwardLinks(Token *tok) 
{
	ForwardLink *l = tok->links, *m;
	while (l != NULL) 
	{
		m = l->next;
		delete l;
		l = m;
	}
	tok->links = NULL;
}

void DestroyDecode()
{
	DeleteElems(ClearHashList());
	ClearActiveTokens();
	cost_offsets_.clear();
	queue_.clear();
	tmp_array_.clear();
	final_costs_.clear();

	DeleteFstStates(raw_fst_);
	DeleteFstStates(best_fst_);
}