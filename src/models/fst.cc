/*
* This File hclg traversing operations
*
* author: xiangxyq
* E-mail: xiangxyq@foxmail.com
*/

#include <vector>
#include <cassert>
#include <math.h>
#include <cstdlib>
#include "hclg.h"
#include "fst.h"

static Fst_Strcut_Def hclg_, raw_fst_, best_fst_;
static State fst_[] = HCLG_GRAPH;

void InitFst()
{
    hclg_.type = HCLG_FST;
    hclg_.fst_start_ = START_STATE;
    hclg_.num_states_ = NUM_STATE_NODE;

    raw_fst_.type = RAW_FST;
    raw_fst_.fst_start_ = kNoStateId;
    raw_fst_.num_states_ = 0;

    best_fst_.type = BEST_FST;
    best_fst_.fst_start_ = kNoStateId;
    best_fst_.num_states_ = 0;
}

Fst_Strcut_Def *GetHclg()
{
	return &hclg_;
}

Fst_Strcut_Def *GetRawFst()
{
	return &raw_fst_;
}

Fst_Strcut_Def *GetBestFst()
{
	return &best_fst_;
}

State *GetFstState(StateId s, Fst_Strcut_Def *fst)
{
    assert(s >= fst->fst_start_ && s < fst->num_states_);

    if (fst->type == HCLG_FST)
	    return &(fst_[s]);
    else
        return fst->fst_states_[s];
}

void SetFstStartState(StateId s, Fst_Strcut_Def *fst)
{
    fst->fst_start_ = s;
}

StateId GetFstStartState(Fst_Strcut_Def *fst) 
{
	return fst->fst_start_;
}

void SetFstFinalWeight(StateId s, Weight weight, Fst_Strcut_Def *fst)
{
	fst->fst_states_[s]->weight = weight;
}

Weight GetFstFinalWeight(StateId s, Fst_Strcut_Def *fst) 
{
	State state = *(GetFstState(s, fst));
	return state.weight;
}

int32 GetFstNumInputEpsilons(StateId s, Fst_Strcut_Def *fst)
{
    State state = *(GetFstState(s, fst));
	return state.niepsilons;
}

int32 GetFstNumStates(Fst_Strcut_Def *fst)
{
	return fst->num_states_;
}

int32 AddFstState(Fst_Strcut_Def *fst)
{
	State *fst_state= (State *)calloc(1, sizeof(State));
	fst->fst_states_.push_back(fst_state);
	fst->num_states_ = fst->fst_states_.size();
	
	return fst->fst_states_.size() - 1;
}

void AddFstArc(StateId s, Arc arc, Fst_Strcut_Def *fst)
{
	fst->fst_states_[s]->num_arc++;
	int index = fst->fst_states_[s]->num_arc - 1;
	if (0 == index)
		fst->fst_states_[s]->arc = (Arc *)malloc(sizeof(Arc) * fst->fst_states_[s]->num_arc);
	else 
		fst->fst_states_[s]->arc = (Arc *)realloc(fst->fst_states_[s]->arc, sizeof(Arc) * fst->fst_states_[s]->num_arc);

	fst->fst_states_[s]->arc[index].ilabel = arc.ilabel;
	fst->fst_states_[s]->arc[index].olabel = arc.olabel;
	fst->fst_states_[s]->arc[index].arc_weight = arc.arc_weight;
	fst->fst_states_[s]->arc[index].next_state = arc.next_state;

	if (0 == fst->fst_states_[s]->arc[index].ilabel)
		fst->fst_states_[s]->niepsilons++;

	if (0 == fst->fst_states_[s]->arc[index].olabel)
		fst->fst_states_[s]->noepsilons++;

	fst->fst_states_[s]->weight = INFINITY;
}

int32 GetNumFstArcs(StateId s, Fst_Strcut_Def *fst)
{
	return fst->fst_states_[s]->num_arc;
}

Weight Times(BaseFloat a, BaseFloat b)
{
	if (1 == isinf(a)) 
		return a;
	else if (1 == isinf(b))  
		return b;
	else 
		return a + b;
}

void DeleteFstStates(Fst_Strcut_Def *fst)
{
	for (int i = 0; i < fst->num_states_; ++i)
	{
		for (int j = 0; j < fst->fst_states_[i]->num_arc; ++j)
		{
			free(fst->fst_states_[i]->arc);
			fst->fst_states_[i]->arc = NULL;
		}
		free(fst->fst_states_[i]);
		fst->fst_states_[i] = NULL;
	}
	fst->fst_start_ = kNoStateId;
	fst->fst_states_.clear();
}
