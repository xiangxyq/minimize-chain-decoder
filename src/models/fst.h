/*
* This File hclg traversing operations
*
* author: xiangxyq
* E-mail: xiangxyq@foxmail.com
*/

#ifndef FST_H
#define FST_H

#include <vector>
#include "../utils/asr-types.h"

#define kNoStateId (-1)
typedef int32 StateId;
typedef BaseFloat Weight;

typedef struct
{
	int32 ilabel;
	int32 olabel;
	Weight arc_weight;
	StateId next_state;
} Arc;

typedef struct
{
	Weight weight;  /* e.g. , final weight */
	int32 num_arc;
	int32 niepsilons;
	int32 noepsilons;
	Arc *arc;
} State;

typedef enum
{
    HCLG_FST,
    RAW_FST,
    BEST_FST
} FstType;

typedef struct
{
    FstType type;
    StateId fst_start_;
    int32 num_states_;
    std::vector<State *> fst_states_;
} Fst_Strcut_Def;

extern void InitFst();

extern Fst_Strcut_Def *GetHclg();

extern Fst_Strcut_Def *GetRawFst();

extern Fst_Strcut_Def *GetBestFst();

extern State *GetFstState(StateId s, Fst_Strcut_Def *fst);

extern void SetFstStartState(StateId s, Fst_Strcut_Def *fst);

extern StateId GetFstStartState(Fst_Strcut_Def *fst);

extern void SetFstFinalWeight(StateId s, Weight weight, Fst_Strcut_Def *fst);

extern Weight GetFstFinalWeight(StateId s, Fst_Strcut_Def *fst);

extern int32 GetFstNumInputEpsilons(StateId s, Fst_Strcut_Def *fst);

extern int32 GetFstNumStates(Fst_Strcut_Def *fst);

extern int32 AddFstState(Fst_Strcut_Def *fst);

extern void AddFstArc(StateId s, Arc arc, Fst_Strcut_Def *fst);

extern int32 GetNumFstArcs(StateId s, Fst_Strcut_Def *fst);

//Tropical SemiRing , and multiply : plus
extern Weight Times(BaseFloat a, BaseFloat b);

extern void DeleteFstStates(Fst_Strcut_Def *fst);


#endif // FST_H