/*
* Feature extraction header
*
* author: xiangxyq
* E-mail: xiangxyq@foxmail.com
*/

#ifndef COMPUTE_FEAT_H
#define COMPUTE_FEAT_H

#include <vector>
#include "../utils/asr-types.h"
using namespace std;

// feature extraction init itf
extern void ComputeFeatInit();

// compute mfcc features and saved in vector
extern void ComputeFeatures(const vector<BaseFloat> &waveform, vector<BaseFloat* > &feats);

// return total mfcc features' frames
extern int32 NumFramesReady(const vector<BaseFloat* > &feats);

//  last frame or not
extern bool IsLastFrame(const vector<BaseFloat* > &feats, const int32 frame);

// return current frame index feature
extern BaseFloat *GetFrame(const vector<BaseFloat* > &feats, const int32 frame);

// free current waves alloc memory
extern void FeatsFree(vector<BaseFloat* > &feats);

// free all feature extract memory, contains mfcc features, mel bank, srfft malloc memory
extern void FeatDestroy(vector<BaseFloat* > &feats);


#endif // COMPUTE_FEAT_H
