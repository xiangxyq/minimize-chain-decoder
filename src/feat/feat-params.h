/*
* This file define some feature extraction parameters
*
* Notice Only Support MFCC
*
* author: xiangxyq
* E-mail: xiangxyq@foxmail.com
*/


#ifndef FEAT_PARAMS_H
#define FEAT_PARAMS_H

#include <cmath>
#include "../utils/asr-types.h"


// Define mfcc sample rate, frame shift, frame length and preemph coeff
#define SAMPLEING_RATE        (16000.0f)
#define FRAME_SHIFT_MS        (10.0f)
#define FRAME_LENGHT_MS       (25.0f)
#define PREEMPH_COEFF         (0.97f)

// Define mfcc window type, default is povey(WINDOW_TYPE_POVEY)
// if specified, you can choose below windows:
// 1. hamming(WINDOW_TYPE_HAMMING) 2. hanning(WINDOW_TYPE_HANNING)
// 3. rectangular((WINDOW_TYPE_RECTANGULAR)) 4. blackman(WINDOW_TYPE_BLACKMAN)
#define WINDOW_TYPE_POVEY
#ifdef WINDOW_TYPE_BLACKMAN
  #define BLACKMAN_COEFF        (0.42f)
#endif

// Define window shift steps, window size and padding window size
#if 0
static inline int32 RoundUpToNearestPowerOfTwo(int32 n)
{
  n--;
  n |= n >> 1;
  n |= n >> 2;
  n |= n >> 4;
  n |= n >> 8;
  n |= n >> 16;
  return n+1;
}
#endif

#define WINDOW_SHIFT          (static_cast<int32> ((SAMPLEING_RATE) * 0.001f * (FRAME_SHIFT_MS)))
#define WINDOW_SIZE           (static_cast<int32> ((SAMPLEING_RATE) * 0.001f * (FRAME_LENGHT_MS)))
//#define PADDING_WINDOWS_SIZE  ((RoundUpToNearestPowerOfTwo(WINDOW_SIZE)))
#define PADDING_WINDOWS_SIZE  (512)

// Define mfcc dims if dim is 13, you need to define macro MFCC_13_DIM
// default is 40 dims
#ifdef MFCC_13_DIM
  #define NUM_BINS               (23)
  #define NUM_CEPS               (13)
  #define CEPSTRAL_LIFTER        (22)
  #define STRIDE                 (24)
#else
  #define NUM_BINS               (40)
  #define NUM_CEPS               (40)
  #define CEPSTRAL_LIFTER        (22)
  #define STRIDE                 (40)
#endif

#define FirstSampleOfFrame(frame) ((frame) * (WINDOW_SHIFT))

// Define mel-banks params
#define MelScale(freq)             (1127.0f * logf(1.0f + (freq) / 700.0f))
#define InverseMelScale(mel_freq)  (700.0f * (expf((mel_freq) / 1127.0f) - 1.0f))

// Define fft params: number fft bins, nyquist freq, low freq, high freq and fft bin width
#define NUM_FFT_BINS           (static_cast<int32> (PADDING_WINDOWS_SIZE/2))
#define NYQUIST                (static_cast<int32> (0.5f * (SAMPLEING_RATE)))
#define LOW_FREQ               (20)
#define HIGH_FREQ              (static_cast<int32>((NYQUIST) - 400))
#define FFT_BIN_WIDTH          (static_cast<BaseFloat> (static_cast<BaseFloat> (SAMPLEING_RATE) / static_cast<BaseFloat> (PADDING_WINDOWS_SIZE)))

// Define mel params: mel low freq, mel high freq and vtln params
#define MEL_LOW_FREQ           (static_cast<BaseFloat> (MelScale(LOW_FREQ)))
#define MEL_HIGH_FREQ          (static_cast<BaseFloat> (MelScale(HIGH_FREQ)))
#define MEL_FREQ_DELTA         (static_cast<BaseFloat> (((MEL_HIGH_FREQ) - (MEL_LOW_FREQ)) / (static_cast<BaseFloat> (NUM_BINS) + 1)))
#define VTLN_WARP_FACTOR       (1.0f)
#define VTLN_LOW               (100)
#define VTLN_HIGH              ((NYQUIST) - 500)

// Define real and complex size
#define REAL_N                 (PADDING_WINDOWS_SIZE)
#define COMPLEX_N              (static_cast<int32>((PADDING_WINDOWS_SIZE) / 2))


#endif // FEAT_PARAMS_H
