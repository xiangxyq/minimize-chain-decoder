/*
* The main feat extraction file
*
* author: xiangxyq
* E-mail: xiangxyq@foxmail.com
*/

#include <cstdlib>
#include <cfloat>
#include <cassert>
#include <cstring>
#include "compute-feat.h"
#include "../configs/feat-conf.h"
#include "mel-banks.h"
#include "srfft.h"
#include "../utils/asr-math.h"


static BaseFloat dct_matrix_rows[NUM_CEPS * NUM_BINS] = {0};
static BaseFloat lifter_coeffs[NUM_CEPS] = {0};
static BaseFloat window[PADDING_WINDOWS_SIZE] = {0};
static BaseFloat mel_energies[NUM_BINS] = {0};
static BaseFloat window_func[WINDOW_SIZE] = {0};
static BaseFloat power_spectrum[PADDING_WINDOWS_SIZE / 2 + 1] = {0};

//  pre-emphasis
static void Preemphasize(BaseFloat *waveform, int32 len);
// process window function
static void ProcessWindow(BaseFloat *window, int32 len);
// choose window function, defaut is provey
static void FeatureWindowFunction(BaseFloat *window_func);
// framing
static void ExtractWindow(const std::vector<BaseFloat> &wave, int32 f, BaseFloat *window);
// Discrete Cosine Transformation (dct)
static void ComputeDctMatrix(BaseFloat *M, int32 row, int32 col);
// Scaling factor on cepstra for HTK compatibility.
static void ComputeLifterCoeffs(BaseFloat *coeffs, int32 len);
// compute power spectrum
static void ComputePowerSpectrum(BaseFloat *waveform, int32 len);


void ComputeFeatInit()
{
	FeatureWindowFunction(window_func);
	SplitRadixFftInit();
	MelBanksInit();

	BaseFloat dct_matrix[NUM_BINS * NUM_BINS];
	memset(dct_matrix, 0, sizeof(BaseFloat) * NUM_BINS * NUM_BINS);
	ComputeDctMatrix(dct_matrix, NUM_BINS, NUM_BINS);
	memcpy((BaseFloat*)(dct_matrix_rows), (BaseFloat*)(dct_matrix), sizeof(BaseFloat) * (NUM_CEPS * NUM_BINS));

	ComputeLifterCoeffs(lifter_coeffs, NUM_CEPS);
}


void ComputeFeatures(const vector<BaseFloat> &waveform, vector<BaseFloat* > &feats)
{
	if (waveform.size() == 0)
		return;  //nothing to do

	int32 num_frames_total = waveform.size();
	int32 num_frames_new = (num_frames_total < WINDOW_SIZE) ? 0 : \
                           (1 + ((num_frames_total - WINDOW_SIZE) / WINDOW_SHIFT));

	feats.resize(num_frames_new, NULL);

	for (int32 frame = 0; frame < num_frames_new; ++frame)
	{
		ExtractWindow(waveform, frame, window);

		SplitRadixFftCompute(window);
		// Convert the FFT into a power spectrum.
		ComputePowerSpectrum(window, PADDING_WINDOWS_SIZE);
		memcpy((BaseFloat*)(power_spectrum), (BaseFloat*)(window), sizeof(BaseFloat) * (PADDING_WINDOWS_SIZE / 2 + 1));

		MelBanksCompute(power_spectrum, mel_energies);
		// avoid log of zero (which should be prevented anyway by dithering).
		ApplyFloor(mel_energies,FLT_EPSILON, NUM_BINS);  //numeric_limits<float>::epsilon()
		ApplyLog(mel_energies, NUM_BINS);  // take the log.

		// feature = dct_matrix_ * mel_energies [which now have log]
		BaseFloat *this_feature = (BaseFloat *)malloc(sizeof(BaseFloat) * (NUM_CEPS));
		AddMatVec(1.0f, dct_matrix_rows, NUM_CEPS, NUM_BINS, STRIDE, kNoTrans, mel_energies, NUM_BINS, 0.0f, this_feature, NUM_CEPS);
		VectorMulElements(this_feature, lifter_coeffs, NUM_CEPS);

		feats[frame] = this_feature;
	}
	//dropout last frames, no need process
}

int32 NumFramesReady(const vector<BaseFloat* > &feats)
{
	return feats.size();
}

bool IsLastFrame(const vector<BaseFloat* > &feats, const int32 frame)
{
	if (frame == NumFramesReady(feats) - 1)
		return true;
	else
		return false;
}

BaseFloat *GetFrame(const vector<BaseFloat* > &feats, const int32 frame)
{
	if (frame < feats.size())
		return feats[frame];

	return NULL;
}


static void FeatureWindowFunction(BaseFloat *window_func)
{
	double a = M_2PI / (WINDOW_SIZE - 1);

	for (int32 i = 0; i < WINDOW_SIZE; ++i)
	{
		double i_fl = (double)(i);
#ifdef WINDOW_TYPE_HANNING
		window_func[i] = 0.5f - 0.5f * cos(a * i_fl);
#endif
#ifdef WINDOW_TYPE_HAMMING
		window_func[i] = 0.54f - 0.46f * cos(a * i_fl);
#endif
#ifdef WINDOW_TYPE_POVEY
  // like hamming but goes to zero at edges.
		window_func[i] = pow(0.5f - 0.5f * cos(a * i_fl), 0.85f);
#endif
#ifdef WINDOW_TYPE_RECTANGULAR
		window_func[i] = 1.0;
#endif
#ifdef WINDOW_TYPE_BLACKMAN
		window_func[i] = BLACKMAN_COEFF - 0.5f * cos(a * i_fl) + (0.5f - BLACKMAN_COEFF) * cos(2 * a * i_fl);
#endif
	}
}

static void Preemphasize(BaseFloat *waveform, int32 len)
{
	for (int32 i = len - 1; i > 0; --i)
			waveform[i] -= PREEMPH_COEFF * waveform[i - 1];
	waveform[0] -= PREEMPH_COEFF * waveform[0];
}

static void ProcessWindow(BaseFloat *window, int32 len)
{
	BaseFloat window_sum = VectorSum(window, len);
	VectorAdd(window, len, -window_sum / WINDOW_SIZE);
	Preemphasize(window, len);

	VectorMulElements(window, window_func, WINDOW_SIZE);
}

static void ExtractWindow(const std::vector<BaseFloat> &wave, int32 f, BaseFloat *window)
{
	int32 num_samples = wave.size(), start_sample = FirstSampleOfFrame(f), end_sample = start_sample + WINDOW_SIZE;
	assert(start_sample >= 0 && end_sample <= num_samples);

	memset(window, 0, sizeof(BaseFloat) * PADDING_WINDOWS_SIZE);
	// copy window length size data to window buffer, the remaining (PADDING_WINDOWS_SIZE-WINDOW_SIZE) filled with zero
	if (start_sample >= 0 && end_sample <= wave.size())
		memcpy((BaseFloat*)(window), (const BaseFloat*)(&wave[0] + start_sample), sizeof(BaseFloat) * WINDOW_SIZE);
	ProcessWindow(window, WINDOW_SIZE);
}

static void ComputeDctMatrix(BaseFloat *M, int32 row, int32 col)
{
	BaseFloat normalizer = sqrt(1.0f / (BaseFloat)(col));  // normalizer for X_0.
	for (int32 j = 0; j < col; ++j) M[j] = normalizer;
	normalizer = sqrt(2.0 / (BaseFloat)(col));  // normalizer for other elements.
	for (int32 k = 1; k < row; ++k)
		for (int32 n = 0; n < col; ++n)
			M[k * col + n] = normalizer * cos((double)(M_PI) / col * (n + 0.5f) * k);
}

static void ComputeLifterCoeffs(BaseFloat *coeffs, int32 len)
{
	// Compute liftering coefficients (scaling on cepstral coeffs)
	// coeffs are numbered slightly differently from HTK: the zeroth
	// index is C0, which is not affected.
	for (int32 i = 0; i < len; ++i)
		coeffs[i] = 1.0f + 0.5f * CEPSTRAL_LIFTER * sin(M_PI * i / CEPSTRAL_LIFTER);
}


static void ComputePowerSpectrum(BaseFloat *waveform, int32 len)
{
	int32 dim = len;

	// no, letting it be non-power-of-two for now.
	// KALDI_ASSERT(dim > 0 && (dim & (dim-1) == 0));  // make sure a power of two.. actually my FFT code
	// does not require this (dan) but this is better in case we use different code [dan].

	// RealFft(waveform, true);  // true == forward (not inverse) FFT; makes no difference here,
	// as we just want power spectrum.

	// now we have in waveform, first half of complex spectrum
	// it's stored as [real0, realN/2-1, real1, im1, real2, im2, ...]
	int32 half_dim = dim / 2;
	BaseFloat first_energy = waveform[0] * waveform[0],
		last_energy = waveform[1] * waveform[1];  // handle this special case
	for (int32 i = 1; i < half_dim; ++i) {
		BaseFloat real = waveform[i * 2], im = waveform[i * 2 + 1];
		waveform[i] = real * real + im * im;
	}
	waveform[0] = first_energy;
	waveform[half_dim] = last_energy;  // Will actually never be used, and anyway
	// if the signal has been bandlimited sensibly this should be zero.
}

void FeatsFree(vector<BaseFloat* > &feats)
{
	for (int32 i = 0; i < feats.size(); ++i)
	{
		free(feats[i]);
		feats[i] = NULL;
	}
	feats.clear();
}

void FeatDestroy(vector<BaseFloat* > &feats)
{
	FreeMelBanksMemory();
	FreeSrfftMemory();
	FeatsFree(feats);
}
