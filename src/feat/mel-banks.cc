/*
* mel filter banks
*
* author: xiangxyq
* E-mail: xiangxyq@foxmail.com
*/

#include <cstdlib>
#include <cstdlib>
#include <cstring>
#include "mel-banks.h"
#include "../utils/asr-math.h"
#include "../configs/feat-conf.h"

static MelBins mel_bins[NUM_BINS] = {0};
static BaseFloat center_freqs[NUM_BINS] = {0};

void MelBanksInit()
{
	for (int32 bin = 0; bin < NUM_BINS; bin++)
	{
		BaseFloat left_mel = MEL_LOW_FREQ + bin * MEL_FREQ_DELTA,
			center_mel = MEL_LOW_FREQ + (bin + 1) * MEL_FREQ_DELTA,
			right_mel = MEL_LOW_FREQ + (bin + 2) * MEL_FREQ_DELTA;

		center_freqs[bin] = InverseMelScale(center_mel);

		BaseFloat this_bin[NUM_FFT_BINS] = {0};
		int32 first_index = -1, last_index = -1;
		for (int32 i = 0; i < NUM_FFT_BINS; ++i)
		{
			BaseFloat freq = (FFT_BIN_WIDTH * i);  // Center frequency of this fft bin.
			BaseFloat mel = MelScale(freq);
			if (mel > left_mel && mel < right_mel)
			{
				BaseFloat weight;
				if (mel <= center_mel)
					weight = (mel - left_mel) / (center_mel - left_mel);
				else
					weight = (right_mel - mel) / (right_mel - center_mel);
				this_bin[i] = weight;
				if (first_index == -1)
					first_index = i;
				last_index = i;
			}
		}
		mel_bins[bin].offset = first_index;
		int32 size = last_index + 1 - first_index;
		mel_bins[bin].len = size;
		mel_bins[bin].data = (BaseFloat *)malloc(size * sizeof(BaseFloat));

		memcpy(mel_bins[bin].data, this_bin + first_index, sizeof(BaseFloat) * size);
	}
}

// "power_spectrum" contains fft energies.
void MelBanksCompute(BaseFloat *power_spectrum, BaseFloat *mel_energies_out)
{
	//mel_energies_out size must be NUM_BINS
	for (int32 i = 0; i < NUM_BINS; ++i)
	{
		int32 offset = mel_bins[i].offset;
		int32 mel_bins_dim = mel_bins[i].len;

		BaseFloat v[mel_bins_dim];
		memset(v, 0, sizeof(BaseFloat) * mel_bins_dim);
		BaseFloat power_spectrum_offset[mel_bins_dim];
		memset(power_spectrum_offset, 0, sizeof(BaseFloat) * mel_bins_dim);

		memcpy((BaseFloat*)(power_spectrum_offset), (BaseFloat*)(power_spectrum + offset), sizeof(BaseFloat) * mel_bins_dim);
		memcpy((BaseFloat*)(v), (BaseFloat*)(mel_bins[i].data), sizeof(BaseFloat) * mel_bins_dim);
		BaseFloat energy = VecVec(v, power_spectrum_offset, mel_bins_dim);

		mel_energies_out[i] = energy;
	}
}

void FreeMelBanksMemory()
{
	for (int32 i = 0; i < NUM_BINS; ++i)
	{
		free(mel_bins[i].data);
		mel_bins[i].data = NULL;
	}
}
