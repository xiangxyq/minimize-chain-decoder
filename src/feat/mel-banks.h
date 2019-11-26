/*
* mel filter banks
*
* author: xiangxyq
* E-mail: xiangxyq@foxmail.com
*/

#ifndef MEL_BANKS_H
#define MEL_BANKS_H

#include "../utils/asr-types.h"

// mel bins struct
typedef struct
{
	int32 offset;
	BaseFloat *data;
	int32 len;
} MelBins;

// mel filter banks init
extern void MelBanksInit();

// compute mel_energies
extern void MelBanksCompute(BaseFloat *power_spectrum, BaseFloat *mel_energies_out);

// free mel filter banks alloc memory
extern void FreeMelBanksMemory();


#endif // MEL_BANKS_H
