/*
* fft compute
*
* author: xiangxyq
* E-mail: xiangxyq@foxmail.com
*/

#include <iostream>
#include <cstdlib>
#include <cstring>
#include "srfft.h"
#include "../utils/asr-math.h"
#include "../configs/feat-conf.h"

using namespace std;

static int32 logn_;
static BaseFloat **tab;
static int32 *brseed;

static void ComputeTables();
static void ComputeRecursive(BaseFloat *xr, BaseFloat *xi, int32 logn);
static void BitReversePermute(BaseFloat *x, int32 logn);
static void Compute(BaseFloat *xr, BaseFloat *xi);

void SplitRadixFftInit()
{
	if ((COMPLEX_N & (COMPLEX_N - 1)) != 0 || COMPLEX_N <= 1)
	{
		cout << "SplitRadixComplexFft called with invalid number of points " << COMPLEX_N << endl;
		exit(1);
	}

	int32 complex_N = COMPLEX_N;
	logn_ = 0;
	while (complex_N > 1) {
		complex_N >>= 1;
		logn_++;
	}
	ComputeTables();
}

static void ComputeTables()
{
	int32    imax, lg2, i, j;
	int32     m, m2, m4, m8, nel, n;
	BaseFloat    *cn, *spcn, *smcn, *c3n, *spc3n, *smc3n;
	BaseFloat    ang, c, s;

	lg2 = logn_ >> 1;
	if (logn_ & 1) lg2++;
	brseed = (int32 *) malloc (sizeof(int32) * (1 << lg2)); // brseed不要忘记释放

	brseed[0] = 0;
	brseed[1] = 1;
	for (j = 2; j <= lg2; ++j) {
		imax = 1 << (j - 1);
		for (i = 0; i < imax; ++i) {
			brseed[i] <<= 1;
			brseed[i + imax] = brseed[i] + 1;
		}
	}

	if (logn_ < 4) {
		tab = NULL;
	}
	else {
		tab = (BaseFloat **) malloc (sizeof(BaseFloat *) * (logn_ - 3));

		for (i = logn_; i >= 4; --i) {
			/* Compute a few constants */
			m = 1 << i; m2 = m / 2; m4 = m2 / 2; m8 = m4 / 2;
			/* Allocate memory for tables */
			nel = m4 - 2;
			tab[i - 4] = (BaseFloat *)malloc(sizeof(BaseFloat) * (6 * nel));

			/* Initialize pointers */
			cn = tab[i - 4]; spcn = cn + nel;  smcn = spcn + nel;
			c3n = smcn + nel;  spc3n = c3n + nel; smc3n = spc3n + nel;

			/* Compute tables */
			for (n = 1; n < m4; n++) {
				if (n == m8) continue;
				ang = n * M_2PI / m;
				c = cos(ang); s = sin(ang);
				*cn++ = c; *spcn++ = -(s + c); *smcn++ = s - c;
				ang = 3 * n * M_2PI / m;
				c = cos(ang); s = sin(ang);
				*c3n++ = c; *spc3n++ = -(s + c); *smc3n++ = s - c;
			}
		}
	}
}

static void ComputeRecursive(BaseFloat *xr, BaseFloat *xi, int32 logn)
{
	int32 m, m2, m4, m8, nel, n;
	BaseFloat *xr1, *xr2, *xi1, *xi2;
	BaseFloat *cn = NULL, *spcn = NULL, *smcn = NULL, *c3n = NULL,
		*spc3n = NULL, *smc3n = NULL;
	BaseFloat   tmp1, tmp2;
	BaseFloat   sqhalf = M_SQRT1_2;

	/* Check range of logn */
	if (logn < 0)
		cout << "invalid logn " << logn << endl;

	/* Compute trivial cases */
	if (logn < 3) {
		if (logn == 2) {  /* length m = 4 */
			xr2 = xr + 2;
			xi2 = xi + 2;
			tmp1 = *xr + *xr2;
			*xr2 = *xr - *xr2;
			*xr = tmp1;
			tmp1 = *xi + *xi2;
			*xi2 = *xi - *xi2;
			*xi = tmp1;
			xr1 = xr + 1;
			xi1 = xi + 1;
			xr2++;
			xi2++;
			tmp1 = *xr1 + *xr2;
			*xr2 = *xr1 - *xr2;
			*xr1 = tmp1;
			tmp1 = *xi1 + *xi2;
			*xi2 = *xi1 - *xi2;
			*xi1 = tmp1;
			xr2 = xr + 1;
			xi2 = xi + 1;
			tmp1 = *xr + *xr2;
			*xr2 = *xr - *xr2;
			*xr = tmp1;
			tmp1 = *xi + *xi2;
			*xi2 = *xi - *xi2;
			*xi = tmp1;
			xr1 = xr + 2;
			xi1 = xi + 2;
			xr2 = xr + 3;
			xi2 = xi + 3;
			tmp1 = *xr1 + *xi2;
			tmp2 = *xi1 + *xr2;
			*xi1 = *xi1 - *xr2;
			*xr2 = *xr1 - *xi2;
			*xr1 = tmp1;
			*xi2 = tmp2;
			return;
		}
		else if (logn == 1) {   /* length m = 2 */
			xr2 = xr + 1;
			xi2 = xi + 1;
			tmp1 = *xr + *xr2;
			*xr2 = *xr - *xr2;
			*xr = tmp1;
			tmp1 = *xi + *xi2;
			*xi2 = *xi - *xi2;
			*xi = tmp1;
			return;
		}
		else if (logn == 0) return;   /* length m = 1 */
	}
	/* Compute a few constants */
	m = 1 << logn; m2 = m / 2; m4 = m2 / 2; m8 = m4 / 2;


	/* Step 1 */
	xr1 = xr; xr2 = xr1 + m2;
	xi1 = xi; xi2 = xi1 + m2;
	for (n = 0; n < m2; n++) {
		tmp1 = *xr1 + *xr2;
		*xr2 = *xr1 - *xr2;
		xr2++;
		*xr1++ = tmp1;
		tmp2 = *xi1 + *xi2;
		*xi2 = *xi1 - *xi2;
		xi2++;
		*xi1++ = tmp2;
	}

	/* Step 2 */
	xr1 = xr + m2; xr2 = xr1 + m4;
	xi1 = xi + m2; xi2 = xi1 + m4;
	for (n = 0; n < m4; n++) {
		tmp1 = *xr1 + *xi2;
		tmp2 = *xi1 + *xr2;
		*xi1 = *xi1 - *xr2;
		xi1++;
		*xr2++ = *xr1 - *xi2;
		*xr1++ = tmp1;
		*xi2++ = tmp2;
	}
	/* Steps 3 & 4 */
	xr1 = xr + m2; xr2 = xr1 + m4;
	xi1 = xi + m2; xi2 = xi1 + m4;
	if (logn >= 4) {
		nel = m4 - 2;
		cn = tab[logn - 4]; spcn = cn + nel;  smcn = spcn + nel;
		c3n = smcn + nel;  spc3n = c3n + nel; smc3n = spc3n + nel;
	}
	xr1++; xr2++; xi1++; xi2++;
	for (n = 1; n < m4; n++) {
		if (n == m8) {
			tmp1 = sqhalf * (*xr1 + *xi1);
			*xi1 = sqhalf * (*xi1 - *xr1);
			*xr1 = tmp1;
			tmp2 = sqhalf * (*xi2 - *xr2);
			*xi2 = -sqhalf * (*xr2 + *xi2);
			*xr2 = tmp2;
		}
		else {
			tmp2 = *cn++ * (*xr1 + *xi1);
			tmp1 = *spcn++ * *xr1 + tmp2;
			*xr1 = *smcn++ * *xi1 + tmp2;
			*xi1 = tmp1;
			tmp2 = *c3n++ * (*xr2 + *xi2);
			tmp1 = *spc3n++ * *xr2 + tmp2;
			*xr2 = *smc3n++ * *xi2 + tmp2;
			*xi2 = tmp1;
		}
		xr1++; xr2++; xi1++; xi2++;
	}

	/* Call ssrec again with half DFT length */
	ComputeRecursive(xr, xi, logn - 1);

	/* Call ssrec again twice with one quarter DFT length.
	   Constants have to be recomputed, because they are static! */
	ComputeRecursive(xr + m2, xi + m2, logn - 2);
	m4 = 3 * (m / 4);
	ComputeRecursive(xr + m4, xi + m4, logn - 2);
}

static void BitReversePermute(BaseFloat *x, int32 logn)
{
	int32      i, j, lg2, n;
	int32      off, fj, gno, *brp;
	BaseFloat    tmp, *xp, *xq;

	lg2 = logn >> 1;
	n = 1 << lg2;
	if (logn & 1) lg2++;

	/* Unshuffling loop */
	for (off = 1; off < n; off++) {
		fj = n * brseed[off]; i = off; j = fj;
		tmp = x[i]; x[i] = x[j]; x[j] = tmp;
		xp = &x[i];
		brp = &(brseed[1]);
		for (gno = 1; gno < brseed[off]; gno++) {
			xp += n;
			j = fj + *brp++;
			xq = x + j;
			tmp = *xp; *xp = *xq; *xq = tmp;
		}
	}
}

static void Compute(BaseFloat *xr, BaseFloat *xi)
{
	ComputeRecursive(xr, xi, logn_);

	if (logn_ > 1)
	{
		BitReversePermute(xr, logn_);
		BitReversePermute(xi, logn_);
	}
}

void SplitRadixFftCompute(BaseFloat *data)
{
	int32 N = REAL_N, N2 = REAL_N / 2;
	BaseFloat temp_buffer[COMPLEX_N];

	BaseFloat *temp_ptr = &(temp_buffer[0]);
	for (int32 i = 0; i < COMPLEX_N; i++)
	{
		data[i] = data[i * 2];  // put the real part in the first half of data.
		temp_ptr[i] = data[i * 2 + 1];  // put the imaginary part in temp_buffer.
	}

	// copy the imaginary part back to the second half of data.
	memcpy((BaseFloat*)(data + COMPLEX_N), (BaseFloat*)(temp_ptr), sizeof(BaseFloat) * COMPLEX_N);
	Compute(data, data + COMPLEX_N);

	memcpy((BaseFloat*)(temp_ptr), (BaseFloat*)(data + COMPLEX_N), sizeof(BaseFloat) * COMPLEX_N);

	for (int32 i = COMPLEX_N - 1; i > 0; --i)
	{
		// don't include 0,
  	// in case MatrixIndexT is unsigned, the loop would not terminate.
  	// Treat it as a special case.
		data[i * 2] = data[i];
		data[i * 2 + 1] = temp_ptr[i];
	}
	data[1] = temp_ptr[0];  // special case of i = 0.

	BaseFloat rootN_re, rootN_im;  // exp(-2pi/N), forward; exp(2pi/N), backward
	int32 forward_sign = -1;
	ComplexImExp((BaseFloat)(M_2PI / N * forward_sign), &rootN_re, &rootN_im);
	BaseFloat kN_re = -forward_sign, kN_im = 0.0;  // exp(-2pik/N), forward; exp(-2pik/N), backward
	// kN starts out as 1.0 for forward algorithm but -1.0 for backward.
	for (int32 k = 1; 2 * k <= N2; k++) {
		ComplexMul(rootN_re, rootN_im, &kN_re, &kN_im);

		BaseFloat Ck_re, Ck_im, Dk_re, Dk_im;
		// C_k = 1/2 (B_k + B_{N/2 - k}^*) :
		Ck_re = 0.5 * (data[2 * k] + data[N - 2 * k]);
		Ck_im = 0.5 * (data[2 * k + 1] - data[N - 2 * k + 1]);
		// re(D_k)= 1/2 (im(B_k) + im(B_{N/2-k})):
		Dk_re = 0.5 * (data[2 * k + 1] + data[N - 2 * k + 1]);
		// im(D_k) = -1/2 (re(B_k) - re(B_{N/2-k}))
		Dk_im = -0.5 * (data[2 * k] - data[N - 2 * k]);
		// A_k = C_k + 1^(k/N) D_k:
		data[2 * k] = Ck_re;  // A_k <-- C_k
		data[2 * k + 1] = Ck_im;
		// now A_k += D_k 1^(k/N)
		ComplexAddProduct(Dk_re, Dk_im, kN_re, kN_im, &(data[2 * k]), &(data[2 * k + 1]));

		int32 kdash = N2 - k;
		if (kdash != k) {
			// Next we handle the index k' = N/2 - k.  This is necessary
			// to do now, to avoid invalidating data that we will later need.
			// The quantities C_{k'} and D_{k'} are just the conjugates of C_k
			// and D_k, so the equations are simple modifications of the above,
			// replacing Ck_im and Dk_im with their negatives.
			data[2 * kdash] = Ck_re;  // A_k' <-- C_k'
			data[2 * kdash + 1] = -Ck_im;
			// now A_k' += D_k' 1^(k'/N)
			// We use 1^(k'/N) = 1^((N/2 - k) / N) = 1^(1/2) 1^(-k/N) = -1 * (1^(k/N))^*
			// so it's the same as 1^(k/N) but with the real part negated.
			ComplexAddProduct(Dk_re, -Dk_im, -kN_re, kN_im, &(data[2 * kdash]), &(data[2 * kdash + 1]));
		}
	}

	{// Now handle k = 0.
	 // In simple terms: after the complex fft, data[0] becomes the sum of real
	 // parts input[0], input[2]... and data[1] becomes the sum of imaginary
	 // pats input[1], input[3]...
	 // "zeroth" [A_0] is just the sum of input[0]+input[1]+input[2]..
	 // and "n2th" [A_{N/2}] is input[0]-input[1]+input[2]... .
		BaseFloat zeroth = data[0] + data[1],
			n2th = data[0] - data[1];
		data[0] = zeroth;
		data[1] = n2th;
	}
}

void FreeSrfftMemory()
{
	free(brseed);
	brseed = NULL;

	for (int i = 0; i < logn_ - 3; i++)
	{
		free(tab[i]);
		tab[i] = NULL;
	}
	free(tab);
	tab = NULL;
}
