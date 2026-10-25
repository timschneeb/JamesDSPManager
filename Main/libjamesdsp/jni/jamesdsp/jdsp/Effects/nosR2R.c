#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include "../jdsp_header.h"

/*
 * NOS R2R Simulator — port of foo_nos_r2r v2.2.123 by Roman Kuznetsov (RAA)
 *
 * Full parameters (with absurd ranges for demonstration):
 * 1. Target sampling rate (ZOH resampling for NOS simulation)
 * 2. Resistor tolerance (0-100%, absurd ranges for demo)
 * 3. Deviation growth (1-10, how deviation scales with bit significance)
 * 4. Harmony (Random/Even/Odd harmonic dominance)
 * 5. Serial number (PRNG seed for deviation pattern)
 * 6. Jitter amount (0.0-1.0, per-sample random variation, absurd at high values)
 * 7. Harmonics amount (0.0-1.0, 2nd/3rd harmonic distortion, absurd at high values)
 * 8. Invert distortion phase
 * 9. Bit depth ("bit + sign" convention: N bits = N-1 magnitude + sign)
 */

static uint64_t nosr2r_seed = 1;
static void nosr2r_srand(uint64_t s) { nosr2r_seed = s ? s : 1; }
static double nosr2r_rand()
{
	uint64_t x = nosr2r_seed;
	x ^= x << 13; x ^= x >> 7; x ^= x << 17;
	nosr2r_seed = x;
	return (double)(x >> 11) * (1.0 / 9007199254740992.0);
}

void NosR2RInit(NosR2R *nr, double fs)
{
	nr->srcRate = fs;
	nr->targetRate = fs;
	nr->ratio = 1.0;
	nr->bitDepth = 24;
	nr->resistorTolerance = 0.0;
	nr->deviationGrowth = 1;
	nr->harmony = 0;
	nr->serialNumber = 0;
	nr->jitterAmount = 0.0;
	nr->harmonicsAmount = 0.0;
	nr->invertPhase = 0;
	nr->fracAccum = 0.0;
	nr->lastSampleL = 0.0f;
	nr->lastSampleR = 0.0f;
	memset(nr->bitDeviation, 0, sizeof(nr->bitDeviation));
	memset(nr->bitDeviationR, 0, sizeof(nr->bitDeviationR));
	nr->deviationValid = 0;
}

static void NosR2RGenDeviations(NosR2R *nr)
{
	if (nr->resistorTolerance <= 0.0)
	{
		memset(nr->bitDeviation, 0, sizeof(nr->bitDeviation));
		memset(nr->bitDeviationR, 0, sizeof(nr->bitDeviationR));
		nr->deviationValid = 0;
		return;
	}
	nosr2r_srand((uint64_t)nr->serialNumber);
	double tol = nr->resistorTolerance;
	for (int ch = 0; ch < 2; ch++)
	{
		double *dev = (ch == 0) ? nr->bitDeviation : nr->bitDeviationR;
		for (int bit = 0; bit < 24; bit++)
		{
			double r = (nosr2r_rand() * 2.0 - 1.0) * tol;
			if (nr->deviationGrowth > 1)
			{
				double scale = pow((double)(bit + 1), (double)(nr->deviationGrowth - 1));
				r *= scale * 0.01;
			}
			if (nr->harmony == 1) { r *= (bit % 2 == 0) ? 1.5 : 0.3; }
			else if (nr->harmony == 2) { r *= (bit % 2 == 1) ? 1.5 : 0.3; }
			/* Clamp to prevent NaN/inf at absurd settings */
			if (r > 2.0) r = 2.0;
			if (r < -2.0) r = -2.0;
			dev[bit] = r;
		}
	}
	nr->deviationValid = 1;
}

void NosR2RSetParam(JamesDSPLib *jdsp, double targetRate, int bitDepth,
                     double resistorTolerance, int deviationGrowth,
                     int harmony, long serialNumber,
                     double jitterAmount, double harmonicsAmount, int invertPhase)
{
	NosR2R *nr = &jdsp->nosR2R;
	nr->srcRate = jdsp->fs;
	nr->targetRate = targetRate;
	if (nr->targetRate < nr->srcRate) nr->targetRate = nr->srcRate;
	nr->ratio = nr->targetRate / nr->srcRate;
	if (nr->ratio < 1.0) nr->ratio = 1.0;
	if (bitDepth < 1) bitDepth = 1;
	if (bitDepth > 24) bitDepth = 24;
	nr->bitDepth = bitDepth;
	/* Absurd tolerance: allow up to 100% (1.0) */
	if (resistorTolerance < 0.0) resistorTolerance = 0.0;
	if (resistorTolerance > 1.0) resistorTolerance = 1.0;
	nr->resistorTolerance = resistorTolerance;
	if (deviationGrowth < 1) deviationGrowth = 1;
	if (deviationGrowth > 10) deviationGrowth = 10;
	nr->deviationGrowth = deviationGrowth;
	nr->harmony = harmony;
	nr->serialNumber = serialNumber;
	/* Absurd jitter/harmonics: allow up to 1.0 */
	if (jitterAmount < 0.0) jitterAmount = 0.0;
	if (jitterAmount > 1.0) jitterAmount = 1.0;
	nr->jitterAmount = jitterAmount;
	if (harmonicsAmount < 0.0) harmonicsAmount = 0.0;
	if (harmonicsAmount > 1.0) harmonicsAmount = 1.0;
	nr->harmonicsAmount = harmonicsAmount;
	nr->invertPhase = invertPhase;
	NosR2RGenDeviations(nr);
}

void NosR2REnable(JamesDSPLib *jdsp)
{
	/* Do not call NosR2RInit here — it would wipe parameters set by NosR2RSetParam.
	 * Just reset runtime state and mark enabled. */
	NosR2R *nr = &jdsp->nosR2R;
	nr->fracAccum = 0.0;
	nr->lastSampleL = 0.0f;
	nr->lastSampleR = 0.0f;
	jdsp->nosR2REnabled = 1;
}

void NosR2RDisable(JamesDSPLib *jdsp)
{
	jdsp->nosR2REnabled = 0;
}

static inline float NosR2RApplyDeviation(float sample, const double *dev, int bitDepth,
                                          double jitterAmount, double harmonicsAmount,
                                          int invertPhase, int harmony)
{
	if (bitDepth <= 0) return 0.0f;
	float distortion = 0.0f;
	int magBits = bitDepth - 1;
	if (magBits > 23) magBits = 23;
	if (magBits < 0) magBits = 0;
	if (magBits > 0)
	{
		float sign = (sample >= 0.0f) ? 1.0f : -1.0f;
		float mag = fabsf(sample);
		float scale = (float)((1 << magBits) - 1);
		if (scale < 1.0f) scale = 1.0f;
		int intval = (int)(mag * scale + 0.5f);
		int maxVal = (1 << magBits) - 1;
		if (intval > maxVal) intval = maxVal;
		double reconstructed = 0.0;
		for (int bit = 0; bit < magBits; bit++)
		{
			if (intval & (1 << bit))
			{
				double normContrib = (double)(1 << bit) / (double)((1 << magBits) - 1);
				reconstructed += normContrib * (1.0 + dev[bit]);
			}
		}
		/* Jitter: scale from subtle (0.0001) to absurd (0.1) */
		if (jitterAmount > 0.0)
			reconstructed += (nosr2r_rand() * 2.0 - 1.0) * jitterAmount * 0.1;
		float result = (float)(sign * reconstructed);
		distortion = result - sample;
	}
	/* Harmonics: scale from subtle to absurd */
	if (harmonicsAmount > 0.0)
	{
		if (harmony == 1)
		{
			/* Even: 2nd harmonic, up to 1.5x at full */
			float h2 = sample * sample * (float)(harmonicsAmount * 1.5);
			if (invertPhase) h2 = -h2;
			distortion += h2;
		}
		else
		{
			/* Odd/Random: 3rd harmonic, up to 1.0x at full */
			float h3 = sample * sample * sample * (float)harmonicsAmount;
			if (invertPhase) h3 = -h3;
			distortion += h3;
		}
	}
	if (invertPhase && harmonicsAmount <= 0.0) distortion = -distortion;
	return sample + distortion;
}

void NosR2RProcess(JamesDSPLib *jdsp, size_t n)
{
	NosR2R *nr = &jdsp->nosR2R;
	float *bufL = jdsp->tmpBuffer[0];
	float *bufR = jdsp->tmpBuffer[1];
	if (nr->ratio <= 1.0 && nr->bitDepth >= 24 && nr->resistorTolerance <= 0.0
	    && nr->jitterAmount <= 0.0 && nr->harmonicsAmount <= 0.0)
		return;
	int magBits = nr->bitDepth - 1;
	if (magBits > 23) magBits = 23;
	for (size_t i = 0; i < n; i++)
	{
		float sL = bufL[i], sR = bufR[i];
		/* 1. Zero-Order Hold */
		if (nr->ratio > 1.0)
		{
			nr->fracAccum += 1.0;
			if (nr->fracAccum >= nr->ratio)
			{
				nr->fracAccum -= nr->ratio;
				nr->lastSampleL = sL;
				nr->lastSampleR = sR;
			}
			else { sL = nr->lastSampleL; sR = nr->lastSampleR; }
		}
		/* 2. Bit-depth truncation */
		if (magBits >= 0 && magBits < 23)
		{
			float signL = (sL >= 0.0f) ? 1.0f : -1.0f;
			float signR = (sR >= 0.0f) ? 1.0f : -1.0f;
			float magL = fabsf(sL), magR = fabsf(sR);
			float scale = (float)((1 << magBits) - 1);
			if (scale < 1.0f) { sL = 0.0f; sR = 0.0f; }
			else
			{
				int intL = (int)(magL * scale + 0.5f);
				int intR = (int)(magR * scale + 0.5f);
				int maxVal = (1 << magBits) - 1;
				if (intL > maxVal) intL = maxVal;
				if (intR > maxVal) intR = maxVal;
				sL = signL * (float)intL / scale;
				sR = signR * (float)intR / scale;
			}
		}
		/* 3. Resistor deviation + harmonics */
		if (nr->deviationValid)
		{
			sL = NosR2RApplyDeviation(sL, nr->bitDeviation, nr->bitDepth,
			                           nr->jitterAmount, nr->harmonicsAmount, nr->invertPhase, nr->harmony);
			sR = NosR2RApplyDeviation(sR, nr->bitDeviationR, nr->bitDepth,
			                           nr->jitterAmount, nr->harmonicsAmount, nr->invertPhase, nr->harmony);
		}
		else if (nr->harmonicsAmount > 0.0)
		{
			if (nr->harmony == 1)
			{
				float h2 = sL * sL * (float)(nr->harmonicsAmount * 1.5); if (nr->invertPhase) h2 = -h2; sL += h2;
				float h2r = sR * sR * (float)(nr->harmonicsAmount * 1.5); if (nr->invertPhase) h2r = -h2r; sR += h2r;
			}
			else
			{
				float h3 = sL * sL * sL * (float)nr->harmonicsAmount; if (nr->invertPhase) h3 = -h3; sL += h3;
				float h3r = sR * sR * sR * (float)nr->harmonicsAmount; if (nr->invertPhase) h3r = -h3r; sR += h3r;
			}
		}
		/* Clamp output to prevent NaN/inf at absurd settings */
		if (sL > 2.0f) sL = 2.0f;
		if (sL < -2.0f) sL = -2.0f;
		if (sR > 2.0f) sR = 2.0f;
		if (sR < -2.0f) sR = -2.0f;
		bufL[i] = sL;
		bufR[i] = sR;
	}
}
