#include <math.h>
#include <string.h>
#include "../jdsp_header.h"

static const float SPECTRUM_HARMONICS[10] = {
	0.02f, 0.0f, 0.02f, 0.0f, 0.02f,
	0.0f, 0.02f, 0.0f, 0.02f, 0.0f,
};
static const float SPECTRUM_DEFAULT_Q = 0.717f;
static const int SPECTRUM_DEFAULT_LP_OFFSET_HZ = 2000;
static const float SPECTRUM_DEFAULT_WET_MIX = 1.0f;
static const float SPECTRUM_DEFAULT_POST_GAIN = 1.0f;
static const int SPECTRUM_DEFAULT_SAMPLE_RATE = 44100;
static const int SPECTRUM_REFERENCE_NYQUIST_MARGIN_HZ = 100;
static const float SPECTRUM_MIN_CUTOFF_HZ = 20.0f;

static inline int spectrumClampSamplingRate(int samplingRate)
{
	return samplingRate > 1 ? samplingRate : SPECTRUM_DEFAULT_SAMPLE_RATE;
}

static inline float spectrumClampFilterFrequency(float freq, int samplingRate)
{
	const float minFreq = SPECTRUM_MIN_CUTOFF_HZ;
	float maxFreq = (float)samplingRate * 0.5f - 1.0f;
	if (maxFreq < minFreq)
		maxFreq = minFreq;
	if (!isfinite(freq))
		return minFreq;
	if (freq < minFreq)
		return minFreq;
	if (freq > maxFreq)
		return maxFreq;
	return freq;
}

static inline float spectrumSanitizeFloat(float value, float fallback)
{
	return isfinite(value) ? value : fallback;
}

static inline float spectrumLinearFromDb(float db)
{
	return powf(10.0f, db / 20.0f);
}

static inline double spectrumClampSafety(double sample)
{
	if (sample > 1.0)
		return 1.0;
	if (sample < -1.0)
		return -1.0;
	return sample;
}

static inline void spectrumBiquadClear(SpectrumExtBiquad *bq)
{
	bq->x1 = 0.0;
	bq->x2 = 0.0;
	bq->y1 = 0.0;
	bq->y2 = 0.0;
}

static inline double spectrumBiquadProcessSample(SpectrumExtBiquad *bq, double sample)
{
	double out = sample * bq->b0 + bq->x1 * bq->b1 + bq->x2 * bq->b2 + bq->y1 * bq->a1 + bq->y2 * bq->a2;
	bq->x2 = bq->x1;
	bq->x1 = sample;
	bq->y2 = bq->y1;
	bq->y1 = out;
	return out;
}

static void spectrumBiquadSet(SpectrumExtBiquad *bq, int samplingRate, float freq, float q, char highPass)
{
	samplingRate = spectrumClampSamplingRate(samplingRate);
	freq = spectrumClampFilterFrequency(freq, samplingRate);
	q = spectrumSanitizeFloat(q, SPECTRUM_DEFAULT_Q);
	if (q <= 0.0f)
		q = SPECTRUM_DEFAULT_Q;

	const double pi = 3.14159265358979323846;
	double omega = 2.0 * pi * (double)freq / (double)samplingRate;
	double sinOmega = sin(omega);
	double cosOmega = cos(omega);
	double y = sinOmega / (2.0 * (double)q);
	double a0 = 1.0 + y;
	double a1 = -2.0 * cosOmega;
	double a2 = 1.0 - y;
	double b0;
	double b1;
	double b2;

	if (highPass)
	{
		b0 = (1.0 + cosOmega) / 2.0;
		b1 = -(1.0 + cosOmega);
		b2 = (1.0 + cosOmega) / 2.0;
	}
	else
	{
		b0 = (1.0 - cosOmega) / 2.0;
		b1 = 1.0 - cosOmega;
		b2 = (1.0 - cosOmega) / 2.0;
	}

	bq->a1 = -(a1 / a0);
	bq->a2 = -(a2 / a0);
	bq->b0 = b0 / a0;
	bq->b1 = b1 / a0;
	bq->b2 = b2 / a0;
	spectrumBiquadClear(bq);
}

static inline void spectrumHarmonicResetState(SpectrumExtHarmonic *harmonic)
{
	harmonic->prevOut = 0.0;
	harmonic->prevLast = 0.0;
	harmonic->sampleCounter = 0;
}

static void spectrumHarmonicUpdateCoeffs(SpectrumExtHarmonic *harmonic, const float in[10])
{
	float unkarr1[11] = { 0.0f };
	float unkarr2[11] = { 0.0f };
	float biggestCoeffVal = 0.0f;
	float absCoeffSum = 0.0f;

	for (uint32_t i = 0; i < 10; i++)
	{
		float a = fabsf(in[i]);
		absCoeffSum += a;
		if (a > biggestCoeffVal)
			biggestCoeffVal = a;
	}
	harmonic->biggestCoeff = (uint32_t)(biggestCoeffVal * 10000.0f);

	memcpy(unkarr1 + 1, in, 10 * sizeof(float));

	if (absCoeffSum > 1.0f)
	{
		float scale = 1.0f / absCoeffSum;
		for (uint32_t i = 1; i < 11; i++)
			unkarr1[i] *= scale;
	}

	memset(harmonic->coeffs, 0, sizeof(harmonic->coeffs));
	harmonic->coeffs[10] = unkarr1[10];

	for (uint32_t i = 2; i < 11; i++)
	{
		for (uint32_t j = 0; j < i; j++)
		{
			float tmp = unkarr2[i - j];
			unkarr2[i - j] = (float)harmonic->coeffs[i - j];
			harmonic->coeffs[i - j] = harmonic->coeffs[i - j - 1] * 2.0f - tmp;
		}
		float tmp = unkarr1[10 - i + 1] - unkarr2[0];
		unkarr2[0] = (float)harmonic->coeffs[0];
		harmonic->coeffs[0] = tmp;
	}

	for (uint32_t i = 1; i < 11; i++)
		harmonic->coeffs[10 - i + 1] = harmonic->coeffs[10 - i] - unkarr2[10 - i + 1];

	harmonic->coeffs[0] = unkarr1[0] / 2.0f - unkarr2[0];
}

static inline void spectrumHarmonicSet(SpectrumExtHarmonic *harmonic, const float harmonics[10])
{
	spectrumHarmonicUpdateCoeffs(harmonic, harmonics);
	spectrumHarmonicResetState(harmonic);
}

static inline double spectrumHarmonicProcessSample(SpectrumExtHarmonic *harmonic, double sample)
{
	double poly = harmonic->coeffs[10];
	for (int i = 9; i >= 0; i--)
		poly = harmonic->coeffs[i] + sample * poly;

	double prevOut = (poly + harmonic->prevOut * 0.999) - harmonic->prevLast;
	harmonic->prevLast = poly;
	harmonic->prevOut = prevOut;

	if (harmonic->sampleCounter < harmonic->biggestCoeff)
	{
		harmonic->sampleCounter++;
		return 0.0;
	}

	return prevOut;
}

static inline int spectrumClampReferenceFrequency(const SpectrumExtension *ext, int referenceFreq)
{
	int maxReferenceFreq = spectrumClampSamplingRate(ext->samplingRate) / 2 - SPECTRUM_REFERENCE_NYQUIST_MARGIN_HZ;
	if (referenceFreq > maxReferenceFreq)
		return maxReferenceFreq;
	return referenceFreq;
}

static inline int spectrumClampLpOffsetHz(const SpectrumExtension *ext, int lpOffsetHz)
{
	int maxOffset = spectrumClampSamplingRate(ext->samplingRate) / 2 - (int)SPECTRUM_MIN_CUTOFF_HZ;
	if (maxOffset < 0)
		maxOffset = 0;
	if (lpOffsetHz < 0)
		return 0;
	if (lpOffsetHz > maxOffset)
		return maxOffset;
	return lpOffsetHz;
}

static void spectrumResetLocked(SpectrumExtension *ext)
{
	ext->samplingRate = spectrumClampSamplingRate(ext->samplingRate);
	ext->referenceFreq = spectrumClampReferenceFrequency(ext, ext->referenceFreq);
	ext->lpOffsetHz = spectrumClampLpOffsetHz(ext, ext->lpOffsetHz);

	float hpQ = ext->hpQ > 0.0f ? ext->hpQ : SPECTRUM_DEFAULT_Q;
	float lpQ = ext->lpQ > 0.0f ? ext->lpQ : SPECTRUM_DEFAULT_Q;
	float lowPassFreq = spectrumClampFilterFrequency(
		(float)ext->samplingRate / 2.0f - (float)ext->lpOffsetHz,
		ext->samplingRate
	);
	for (int ch = 0; ch < 2; ch++)
	{
		spectrumBiquadSet(&ext->highpass[ch], ext->samplingRate, (float)ext->referenceFreq, hpQ, 1);
		spectrumBiquadSet(&ext->lowpass[ch], ext->samplingRate, lowPassFreq, lpQ, 0);
		spectrumHarmonicSet(&ext->harmonics[ch], ext->harmonicSeed);
	}
}

static inline void spectrumSetSamplingRateLocked(SpectrumExtension *ext, int samplingRate)
{
	samplingRate = spectrumClampSamplingRate(samplingRate);
	if (ext->samplingRate == samplingRate)
		return;

	ext->samplingRate = samplingRate;
	ext->referenceFreq = spectrumClampReferenceFrequency(ext, ext->referenceFreq);
	ext->lpOffsetHz = spectrumClampLpOffsetHz(ext, ext->lpOffsetHz);
	spectrumResetLocked(ext);
}

static inline void spectrumSetReferenceFrequencyLocked(SpectrumExtension *ext, int referenceFreq)
{
	ext->referenceFreq = spectrumClampReferenceFrequency(ext, referenceFreq);
	spectrumResetLocked(ext);
}

void SpectrumExtensionConstructor(JamesDSPLib *jdsp)
{
	SpectrumExtension *ext = &jdsp->spectrumExt;
	memset(ext, 0, sizeof(*ext));
	ext->samplingRate = 44100;
	ext->referenceFreq = 7600;
	ext->enabled = 0;
	ext->exciter = 0.0f;
	ext->wetMix = SPECTRUM_DEFAULT_WET_MIX;
	ext->postGain = SPECTRUM_DEFAULT_POST_GAIN;
	ext->safetyEnabled = 0;
	ext->hpQ = SPECTRUM_DEFAULT_Q;
	ext->lpQ = SPECTRUM_DEFAULT_Q;
	ext->lpOffsetHz = SPECTRUM_DEFAULT_LP_OFFSET_HZ;
	memcpy(ext->harmonicSeed, SPECTRUM_HARMONICS, sizeof(SPECTRUM_HARMONICS));
	spectrumResetLocked(ext);
	if ((int)jdsp->fs != ext->samplingRate)
		spectrumSetSamplingRateLocked(ext, (int)jdsp->fs);
}

void SpectrumExtensionSetSamplingRate(JamesDSPLib *jdsp, int samplingRate)
{
	jdsp_lock(jdsp);
	spectrumSetSamplingRateLocked(&jdsp->spectrumExt, samplingRate);
	jdsp_unlock(jdsp);
}

void SpectrumExtensionSetReferenceFrequency(JamesDSPLib *jdsp, int referenceFreq)
{
	jdsp_lock(jdsp);
	spectrumSetReferenceFrequencyLocked(&jdsp->spectrumExt, referenceFreq);
	jdsp_unlock(jdsp);
}

void SpectrumExtensionSetExciter(JamesDSPLib *jdsp, float exciter)
{
	jdsp_lock(jdsp);
	jdsp->spectrumExt.exciter = spectrumSanitizeFloat(exciter, 0.0f);
	jdsp_unlock(jdsp);
}

void SpectrumExtensionSetParam(JamesDSPLib *jdsp, float strengthLinear, int referenceFreq, float wetMix, float postGainDb, char safetyEnabled, float hpQ, float lpQ, int lpOffsetHz, const double harmonics[10])
{
	jdsp_lock(jdsp);
	SpectrumExtension *ext = &jdsp->spectrumExt;
	ext->referenceFreq = spectrumClampReferenceFrequency(ext, referenceFreq);
	ext->exciter = spectrumSanitizeFloat(strengthLinear, 0.0f);
	ext->wetMix = spectrumSanitizeFloat(wetMix, SPECTRUM_DEFAULT_WET_MIX);
	if (ext->wetMix < 0.0f)
		ext->wetMix = 0.0f;
	if (ext->wetMix > 1.0f)
		ext->wetMix = 1.0f;
	postGainDb = spectrumSanitizeFloat(postGainDb, 0.0f);
	if (postGainDb < -24.0f)
		postGainDb = -24.0f;
	if (postGainDb > 24.0f)
		postGainDb = 24.0f;
	ext->postGain = spectrumLinearFromDb(postGainDb);
	ext->safetyEnabled = safetyEnabled ? 1 : 0;
	ext->hpQ = spectrumSanitizeFloat(hpQ, SPECTRUM_DEFAULT_Q);
	if (ext->hpQ <= 0.0f)
		ext->hpQ = SPECTRUM_DEFAULT_Q;
	ext->lpQ = spectrumSanitizeFloat(lpQ, SPECTRUM_DEFAULT_Q);
	if (ext->lpQ <= 0.0f)
		ext->lpQ = SPECTRUM_DEFAULT_Q;
	ext->lpOffsetHz = spectrumClampLpOffsetHz(ext, lpOffsetHz);
	for (int i = 0; i < 10; i++)
	{
		if (harmonics && isfinite(harmonics[i]))
			ext->harmonicSeed[i] = (float)harmonics[i];
		else
			ext->harmonicSeed[i] = SPECTRUM_HARMONICS[i];
	}
	spectrumResetLocked(ext);
	jdsp_unlock(jdsp);
}

void SpectrumExtensionRefresh(JamesDSPLib *jdsp)
{
	SpectrumExtensionSetSamplingRate(jdsp, (int)jdsp->fs);
}

void SpectrumExtensionEnable(JamesDSPLib *jdsp)
{
	jdsp_lock(jdsp);
	SpectrumExtension *ext = &jdsp->spectrumExt;
	if (!ext->enabled)
	{
		ext->enabled = 1;
		spectrumResetLocked(ext);
	}
	jdsp_unlock(jdsp);
}

void SpectrumExtensionDisable(JamesDSPLib *jdsp)
{
	jdsp_lock(jdsp);
	jdsp->spectrumExt.enabled = 0;
	jdsp_unlock(jdsp);
}

void SpectrumExtensionProcess(JamesDSPLib *jdsp, size_t n)
{
	SpectrumExtension *ext = &jdsp->spectrumExt;
	if (!ext->enabled)
		return;

	for (size_t i = 0; i < n; i++)
	{
		double inL = jdsp->tmpBuffer[0][i];
		double inR = jdsp->tmpBuffer[1][i];

		double hpL = spectrumBiquadProcessSample(&ext->highpass[0], inL);
		double hpR = spectrumBiquadProcessSample(&ext->highpass[1], inR);

		double harmonicL = spectrumHarmonicProcessSample(&ext->harmonics[0], hpL);
		double harmonicR = spectrumHarmonicProcessSample(&ext->harmonics[1], hpR);

		double lpL = spectrumBiquadProcessSample(&ext->lowpass[0], harmonicL * ext->exciter);
		double lpR = spectrumBiquadProcessSample(&ext->lowpass[1], harmonicR * ext->exciter);
		double outL = (inL + lpL * ext->wetMix) * ext->postGain;
		double outR = (inR + lpR * ext->wetMix) * ext->postGain;
		if (ext->safetyEnabled)
		{
			outL = spectrumClampSafety(outL);
			outR = spectrumClampSafety(outR);
		}

		jdsp->tmpBuffer[0][i] = (float)outL;
		jdsp->tmpBuffer[1][i] = (float)outR;
	}
}
