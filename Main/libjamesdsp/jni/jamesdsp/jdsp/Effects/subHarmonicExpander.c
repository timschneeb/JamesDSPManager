/**
 * Sub-Harmonic Expander — synthesises sub-harmonics f/2, f/3, … f/10 of the
 * bass fundamental using a true frequency-divider approach.
 *
 * Algorithm:
 *   1. Low-pass filter (4th-order Linkwitz-Riley) to extract the bass band.
 *   2. Track the fundamental frequency via zero-crossing interval measurement.
 *   3. For each sub-harmonic k (f/(k+1)), generate a sine via phase accumulator
 *      at frequency f/(k+1), with amplitude scaled by the detected envelope
 *      and the user gain (0–100 %).
 *   4. Low-pass the synthetic output to keep only sub-bass content.
 *   5. Mix with original.
 *
 * This is a true sub-harmonic synthesiser (not a rectifier), so each
 * sub-harmonic is independently controllable.
 */
#include <math.h>
#include <string.h>
#include <float.h>
#include "../jdsp_header.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_2PI
#define M_2PI (2.0 * M_PI)
#endif

/* ---- 2nd-order Butterworth low-pass ---- */
static void shbiquadLP_init(SHBiquadStage *f)
{
	memset(f, 0, sizeof(*f));
}

static void shbiquadLP_set(SHBiquadStage *f, double fs, double fc)
{
	if (fc >= fs * 0.49) fc = fs * 0.49;
	if (fc < 1.0) fc = 1.0;
	double w0 = 2.0 * M_PI * fc / fs;
	double c = cos(w0);
	double s = sin(w0);
	double alpha = s / M_SQRT2; /* Q = 0.707 */
	double a0 = 1.0 + alpha;
	f->b0 = (float)((1.0 - c) / 2.0 / a0);
	f->b1 = (float)((1.0 - c) / a0);
	f->b2 = (float)((1.0 - c) / 2.0 / a0);
	f->a1 = (float)(-2.0 * c / a0);
	f->a2 = (float)((1.0 - alpha) / a0);
}

static inline float shbiquadLP_process(SHBiquadStage *f, float x)
{
	float y = f->b0 * x + f->b1 * f->x1 + f->b2 * f->x2 - f->a1 * f->y1 - f->a2 * f->y2;
	f->x2 = f->x1; f->x1 = x;
	f->y2 = f->y1; f->y1 = y;
	return y;
}

/* ---- 4th-order Linkwitz-Riley low-pass (cascade of two 2nd-order) ---- */
static void shbiquadLP4_init(SHBiquadLP4 *f)
{
	shbiquadLP_init(&f->stage1);
	shbiquadLP_init(&f->stage2);
}

static void shbiquadLP4_set(SHBiquadLP4 *f, double fs, double fc)
{
	shbiquadLP_set(&f->stage1, fs, fc);
	shbiquadLP_set(&f->stage2, fs, fc);
}

static inline float shbiquadLP4_process(SHBiquadLP4 *f, float x)
{
	return shbiquadLP_process(&f->stage2, shbiquadLP_process(&f->stage1, x));
}

/* ---- Envelope follower (peak detector with attack/release) ---- */
static inline float envelopeFollow(float env, float input, float attackCoef, float releaseCoef)
{
	float absIn = fabsf(input);
	if (absIn > env)
		return env + attackCoef * (absIn - env);
	else
		return env + releaseCoef * (absIn - env);
}

/* ---- Parameter setup ---- */
void SubHarmonicExpanderSetParam(JamesDSPLib *jdsp, float subHarmonicGains[9],
                                 float crossoverFreq, float mix)
{
	SubHarmonicExpander *sh = &jdsp->subHarmonicExpander;
	for (int i = 0; i < 9; i++)
		sh->subHarmonicGain[i] = subHarmonicGains[i];
	sh->mix = mix;
	sh->crossoverFreq = crossoverFreq;

	/* Band extraction low-pass */
	shbiquadLP4_set(&sh->bandLP[0], jdsp->fs, crossoverFreq);
	shbiquadLP4_set(&sh->bandLP[1], jdsp->fs, crossoverFreq);

	/* Output low-pass at 3/4 of crossover to keep sub-bass but allow some harmonics through */
	double outFreq = crossoverFreq * 0.75;
	if (outFreq < 20.0) outFreq = 20.0;
	shbiquadLP4_set(&sh->outLP[0], jdsp->fs, outFreq);
	shbiquadLP4_set(&sh->outLP[1], jdsp->fs, outFreq);

	/* Envelope follower coefficients */
	double attackMs = 2.0;   /* fast attack */
	double releaseMs = 50.0; /* moderate release */
	sh->envAttackCoef = (float)(1.0 - exp(-1.0 / (attackMs * 0.001 * jdsp->fs)));
	sh->envReleaseCoef = (float)(1.0 - exp(-1.0 / (releaseMs * 0.001 * jdsp->fs)));

	/* Reset state */
	sh->phase[0] = sh->phase[1] = 0.0f;
	sh->lastSign[0] = sh->lastSign[1] = 0;
	sh->samplesSinceLastCrossing[0] = sh->samplesSinceLastCrossing[1] = 0;
	sh->detectedFreq[0] = sh->detectedFreq[1] = 0.0f;
	sh->env[0] = sh->env[1] = 0.0f;
}

/* ---- Constructor / Enable / Disable ---- */
void SubHarmonicExpanderConstructor(JamesDSPLib *jdsp)
{
	SubHarmonicExpander *sh = &jdsp->subHarmonicExpander;
	memset(sh, 0, sizeof(*sh));
	for (int i = 0; i < 9; i++)
		sh->subHarmonicGain[i] = 0.0f;
	sh->mix = 1.0f;
	sh->crossoverFreq = 120.0f;
	sh->envAttackCoef = 0.001f;
	sh->envReleaseCoef = 0.0001f;
	shbiquadLP4_init(&sh->bandLP[0]);
	shbiquadLP4_init(&sh->bandLP[1]);
	shbiquadLP4_init(&sh->outLP[0]);
	shbiquadLP4_init(&sh->outLP[1]);
}

void SubHarmonicExpanderEnable(JamesDSPLib *jdsp)
{
	jdsp->subHarmonicExpanderEnabled = 1;
}

void SubHarmonicExpanderDisable(JamesDSPLib *jdsp)
{
	jdsp->subHarmonicExpanderEnabled = 0;
}

/* ---- Processing ----
 *
 * True sub-harmonic synthesis via frequency division:
 *
 * 1. Low-pass extract bass band → signal contains fundamental f.
 * 2. Detect fundamental frequency by counting samples between zero crossings.
 *    Each zero crossing = half period. freq = fs / (2 * crossingInterval).
 * 3. Track envelope of bass band.
 * 4. For each sub-harmonic k (k=1..10), synthesise sine at f/(k+1) using
 *    a phase accumulator: phase += 2*pi*freq/(k+1) / fs.
 *    Amplitude = envelope * gain_k.
 * 5. Sum all sub-harmonic sines → synthetic sub-bass.
 * 6. Low-pass output, mix with original.
 */
void SubHarmonicExpanderProcess(JamesDSPLib *jdsp, size_t n)
{
	SubHarmonicExpander *sh = &jdsp->subHarmonicExpander;
	float *xL = jdsp->tmpBuffer[0];
	float *xR = jdsp->tmpBuffer[1];
	float fs = jdsp->fs;

	for (size_t i = 0; i < n; i++)
	{
		/* === Process left channel === */
		float bandL = shbiquadLP4_process(&sh->bandLP[0], xL[i]);
		float bandR = shbiquadLP4_process(&sh->bandLP[1], xR[i]);

		/* Envelope follow the bass band */
		sh->env[0] = envelopeFollow(sh->env[0], bandL, sh->envAttackCoef, sh->envReleaseCoef);
		sh->env[1] = envelopeFollow(sh->env[1], bandR, sh->envAttackCoef, sh->envReleaseCoef);

		/* Zero-crossing detection for frequency tracking */
		for (int ch = 0; ch < 2; ch++)
		{
			float sample = (ch == 0) ? bandL : bandR;
			int sign = (sample >= 0.0f) ? 1 : -1;
			sh->samplesSinceLastCrossing[ch]++;

			if (sign != sh->lastSign[ch] && sh->lastSign[ch] != 0)
			{
				/* Zero crossing detected — update frequency estimate */
				int interval = sh->samplesSinceLastCrossing[ch];
				if (interval > 2 && interval < (int)(fs / 20.0f))
				{
					/* Smooth frequency estimate */
					float newFreq = fs / (2.0f * (float)interval);
					if (newFreq > 20.0f && newFreq < sh->crossoverFreq * 2.0f)
					{
						/* One-pole smoothing */
						sh->detectedFreq[ch] = 0.3f * newFreq + 0.7f * sh->detectedFreq[ch];
					}
				}
				sh->samplesSinceLastCrossing[ch] = 0;
			}
			sh->lastSign[ch] = sign;
		}

		/* Synthesise sub-harmonics */
		float synthL = 0.0f, synthR = 0.0f;
		float avgFreq = 0.5f * (sh->detectedFreq[0] + sh->detectedFreq[1]);

		if (avgFreq > 20.0f)
		{
			for (int k = 0; k < 9; k++)
			{
				if (sh->subHarmonicGain[k] > 0.0f)
				{
					/* Sub-harmonic k+1: frequency = f / (k+2)
					 * k=0 → f/2, k=1 → f/3, ..., k=8 → f/10 */
					float subFreq = avgFreq / (float)(k + 2);
					if (subFreq < 10.0f) continue; /* skip sub-sonic */

					/* Phase accumulator */
					float phaseInc = M_2PI * subFreq / fs;
					sh->phase[0] += phaseInc;
					sh->phase[1] += phaseInc;
					if (sh->phase[0] > M_PI) sh->phase[0] -= M_2PI;
					if (sh->phase[1] > M_PI) sh->phase[1] -= M_2PI;

					/* Gain: user gain (0-100%) * envelope * normalisation
					 * Higher sub-harmonics are quieter, so boost by (k+2)^0.5 */
					float gainScale = sh->subHarmonicGain[k] * 0.01f;
					float ampL = sh->env[0] * gainScale;
					float ampR = sh->env[1] * gainScale;

					synthL += ampL * sinf(sh->phase[0]);
					synthR += ampR * sinf(sh->phase[1]);
				}
			}
		}

		/* Low-pass the synthetic output */
		synthL = shbiquadLP4_process(&sh->outLP[0], synthL);
		synthR = shbiquadLP4_process(&sh->outLP[1], synthR);

		/* Mix with original */
		xL[i] = xL[i] + synthL * sh->mix;
		xR[i] = xR[i] + synthR * sh->mix;
	}
}
