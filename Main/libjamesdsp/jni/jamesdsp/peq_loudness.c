/*
 * peq_loudness.c — Pure-C PEQ biquad cascade and loudness correction.
 *
 * PEQ: RBJ Audio EQ Cookbook biquad coefficients, cascaded per-band
 *      with per-band L/R channel routing. Ported from BiQuad.h /
 *      ParametricEqProcessor.cpp.
 *
 * Loudness: Fletcher-Munson compensation with low-shelf (75 Hz) and
 *           high-shelf (10 kHz). Ported from LoudnessCorrectionProcessor.cpp
 *           (which was ported from EqualizerAPO).
 *
 * See peq_loudness.h for the public API and parameter ID documentation.
 */

#include "peq_loudness.h"
#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_LN2
#define M_LN2 0.69314718055994530942
#endif

/* ================================================================ */
/* Biquad                                                           */
/* ================================================================ */

void biquad_init(biquad_t *bq, peq_filter_type_t type, double dbGain,
                 double freq, double srate, double bandwidthOrQOrS,
                 bool isBandwidthOrS)
{
    memset(bq, 0, sizeof(*bq));

    /* Clamp frequency */
    double nyquist = srate * 0.5;
    if (freq < 1.0) freq = 1.0;
    if (freq > nyquist * 0.98) freq = nyquist * 0.98;

    double A;
    if (type == PEQ_TYPE_PEAKING || type == PEQ_TYPE_LOW_SHELF || type == PEQ_TYPE_HIGH_SHELF)
        A = pow(10.0, dbGain / 40.0);
    else
        A = pow(10.0, dbGain / 20.0);

    double omega = 2.0 * M_PI * freq / srate;
    double sn = sin(omega);
    double cs = cos(omega);
    double alpha;

    if (!isBandwidthOrS) {
        /* Q mode */
        double q = bandwidthOrQOrS;
        if (q < 0.1) q = 0.1;
        if (q > 24.0) q = 24.0;
        alpha = sn / (2.0 * q);
    } else if (type == PEQ_TYPE_LOW_SHELF || type == PEQ_TYPE_HIGH_SHELF) {
        /* S mode */
        alpha = sn / 2.0 * sqrt((A + 1.0 / A) * (1.0 / bandwidthOrQOrS - 1.0) + 2.0);
    } else {
        /* BW mode */
        alpha = sn * sinh(M_LN2 / 2.0 * bandwidthOrQOrS * omega / sn);
    }

    double beta = 2.0 * sqrt(A) * alpha;

    double b0, b1, b2, a0, a1, a2;

    switch (type) {
    case PEQ_TYPE_LOW_PASS:
        b0 = (1.0 - cs) / 2.0;
        b1 = 1.0 - cs;
        b2 = (1.0 - cs) / 2.0;
        a0 = 1.0 + alpha;
        a1 = -2.0 * cs;
        a2 = 1.0 - alpha;
        break;
    case PEQ_TYPE_HIGH_PASS:
        b0 = (1.0 + cs) / 2.0;
        b1 = -(1.0 + cs);
        b2 = (1.0 + cs) / 2.0;
        a0 = 1.0 + alpha;
        a1 = -2.0 * cs;
        a2 = 1.0 - alpha;
        break;
    case PEQ_TYPE_BAND_PASS:
        b0 = alpha;
        b1 = 0.0;
        b2 = -alpha;
        a0 = 1.0 + alpha;
        a1 = -2.0 * cs;
        a2 = 1.0 - alpha;
        break;
    case PEQ_TYPE_NOTCH:
        b0 = 1.0;
        b1 = -2.0 * cs;
        b2 = 1.0;
        a0 = 1.0 + alpha;
        a1 = -2.0 * cs;
        a2 = 1.0 - alpha;
        break;
    case PEQ_TYPE_ALL_PASS:
        b0 = 1.0 - alpha;
        b1 = -2.0 * cs;
        b2 = 1.0 + alpha;
        a0 = 1.0 + alpha;
        a1 = -2.0 * cs;
        a2 = 1.0 - alpha;
        break;
    case PEQ_TYPE_PEAKING:
        b0 = 1.0 + (alpha * A);
        b1 = -2.0 * cs;
        b2 = 1.0 - (alpha * A);
        a0 = 1.0 + (alpha / A);
        a1 = -2.0 * cs;
        a2 = 1.0 - (alpha / A);
        break;
    case PEQ_TYPE_LOW_SHELF:
        b0 = A * ((A + 1.0) - (A - 1.0) * cs + beta);
        b1 = 2.0 * A * ((A - 1.0) - (A + 1.0) * cs);
        b2 = A * ((A + 1.0) - (A - 1.0) * cs - beta);
        a0 = (A + 1.0) + (A - 1.0) * cs + beta;
        a1 = -2.0 * ((A - 1.0) + (A + 1.0) * cs);
        a2 = (A + 1.0) + (A - 1.0) * cs - beta;
        break;
    case PEQ_TYPE_HIGH_SHELF:
        b0 = A * ((A + 1.0) + (A - 1.0) * cs + beta);
        b1 = -2.0 * A * ((A - 1.0) + (A + 1.0) * cs);
        b2 = A * ((A + 1.0) + (A - 1.0) * cs - beta);
        a0 = (A + 1.0) - (A - 1.0) * cs + beta;
        a1 = 2.0 * ((A - 1.0) - (A + 1.0) * cs);
        a2 = (A + 1.0) - (A - 1.0) * cs - beta;
        break;
    default:
        b0 = 1.0; b1 = 0.0; b2 = 0.0;
        a0 = 1.0; a1 = 0.0; a2 = 0.0;
        break;
    }

    /* Normalize by a0 */
    bq->b0 = b0 / a0;
    bq->b1 = b1 / a0;
    bq->b2 = b2 / a0;
    bq->a1 = a1 / a0;
    bq->a2 = a2 / a0;
}

/* ================================================================ */
/* PEQ cascade                                                      */
/* ================================================================ */

void peq_cascade_init(peq_cascade_t *peq)
{
    memset(peq, 0, sizeof(*peq));
    peq->sampleRate = 48000.0;
    peq->preampLinear = 1.0;
}

void peq_cascade_configure(peq_cascade_t *peq, double sampleRate,
                           const float *data, int numFloats)
{
    peq->sampleRate = sampleRate;
    peq->bandCount = 0;

    if (numFloats < 2)
        return;

    double preampDb = (double)data[1];
    peq->preampLinear = pow(10.0, preampDb / 20.0);

    /* Each band uses 5 floats: freq, gain, q, filterType, channelMode */
    int maxBands = (numFloats - 2) / 5;
    if (maxBands > PEQ_MAX_BANDS)
        maxBands = PEQ_MAX_BANDS;

    for (int i = 0; i < maxBands; i++) {
        int base = 2 + i * 5;
        double freq = (double)data[base + 0];
        double gain = (double)data[base + 1];
        double q    = (double)data[base + 2];
        int ftype   = (int)roundf(data[base + 3]);
        int cmode   = (int)roundf(data[base + 4]);

        /* freq <= 0 means band is disabled/unused */
        if (freq <= 0.0)
            continue;

        /* Clamp Q */
        if (q < 0.1) q = 0.1;
        if (q > 24.0) q = 24.0;

        /* Clamp filter type */
        if (ftype < 0 || ftype > 8)
            ftype = 0; /* default to peaking */

        /* Clamp channel mode */
        if (cmode < 0 || cmode > 2)
            cmode = 0;

        peq_band_t *band = &peq->bands[peq->bandCount];
        band->channelMode = cmode;
        band->enabled = true;
        band->isPreamp = false;
        band->preampLinear = 1.0;

        /* Preamp band: flat gain, no biquad */
        if (ftype == PEQ_TYPE_PREAMP) {
            band->isPreamp = true;
            band->preampLinear = pow(10.0, gain / 20.0);
            peq->bandCount++;
            continue;
        }

        /* Q mode (isBandwidthOrS = false) */
        biquad_init(&band->bqL, (peq_filter_type_t)ftype, gain, freq, sampleRate, q, false);
        biquad_init(&band->bqR, (peq_filter_type_t)ftype, gain, freq, sampleRate, q, false);

        peq->bandCount++;
    }
}

void peq_cascade_set_enabled(peq_cascade_t *peq, bool enabled)
{
    peq->enabled = enabled;
}

void peq_cascade_process_interleaved(peq_cascade_t *peq, float *data, size_t numFrames)
{
    if (!peq->enabled || peq->bandCount == 0)
        return;

    /* Apply preamp */
    if (peq->preampLinear != 1.0) {
        for (size_t i = 0; i < numFrames * 2; i++)
            data[i] = (float)(data[i] * peq->preampLinear);
    }

    /* Process through biquad cascade */
    for (size_t i = 0; i < numFrames; i++) {
        double sampleL = data[i * 2];
        double sampleR = data[i * 2 + 1];

        for (int b = 0; b < peq->bandCount; b++) {
            peq_band_t *band = &peq->bands[b];
            if (!band->enabled)
                continue;

            /* Preamp band: flat gain, no biquad processing */
            if (band->isPreamp) {
                const double g = band->preampLinear;
                switch (band->channelMode) {
                case PEQ_CHAN_BOTH:
                    sampleL *= g;
                    sampleR *= g;
                    break;
                case PEQ_CHAN_LEFT:
                    sampleL *= g;
                    break;
                case PEQ_CHAN_RIGHT:
                    sampleR *= g;
                    break;
                }
                continue;
            }

            biquad_remove_denormals(&band->bqL);
            biquad_remove_denormals(&band->bqR);

            switch (band->channelMode) {
            case PEQ_CHAN_BOTH:
                sampleL = biquad_process(&band->bqL, sampleL);
                sampleR = biquad_process(&band->bqR, sampleR);
                break;
            case PEQ_CHAN_LEFT:
                sampleL = biquad_process(&band->bqL, sampleL);
                break;
            case PEQ_CHAN_RIGHT:
                sampleR = biquad_process(&band->bqR, sampleR);
                break;
            }
        }

        data[i * 2]     = (float)sampleL;
        data[i * 2 + 1] = (float)sampleR;
    }
}

/* ================================================================ */
/* Loudness correction                                              */
/* ================================================================ */

void loudness_init(loudness_t *l)
{
    memset(l, 0, sizeof(*l));
    l->sampleRate = 48000.0;
    l->attenuation = 1.0;
    l->attFactor = 1.0;
    l->neutral = true;
    l->mode = 0;
    /* Default tuning parameters (match Kotlin defaults) */
    l->lsFreq = 75.0;
    l->lsSlope = 0.52;
    l->lsRatio = 0.55;
    l->hsFreq = 10000.0;
    l->hsSlope = 0.9;
    l->hsRatio = 0.225;
    l->isoBasePhon = 80.0;
    l->isoQ = 4.32;
    l->lastComputedVolume = 1e30; /* impossible → forces recompute */
    l->coeffsDirty = true;
    /* Default subsonic filter parameters */
    l->subsonicEnabled = false;
    l->subsonicFreq = 20.0;
    l->subsonicOrder = 4;
    l->subsonicQ = 0.707;
    /* ISO 226:2023 mode state */
    l->isoPreampLinear = 1.0;
    l->isoActiveBandCount = 0;
    for (int i = 0; i < ISO226_NUM_BANDS; i++) {
        l->isoGains[i] = 0.0;
        /* Set passthrough (b0=1) so unconfigured bands don't silence audio */
        l->isoBandsL[i].b0 = 1.0; l->isoBandsL[i].b1 = 0; l->isoBandsL[i].b2 = 0;
        l->isoBandsL[i].a1 = 0; l->isoBandsL[i].a2 = 0;
        l->isoBandsR[i].b0 = 1.0; l->isoBandsR[i].b1 = 0; l->isoBandsR[i].b2 = 0;
        l->isoBandsR[i].a1 = 0; l->isoBandsR[i].a2 = 0;
    }
}

static void loudness_get_low_shelf_params(loudness_t *l, double volume,
                                          double *freq, double *s,
                                          double *gain, double *preAmp)
{
    *freq = l->lsFreq;
    *s = l->lsSlope;
    double volDiff = l->referenceLevel - l->referenceOffset - volume;
    if (volDiff > 0.0) {
        *gain = volDiff * l->lsRatio / (1.0 - l->lsRatio) * l->attenuation;
        *preAmp = -*gain;
    } else if (volDiff < 0.0) {
        *preAmp = 0.0;
        *gain = volDiff * l->lsRatio * exp(volDiff / 90.0) * l->attenuation;
    } else {
        *gain = 0.0;
        *preAmp = 0.0;
    }
}

static void loudness_get_high_shelf_params(loudness_t *l, double volume,
                                           double *freq, double *s, double *gain)
{
    *freq = l->hsFreq;
    *s = l->hsSlope;
    double volDiff = l->referenceLevel - l->referenceOffset - volume;
    if (volDiff > 0.0) {
        *gain = volDiff * l->hsRatio * exp(-volDiff / 100.0) * l->attenuation;
    } else if (volDiff < 0.0) {
        *gain = volDiff * l->hsRatio * 0.778 * exp(volDiff / 80.0) * l->attenuation;
    } else {
        *gain = 0.0;
    }
}

/* ================================================================ */
/* ISO 226:2023 equal-loudness-level contour data and computation   */
/* ================================================================ */

/* ISO 226:2023 Table 1: 29 preferred one-third-octave frequencies (Hz) */
static const double iso226Freqs[ISO226_NUM_BANDS] = {
    20.0, 25.0, 31.5, 40.0, 50.0, 63.0, 80.0, 100.0, 125.0, 160.0,
    200.0, 250.0, 315.0, 400.0, 500.0, 630.0, 800.0, 1000.0, 1250.0, 1600.0,
    2000.0, 2500.0, 3150.0, 4000.0, 5000.0, 6300.0, 8000.0, 10000.0, 12500.0
};

/* ISO 226:2023 Table 1: exponent α_f */
static const double iso226Af[ISO226_NUM_BANDS] = {
    0.635, 0.602, 0.569, 0.537, 0.509, 0.482, 0.456, 0.433, 0.412, 0.391,
    0.373, 0.357, 0.343, 0.330, 0.320, 0.311, 0.303, 0.300, 0.295, 0.292,
    0.290, 0.290, 0.289, 0.289, 0.289, 0.293, 0.303, 0.323, 0.354
};

/* ISO 226:2023 Table 1: magnitude of linear transfer function L_U (dB) */
static const double iso226Lu[ISO226_NUM_BANDS] = {
    -31.5, -27.2, -23.1, -19.3, -16.1, -13.1, -10.4, -8.2, -6.3, -4.6,
    -3.2, -2.1, -1.2, -0.5, 0.0, 0.4, 0.5, 0.0, -2.7, -4.2,
    -1.2, 1.4, 2.3, 1.0, -2.3, -7.2, -11.2, -10.9, -3.5
};

/* ISO 226:2023 Table 1: threshold of hearing T_f (dB) */
static const double iso226Tf[ISO226_NUM_BANDS] = {
    78.1, 68.7, 59.5, 51.1, 44.0, 37.5, 31.5, 26.5, 22.1, 17.9,
    14.4, 11.4, 8.6, 6.2, 4.4, 3.0, 2.2, 2.4, 3.5, 1.7,
    -1.3, -4.2, -6.0, -5.4, -1.5, 6.0, 12.6, 13.9, 12.3
};

/* ISO 226:2023 Formula (1): SPL (dB) from loudness level L_N (phons) */
static double iso226_spl_from_phon(double af, double lu, double tf, double lN)
{
    double term1 = pow(4e-10, 0.3 - af) *
                   (pow(10.0, 0.03 * lN) - pow(10.0, 0.072));
    double term2 = pow(10.0, af * (tf + lu) / 10.0);
    double inner = term1 + term2;
    if (inner <= 0.0)
        return 0.0;
    return (10.0 / af) * log10(inner) - lu;
}

/* Compute ISO 226:2023 loudness correction gains for all 29 bands.
 * Returns the maximum positive gain (for pre-amp calculation). */
static double loudness_compute_iso226_gains(loudness_t *l, double volume,
                                            double gains[ISO226_NUM_BANDS])
{
    double basePhon = l->isoBasePhon;
    double volDiff = l->referenceLevel - l->referenceOffset - volume;
    double refPhon = basePhon;
    double curPhon = basePhon - volDiff;

    /* Clamp phon values to valid ISO 226 range [20, 90] */
    if (refPhon < 20.0) refPhon = 20.0;
    if (refPhon > 90.0) refPhon = 90.0;
    if (curPhon < 20.0) curPhon = 20.0;
    if (curPhon > 90.0) curPhon = 90.0;

    double maxGain = 0.0;
    int i;
    for (i = 0; i < ISO226_NUM_BANDS; i++) {
        double refSpl = iso226_spl_from_phon(iso226Af[i], iso226Lu[i],
                                             iso226Tf[i], refPhon);
        double curSpl = iso226_spl_from_phon(iso226Af[i], iso226Lu[i],
                                             iso226Tf[i], curPhon);
        double gain = curSpl - refSpl + (refPhon - curPhon);
        gains[i] = gain;
        if (gain > maxGain)
            maxGain = gain;
    }
    return maxGain;
}

/* Compute the actual peak gain of the ISO 226 biquad cascade.
 * Evaluates the combined response of all 29 peaking biquads at fine
 * frequency resolution and returns the maximum gain in dB.
 * This accounts for filter overlap between adjacent bands. */
static double loudness_compute_cascade_peak(loudness_t *l)
{
    double peak = 0.0;
    int i;

    /* Evaluate at each ISO band center and halfway between adjacent centers */
    for (i = 0; i < ISO226_NUM_BANDS; i++)
    {
        /* Test at the band center */
        double f = iso226Freqs[i];
        double total = 0.0;
        int j;
        for (j = 0; j < ISO226_NUM_BANDS; j++)
        {
            double gain = l->isoGains[j] * l->attenuation;
            if (fabs(gain) < 0.05)
                continue;
            /* Evaluate peaking biquad magnitude at frequency f */
            double w  = 2.0 * M_PI * f / l->sampleRate;
            double w0 = 2.0 * M_PI * iso226Freqs[j] / l->sampleRate;
            if (iso226Freqs[j] > l->sampleRate * 0.5 * 0.98)
                continue; /* above nyquist, no contribution */
            double A = pow(10.0, gain / 40.0);
            double alpha = sin(w0) / (2.0 * l->isoQ);
            double b0 = 1.0 + alpha * A;
            double b1 = -2.0 * cos(w0);
            double b2 = 1.0 - alpha * A;
            double a0 = 1.0 + alpha / A;
            double a1 = -2.0 * cos(w0);
            double a2 = 1.0 - alpha / A;
            double cosw  = cos(w);
            double cos2w = cos(2.0 * w);
            double sinw  = sin(w);
            double sin2w = sin(2.0 * w);
            double num_re = b0 + b1 * cosw + b2 * cos2w;
            double num_im = b1 * sinw + b2 * sin2w;
            double den_re = a0 + a1 * cosw + a2 * cos2w;
            double den_im = a1 * sinw + a2 * sin2w;
            double num_sq = num_re * num_re + num_im * num_im;
            double den_sq = den_re * den_re + den_im * den_im;
            if (den_sq > 1e-20 && num_sq > 0.0)
            {
                double mag_sq = num_sq / den_sq;
                total += 10.0 * log10(mag_sq);
            }
        }
        if (total > peak)
            peak = total;

        /* Also test at geometric mean between this and next band */
        if (i < ISO226_NUM_BANDS - 1)
        {
            double fMid = sqrt(iso226Freqs[i] * iso226Freqs[i + 1]);
            total = 0.0;
            for (j = 0; j < ISO226_NUM_BANDS; j++)
            {
                double gain = l->isoGains[j] * l->attenuation;
                if (fabs(gain) < 0.05)
                    continue;
                double w  = 2.0 * M_PI * fMid / l->sampleRate;
                double w0 = 2.0 * M_PI * iso226Freqs[j] / l->sampleRate;
                if (iso226Freqs[j] > l->sampleRate * 0.5 * 0.98)
                    continue;
                double A = pow(10.0, gain / 40.0);
                double alpha = sin(w0) / (2.0 * l->isoQ);
                double b0 = 1.0 + alpha * A;
                double b1 = -2.0 * cos(w0);
                double b2 = 1.0 - alpha * A;
                double a0 = 1.0 + alpha / A;
                double a1 = -2.0 * cos(w0);
                double a2 = 1.0 - alpha / A;
                double cosw  = cos(w);
                double cos2w = cos(2.0 * w);
                double sinw  = sin(w);
                double sin2w = sin(2.0 * w);
                double num_re = b0 + b1 * cosw + b2 * cos2w;
                double num_im = b1 * sinw + b2 * sin2w;
                double den_re = a0 + a1 * cosw + a2 * cos2w;
                double den_im = a1 * sinw + a2 * sin2w;
                double num_sq = num_re * num_re + num_im * num_im;
                double den_sq = den_re * den_re + den_im * den_im;
                if (den_sq > 1e-20 && num_sq > 0.0)
                {
                    double mag_sq = num_sq / den_sq;
                    total += 10.0 * log10(mag_sq);
                }
            }
            if (total > peak)
                peak = total;
        }
    }

    /* Ensure peak is at least maxGain (sanity floor) */
    double maxIndiv = 0.0;
    for (i = 0; i < ISO226_NUM_BANDS; i++)
    {
        double gain = l->isoGains[i] * l->attenuation;
        if (gain > maxIndiv)
            maxIndiv = gain;
    }
    if (peak < maxIndiv)
        peak = maxIndiv;

    return peak;
}

static void loudness_recompute(loudness_t *l, double volume)
{
    if (l->mode == 1) {
        /* ---- ISO 226:2023 mode ---- */
        double maxGain = loudness_compute_iso226_gains(l, volume, l->isoGains);

        /* Check if correction is near-zero (bypass for best quality) */
        double maxAbs = 0.0;
        int i;
        for (i = 0; i < ISO226_NUM_BANDS; i++) {
            double absGain = fabs(l->isoGains[i] * l->attenuation);
            if (absGain > maxAbs)
                maxAbs = absGain;
        }
        l->neutral = (maxAbs < 0.2);

        /* Build peaking biquads for each band */
        double qThirdOctave = l->isoQ;
        double nyquist = l->sampleRate * 0.5;

        l->isoActiveBandCount = 0;
        for (i = 0; i < ISO226_NUM_BANDS; i++) {
            double gain = l->isoGains[i] * l->attenuation;
            if (fabs(gain) < 0.05) {
                /* Skip near-zero bands — set passthrough (b0=1) */
                biquad_reset(&l->isoBandsL[i]);
                biquad_reset(&l->isoBandsR[i]);
                l->isoBandsL[i].b0 = 1.0; l->isoBandsL[i].b1 = 0; l->isoBandsL[i].b2 = 0;
                l->isoBandsL[i].a1 = 0; l->isoBandsL[i].a2 = 0;
                l->isoBandsR[i].b0 = 1.0; l->isoBandsR[i].b1 = 0; l->isoBandsR[i].b2 = 0;
                l->isoBandsR[i].a1 = 0; l->isoBandsR[i].a2 = 0;
                continue;
            }

            double freq = iso226Freqs[i];
            if (freq < 1.0) freq = 1.0;
            if (freq > nyquist * 0.98) freq = nyquist * 0.98;

            biquad_init(&l->isoBandsL[i], PEQ_TYPE_PEAKING, gain, freq,
                        l->sampleRate, qThirdOctave, false);
            biquad_init(&l->isoBandsR[i], PEQ_TYPE_PEAKING, gain, freq,
                        l->sampleRate, qThirdOctave, false);
            l->isoActiveBandCount++;
        }

        /* Pre-amp: offset by the ACTUAL peak of the combined cascade.
         * Adjacent peaking filters overlap and sum, so the cascade peak
         * is higher than any individual band gain. Evaluate the combined
         * response at fine frequency resolution to find the true peak. */
        double cascadePeak = 0.0;
        if (l->isoActiveBandCount > 0) {
            cascadePeak = loudness_compute_cascade_peak(l);
        }
        double preAmpDb = -cascadePeak;
        l->attFactor = exp(preAmpDb / 6.0 * log(2.0));

        l->isoPreampLinear = l->attFactor;
        /* Fall through to subsonic filter rebuild below */
    }
    else
    {
    /* ---- Classic mode (Fletcher-Munson two-shelf) ---- */
    double freqLS, sLS, gainLS, preAmp;
    double freqHS, sHS, gainHS;

    loudness_get_low_shelf_params(l, volume, &freqLS, &sLS, &gainLS, &preAmp);
    l->attFactor = exp(preAmp / 6.0 * log(2.0));

    loudness_get_high_shelf_params(l, volume + preAmp, &freqHS, &sHS, &gainHS);

    double absGainLS = fabs(gainLS);
    double absGainHS = fabs(gainHS);
    l->neutral = (absGainLS < 0.2 && absGainHS < 0.2);

    /* S-slope shelf mode (isBandwidthOrS = true) */
    biquad_init(&l->lowShelfL,  PEQ_TYPE_LOW_SHELF,  gainLS, freqLS, l->sampleRate, sLS, true);
    biquad_init(&l->lowShelfR,  PEQ_TYPE_LOW_SHELF,  gainLS, freqLS, l->sampleRate, sLS, true);
    biquad_init(&l->highShelfL, PEQ_TYPE_HIGH_SHELF, gainHS, freqHS, l->sampleRate, sHS, true);
    biquad_init(&l->highShelfR, PEQ_TYPE_HIGH_SHELF, gainHS, freqHS, l->sampleRate, sHS, true);
    }

    /* ---- Subsonic (infrasonic) high-pass filter ---- */
    /* Rebuild regardless of mode; applied before loudness correction. */
    if (l->subsonicEnabled) {
        double nyquist = l->sampleRate * 0.5;
        double f = l->subsonicFreq;
        if (f < 1.0) f = 1.0;
        if (f > nyquist * 0.98) f = nyquist * 0.98;
        int numStages = (l->subsonicOrder + 1) / 2;
        int i;
        for (i = 0; i < numStages && i < MAX_SUBSONIC_STAGES; i++) {
            biquad_init(&l->subsonicL[i], PEQ_TYPE_HIGH_PASS, 0.0, f,
                        l->sampleRate, l->subsonicQ, false);
            biquad_init(&l->subsonicR[i], PEQ_TYPE_HIGH_PASS, 0.0, f,
                        l->sampleRate, l->subsonicQ, false);
        }
        for (; i < MAX_SUBSONIC_STAGES; i++) {
            biquad_reset(&l->subsonicL[i]);
            biquad_reset(&l->subsonicR[i]);
            /* Set flat coefficients */
            l->subsonicL[i].b0 = 1.0; l->subsonicL[i].b1 = 0; l->subsonicL[i].b2 = 0;
            l->subsonicL[i].a1 = 0; l->subsonicL[i].a2 = 0;
            l->subsonicR[i].b0 = 1.0; l->subsonicR[i].b1 = 0; l->subsonicR[i].b2 = 0;
            l->subsonicR[i].a1 = 0; l->subsonicR[i].a2 = 0;
        }
    }
}

void loudness_configure(loudness_t *l, double sampleRate,
                        double refLevel, double refOffset,
                        double attenuation, double volumeDb, int mode)
{
    l->sampleRate = sampleRate;
    l->referenceLevel = refLevel;
    l->referenceOffset = refOffset;
    l->attenuation = (attenuation < 0.0) ? 0.0 : (attenuation > 2.0 ? 2.0 : attenuation);
    l->volumeDb = volumeDb;
    l->mode = mode;
    l->lastComputedVolume = 1e30; /* force recompute */
    l->coeffsDirty = true;
    l->attFactor = 1.0;
    l->neutral = true;
    /* Reset ISO 226 state */
    l->isoPreampLinear = 1.0;
    l->isoActiveBandCount = 0;
    for (int i = 0; i < ISO226_NUM_BANDS; i++) {
        l->isoGains[i] = 0.0;
        biquad_reset(&l->isoBandsL[i]);
        biquad_reset(&l->isoBandsR[i]);
        /* Set passthrough (b0=1) so unconfigured bands don't silence audio */
        l->isoBandsL[i].b0 = 1.0; l->isoBandsL[i].b1 = 0; l->isoBandsL[i].b2 = 0;
        l->isoBandsL[i].a1 = 0; l->isoBandsL[i].a2 = 0;
        l->isoBandsR[i].b0 = 1.0; l->isoBandsR[i].b1 = 0; l->isoBandsR[i].b2 = 0;
        l->isoBandsR[i].a1 = 0; l->isoBandsR[i].a2 = 0;
    }
}

void loudness_set_tuning(loudness_t *l,
                         double lsFreq, double lsSlope, double lsRatio,
                         double hsFreq, double hsSlope, double hsRatio,
                         double isoBasePhon, double isoQ)
{
    l->lsFreq = lsFreq;
    l->lsSlope = lsSlope;
    l->lsRatio = lsRatio;
    l->hsFreq = hsFreq;
    l->hsSlope = hsSlope;
    l->hsRatio = hsRatio;
    l->isoBasePhon = isoBasePhon;
    l->isoQ = isoQ;
    l->coeffsDirty = true;
}

void loudness_set_subsonic(loudness_t *l, bool enable, double freq,
                           int order, double qFactor)
{
    l->subsonicEnabled = enable;
    l->subsonicFreq = (freq < 5.0) ? 5.0 : (freq > 200.0 ? 200.0 : freq);
    l->subsonicOrder = (order < 1) ? 1 : (order > 8 ? 8 : order);
    l->subsonicQ = (qFactor < 0.1) ? 0.1 : (qFactor > 2.0 ? 2.0 : qFactor);
    l->coeffsDirty = true;
}

void loudness_set_volume(loudness_t *l, double volumeDb)
{
    l->volumeDb = volumeDb;
    l->coeffsDirty = true;
}

void loudness_set_enabled(loudness_t *l, bool enabled)
{
    l->enabled = enabled;
}

void loudness_process_interleaved(loudness_t *l, float *data, size_t numFrames)
{
    if (!l->enabled)
        return;

    /* Lazily recompute coefficients if volume or parameters changed */
    if (l->coeffsDirty) {
        loudness_recompute(l, l->volumeDb);
        l->lastComputedVolume = l->volumeDb;
        l->coeffsDirty = false;
    }

    /* Determine subsonic stages */
    int numSubsonicStages = 0;
    if (l->subsonicEnabled)
        numSubsonicStages = (l->subsonicOrder + 1) / 2;

    if (l->neutral && numSubsonicStages == 0)
        return;

    /* If only subsonic is active (neutral loudness) */
    if (l->neutral) {
        for (size_t i = 0; i < numFrames; i++) {
            double sL = data[i * 2];
            double sR = data[i * 2 + 1];
            for (int s = 0; s < numSubsonicStages && s < MAX_SUBSONIC_STAGES; s++) {
                sL = biquad_process(&l->subsonicL[s], sL);
                sR = biquad_process(&l->subsonicR[s], sR);
            }
            data[i * 2]     = (float)sL;
            data[i * 2 + 1] = (float)sR;
        }
        for (int s = 0; s < numSubsonicStages && s < MAX_SUBSONIC_STAGES; s++) {
            biquad_remove_denormals(&l->subsonicL[s]);
            biquad_remove_denormals(&l->subsonicR[s]);
        }
        return;
    }

    if (l->mode == 1) {
        /* ISO 226:2023 mode — 29-band peaking EQ cascade with pre-amp */
        if (l->isoActiveBandCount == 0)
            return;

        double preamp = l->isoPreampLinear;
        for (size_t i = 0; i < numFrames; i++) {
            double sampleL = data[i * 2];
            double sampleR = data[i * 2 + 1];

            /* Subsonic HP filter before loudness correction */
            for (int s = 0; s < numSubsonicStages && s < MAX_SUBSONIC_STAGES; s++) {
                sampleL = biquad_process(&l->subsonicL[s], sampleL);
                sampleR = biquad_process(&l->subsonicR[s], sampleR);
            }

            sampleL *= preamp;
            sampleR *= preamp;

            for (int b = 0; b < ISO226_NUM_BANDS; b++) {
                sampleL = biquad_process(&l->isoBandsL[b], sampleL);
                sampleR = biquad_process(&l->isoBandsR[b], sampleR);
            }

            data[i * 2]     = (float)sampleL;
            data[i * 2 + 1] = (float)sampleR;
        }

        /* Remove denormals */
        for (int b = 0; b < ISO226_NUM_BANDS; b++) {
            biquad_remove_denormals(&l->isoBandsL[b]);
            biquad_remove_denormals(&l->isoBandsR[b]);
        }
        for (int s = 0; s < numSubsonicStages && s < MAX_SUBSONIC_STAGES; s++) {
            biquad_remove_denormals(&l->subsonicL[s]);
            biquad_remove_denormals(&l->subsonicR[s]);
        }
        return;
    }

    /* Classic mode (Fletcher-Munson two-shelf) */
    for (size_t i = 0; i < numFrames; i++) {
        double sampleL = data[i * 2];
        double sampleR = data[i * 2 + 1];

        /* Subsonic HP filter before loudness correction */
        for (int s = 0; s < numSubsonicStages && s < MAX_SUBSONIC_STAGES; s++) {
            sampleL = biquad_process(&l->subsonicL[s], sampleL);
            sampleR = biquad_process(&l->subsonicR[s], sampleR);
        }

        /* Low shelf → pre-amp attenuation → high shelf */
        sampleL = biquad_process(&l->lowShelfL, sampleL);
        sampleR = biquad_process(&l->lowShelfR, sampleR);

        sampleL *= l->attFactor;
        sampleR *= l->attFactor;

        sampleL = biquad_process(&l->highShelfL, sampleL);
        sampleR = biquad_process(&l->highShelfR, sampleR);

        data[i * 2]     = (float)sampleL;
        data[i * 2 + 1] = (float)sampleR;
    }

    /* Remove denormals */
    biquad_remove_denormals(&l->lowShelfL);
    biquad_remove_denormals(&l->lowShelfR);
    biquad_remove_denormals(&l->highShelfL);
    biquad_remove_denormals(&l->highShelfR);
    for (int s = 0; s < numSubsonicStages && s < MAX_SUBSONIC_STAGES; s++) {
        biquad_remove_denormals(&l->subsonicL[s]);
        biquad_remove_denormals(&l->subsonicR[s]);
    }
}
