/**
 * Harmonic Expander — Chebyshev polynomial waveshaper with per-harmonic gain.
 *
 * Generates harmonics 1–10 of the fundamental using Chebyshev polynomials
 * T_k(x). Each harmonic has an independent gain (0.01–100 %).
 * Optionally processes only a selected frequency band (high-pass filtered)
 * so the exciter does not muddy the low end.
 */
#include <math.h>
#include <string.h>
#include <float.h>
#include "../jdsp_header.h"

/* ---- Chebyshev polynomials of the first kind T_n(x) ---- */
static inline float T1(float x) { return x; }
static inline float T2(float x) { return 2.0f*x*x - 1.0f; }
static inline float T3(float x) { return 4.0f*x*x*x - 3.0f*x; }
static inline float T4(float x) { float x2=x*x; return 8.0f*x2*x2 - 8.0f*x2 + 1.0f; }
static inline float T5(float x) { float x2=x*x; return 16.0f*x2*x2*x - 20.0f*x2*x + 5.0f*x; }
static inline float T6(float x) { float x2=x*x, x4=x2*x2; return 32.0f*x4*x2 - 48.0f*x4 + 18.0f*x2 - 1.0f; }
static inline float T7(float x) { float x2=x*x, x4=x2*x2; return 64.0f*x4*x2*x - 112.0f*x4*x + 56.0f*x2*x - 7.0f*x; }
static inline float T8(float x) { float x2=x*x, x4=x2*x2, x6=x4*x2; return 128.0f*x6*x2 - 256.0f*x6 + 160.0f*x4 - 32.0f*x2 + 1.0f; }
static inline float T9(float x) { float x2=x*x, x4=x2*x2, x6=x4*x2; return 256.0f*x6*x2*x - 576.0f*x6*x + 432.0f*x4*x - 120.0f*x2*x + 9.0f*x; }
static inline float T10(float x) { float x2=x*x, x4=x2*x2, x6=x4*x2, x8=x4*x4; return 512.0f*x8*x2 - 1280.0f*x8 + 1120.0f*x6 - 400.0f*x4 + 50.0f*x2 - 1.0f; }

/* ---- DC blocker (1st-order high-pass) ---- */
static void dcblocker_init(HRDCBlocker *d)
{
    d->prevX = 0.0f;
    d->prevY = 0.0f;
}

static inline float dcblocker_process(HRDCBlocker *d, float x)
{
    const float R = 0.995f;
    float y = x - d->prevX + R * d->prevY;
    d->prevX = x;
    d->prevY = y;
    return y;
}

/* ---- 2nd-order Butterworth high-pass for band extraction ---- */
static void hrbiquadHP_init(HRBiquadHP *f)
{
    memset(f, 0, sizeof(*f));
}

static void hrbiquadHP_set(HRBiquadHP *f, double fs, double fc)
{
    /* 2nd-order Butterworth high-pass */
    double w0 = 2.0 * M_PI * fc / fs;
    double c = cos(w0);
    double s = sin(w0);
    double alpha = s / M_SQRT2; /* Q = 0.707 */
    double a0 = 1.0 + alpha;
    f->b0 = (float)((1.0 + c) / 2.0 / a0);
    f->b1 = (float)(-(1.0 + c) / a0);
    f->b2 = (float)((1.0 + c) / 2.0 / a0);
    f->a1 = (float)(-2.0 * c / a0);
    f->a2 = (float)((1.0 - alpha) / a0);
}

static inline float hrbiquadHP_process(HRBiquadHP *f, float x)
{
    float y = f->b0 * x + f->b1 * f->x1 + f->b2 * f->x2 - f->a1 * f->y1 - f->a2 * f->y2;
    f->x2 = f->x1; f->x1 = x;
    f->y2 = f->y1; f->y1 = y;
    return y;
}

/* ---- Parameter setup ---- */
void HarmonicExpanderSetParam(JamesDSPLib *jdsp, float harmonicGains[9],
                              float crossoverFreq, float mix)
{
    HarmonicExpander *he = &jdsp->harmonicExpander;
    for (int i = 0; i < 9; i++)
        he->harmonicGain[i] = harmonicGains[i];
    he->mix = mix;
    he->crossoverFreq = crossoverFreq;
    /* (Re)initialise band-extraction filter */
    hrbiquadHP_set(&he->bandHP[0], jdsp->fs, crossoverFreq);
    hrbiquadHP_set(&he->bandHP[1], jdsp->fs, crossoverFreq);
    dcblocker_init(&he->dc[0]);
    dcblocker_init(&he->dc[1]);
}

/* ---- Constructor / Enable / Disable ---- */
void HarmonicExpanderConstructor(JamesDSPLib *jdsp)
{
    HarmonicExpander *he = &jdsp->harmonicExpander;
    memset(he, 0, sizeof(*he));
    for (int i = 0; i < 9; i++)
        he->harmonicGain[i] = 0.0f;
    he->mix = 1.0f;
    he->crossoverFreq = 2000.0f;
    hrbiquadHP_init(&he->bandHP[0]);
    hrbiquadHP_init(&he->bandHP[1]);
    dcblocker_init(&he->dc[0]);
    dcblocker_init(&he->dc[1]);
}

void HarmonicExpanderEnable(JamesDSPLib *jdsp)
{
    jdsp->harmonicExpanderEnabled = 1;
}

void HarmonicExpanderDisable(JamesDSPLib *jdsp)
{
    jdsp->harmonicExpanderEnabled = 0;
}

/* ---- Processing ---- */
static inline float chebyshevSum(HarmonicExpander *he, float x)
{
    /* Clamp input to [-1, 1] to keep polynomials bounded */
    float xc = x;
    if (xc > 1.0f) xc = 1.0f;
    if (xc < -1.0f) xc = -1.0f;

    float y = 0.0f;
    if (he->harmonicGain[0] > 0.0f)  y += he->harmonicGain[0]  * T2(xc);
    if (he->harmonicGain[1] > 0.0f)  y += he->harmonicGain[1]  * T3(xc);
    if (he->harmonicGain[2] > 0.0f)  y += he->harmonicGain[2]  * T4(xc);
    if (he->harmonicGain[3] > 0.0f)  y += he->harmonicGain[3]  * T5(xc);
    if (he->harmonicGain[4] > 0.0f)  y += he->harmonicGain[4]  * T6(xc);
    if (he->harmonicGain[5] > 0.0f)  y += he->harmonicGain[5]  * T7(xc);
    if (he->harmonicGain[6] > 0.0f)  y += he->harmonicGain[6]  * T8(xc);
    if (he->harmonicGain[7] > 0.0f)  y += he->harmonicGain[7]  * T9(xc);
    if (he->harmonicGain[8] > 0.0f)  y += he->harmonicGain[8]  * T10(xc);
    return y;
}

void HarmonicExpanderProcess(JamesDSPLib *jdsp, size_t n)
{
    HarmonicExpander *he = &jdsp->harmonicExpander;
    float *xL = jdsp->tmpBuffer[0];
    float *xR = jdsp->tmpBuffer[1];
    for (size_t i = 0; i < n; i++)
    {
        /* Band-extract (high-pass) the signal so we excite mids/highs */
        float bandL = hrbiquadHP_process(&he->bandHP[0], xL[i]);
        float bandR = hrbiquadHP_process(&he->bandHP[1], xR[i]);

        /* Apply Chebyshev waveshaper per channel */
        float excL = chebyshevSum(he, bandL);
        float excR = chebyshevSum(he, bandR);

        /* DC-block the exciter output */
        excL = dcblocker_process(&he->dc[0], excL);
        excR = dcblocker_process(&he->dc[1], excR);

        /* Mix exciter output with original */
        xL[i] = xL[i] + excL * he->mix;
        xR[i] = xR[i] + excR * he->mix;
    }
}
