/*
 * peq_loudness.h — Pure-C PEQ biquad cascade and loudness correction
 * for the system JamesDSP audio effect (root/Magisk version).
 *
 * PEQ: RBJ biquad cascade with per-band L/R channel routing, preamp.
 * Loudness: Fletcher-Munson compensation with low/high shelves.
 *
 * Parameter IDs (AudioEffect API):
 *   1214 — PEQ enable/disable (short)
 *   1300 — PEQ config (float array, variable size: 2 + N*5 floats, max 322 = 1288 bytes for 64 bands)
 *   1215 — Loudness enable/disable (short)
 *   1301 — Loudness config (6 floats = 24 bytes: sr, refLevel, refOffset,
 *          attenuation, volume, mode; legacy apps may send 5 floats = 20
 *          bytes, mode defaults to 0)
 *   1302 — Loudness volume update (1 float = 4 bytes)
 *   1303 — Loudness tuning (8 floats = 32 bytes)
 *   1304 — Loudness subsonic filter (4 floats = 16 bytes)
 */

#ifndef PEQ_LOUDNESS_H
#define PEQ_LOUDNESS_H

#include <stdint.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ---- Biquad ---- */

typedef enum {
    PEQ_TYPE_PEAKING    = 0,
    PEQ_TYPE_LOW_SHELF  = 1,
    PEQ_TYPE_HIGH_SHELF = 2,
    PEQ_TYPE_LOW_PASS   = 3,
    PEQ_TYPE_HIGH_PASS  = 4,
    PEQ_TYPE_BAND_PASS  = 5,
    PEQ_TYPE_NOTCH      = 6,
    PEQ_TYPE_ALL_PASS   = 7,
    PEQ_TYPE_PREAMP     = 8, /* Flat frequency-independent gain (no biquad) */
} peq_filter_type_t;

typedef enum {
    PEQ_CHAN_BOTH = 0,
    PEQ_CHAN_LEFT = 1,
    PEQ_CHAN_RIGHT = 2,
} peq_channel_mode_t;

#define PEQ_MAX_BANDS 64

typedef struct {
    /* Normalized coefficients (a0 normalized to 1) */
    double b0, b1, b2;  /* numerator */
    double a1, a2;      /* denominator (a0 = 1) */
    /* State */
    double x1, x2, y1, y2;
} biquad_t;

/* Initialize a biquad with RBJ coefficients.
 * isBandwidthOrS: false = Q mode, true = S/BW mode (S for shelves, BW for others)
 */
void biquad_init(biquad_t *bq, peq_filter_type_t type, double dbGain,
                 double freq, double srate, double bandwidthOrQOrS,
                 bool isBandwidthOrS);

/* Process a single sample through the biquad */
static inline double biquad_process(biquad_t *bq, double sample)
{
    double result = bq->b0 * sample + bq->b1 * bq->x1 + bq->b2 * bq->x2
                   - bq->a1 * bq->y1 - bq->a2 * bq->y2;
    bq->x2 = bq->x1;
    bq->x1 = sample;
    bq->y2 = bq->y1;
    bq->y1 = result;
    return result;
}

/* Clear biquad state */
static inline void biquad_reset(biquad_t *bq)
{
    bq->x1 = bq->x2 = bq->y1 = bq->y2 = 0.0;
}

/* Remove denormals from biquad state */
static inline void biquad_remove_denormals(biquad_t *bq)
{
    const double denorm_limit = 2.2250738585072014e-308;
    if (bq->x1 > -denorm_limit && bq->x1 < denorm_limit) bq->x1 = 0.0;
    if (bq->x2 > -denorm_limit && bq->x2 < denorm_limit) bq->x2 = 0.0;
    if (bq->y1 > -denorm_limit && bq->y1 < denorm_limit) bq->y1 = 0.0;
    if (bq->y2 > -denorm_limit && bq->y2 < denorm_limit) bq->y2 = 0.0;
}

/* ---- PEQ cascade ---- */

typedef struct {
    biquad_t bqL;
    biquad_t bqR;
    int channelMode;  /* PEQ_CHAN_BOTH / LEFT / RIGHT */
    bool enabled;
    bool isPreamp;       /* true for PEQ_TYPE_PREAMP (flat gain, no biquad) */
    double preampLinear; /* linear gain = 10^(gain/20) for preamp bands */
} peq_band_t;

typedef struct {
    peq_band_t bands[PEQ_MAX_BANDS];
    int bandCount;
    double sampleRate;
    double preampLinear;
    bool enabled;
} peq_cascade_t;

void peq_cascade_init(peq_cascade_t *peq);
void peq_cascade_configure(peq_cascade_t *peq, double sampleRate,
                           const float *data, int numFloats);
void peq_cascade_set_enabled(peq_cascade_t *peq, bool enabled);
void peq_cascade_process_interleaved(peq_cascade_t *peq, float *data, size_t numFrames);

/* ---- Loudness correction ---- */

typedef struct {
    /* Parameters */
    double sampleRate;
    double referenceLevel;
    double referenceOffset;
    double attenuation;
    double volumeDb;
    double lastComputedVolume;
    bool coeffsDirty;
    bool enabled;
    int mode; /* 0 = classic (Fletcher-Munson), 1 = ISO 226:2023 */
    /* Tuning parameters (set via cmd 1303) */
    double lsFreq, lsSlope, lsRatio;
    double hsFreq, hsSlope, hsRatio;
    double isoBasePhon, isoQ;
    /* Subsonic filter parameters (set via cmd 1304) */
    bool subsonicEnabled;
    double subsonicFreq;
    int subsonicOrder;
    double subsonicQ;
#define MAX_SUBSONIC_STAGES 4
    biquad_t subsonicL[MAX_SUBSONIC_STAGES];
    biquad_t subsonicR[MAX_SUBSONIC_STAGES];
    /* Internal */
    double attFactor;
    bool neutral;
    biquad_t lowShelfL, lowShelfR;
    biquad_t highShelfL, highShelfR;
    /* ISO 226:2023 mode (29-band peaking EQ cascade) */
#define ISO226_NUM_BANDS 29
    biquad_t isoBandsL[ISO226_NUM_BANDS];
    biquad_t isoBandsR[ISO226_NUM_BANDS];
    double isoGains[ISO226_NUM_BANDS];
    double isoPreampLinear;
    int isoActiveBandCount;
} loudness_t;

void loudness_init(loudness_t *l);
void loudness_configure(loudness_t *l, double sampleRate,
                        double refLevel, double refOffset,
                        double attenuation, double volumeDb, int mode);
void loudness_set_tuning(loudness_t *l,
                         double lsFreq, double lsSlope, double lsRatio,
                         double hsFreq, double hsSlope, double hsRatio,
                         double isoBasePhon, double isoQ);
void loudness_set_subsonic(loudness_t *l, bool enable, double freq,
                           int order, double qFactor);
void loudness_set_volume(loudness_t *l, double volumeDb);
void loudness_set_enabled(loudness_t *l, bool enabled);
void loudness_process_interleaved(loudness_t *l, float *data, size_t numFrames);

#ifdef __cplusplus
}
#endif

#endif /* PEQ_LOUDNESS_H */
