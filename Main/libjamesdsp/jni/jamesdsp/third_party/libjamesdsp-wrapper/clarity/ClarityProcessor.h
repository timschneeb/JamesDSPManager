#pragma once

#include <cstdint>
#include <vector>

namespace clarity {

constexpr uint32_t DEFAULT_SR = 44100;

enum class Mode {
    NATURAL = 0,
    OZONE = 1,
    XHIFI = 2,
};

class IIR1 {
public:
    float b0 = 0.0f;
    float b1 = 0.0f;
    float a1 = 0.0f;
    float prevSample = 0.0f;

    void mute() { prevSample = 0.0f; }
    void setLPF_BW(float frequency, uint32_t samplingRate);
    void setHPF_BW(float frequency, uint32_t samplingRate);
    float process(float sample);
};

class NOrderBW_LH {
public:
    explicit NOrderBW_LH(uint32_t order);
    void mute();
    void setLPF(float frequency, uint32_t samplingRate);
    void setHPF(float frequency, uint32_t samplingRate);
    float process(float sample);

private:
    std::vector<IIR1> filters;
};

class NOrderBW_BP {
public:
    explicit NOrderBW_BP(uint32_t order);
    void mute();
    void setBPF(float lowCut, float highCut, uint32_t samplingRate);
    float process(float sample);

private:
    std::vector<IIR1> lowpass;
    std::vector<IIR1> highpass;
};

class WaveBuffer {
public:
    WaveBuffer(uint32_t channels, uint32_t length);
    void reset();
    float* pushZerosGetBuffer(uint32_t frames);
    void pushZeros(uint32_t frames);
    float* getBuffer() { return buffer.data(); }
    void popSamples(uint32_t frames);

private:
    std::vector<float> buffer;
    uint32_t channels = 2;
    uint32_t index = 0;
};

class NoiseSharpening {
public:
    void setSamplingRate(uint32_t samplingRate);
    void setGain(float gain);
    void setNyquistOffset(float hz);
    void reset();
    void process(float* buffer, uint32_t frames);

private:
    IIR1 filters[2];
    float prevIn[2] = {0.0f, 0.0f};
    uint32_t samplingRate = DEFAULT_SR;
    float gain = 1.0f;
    float nyquistOffsetHz = 1000.0f;
};

class HighShelf {
public:
    void setFrequency(float freq) { frequency = freq; }
    void setGainLinear(float gain);
    void setSamplingRate(uint32_t samplingRate);
    float process(float sample);

private:
    float frequency = 8250.0f;
    double gainDb = 0.0;
    double x1 = 0.0, x2 = 0.0, y1 = 0.0, y2 = 0.0;
    double b0 = 0.0, b1 = 0.0, b2 = 0.0;
    double a0 = 1.0, a1 = 0.0, a2 = 0.0;
};

class HiFi {
public:
    HiFi();
    void setSamplingRate(uint32_t samplingRate);
    void setGainLinear(float gain);
    void setLowCutHz(float hz);
    void setHighCutHz(float hz);
    void setHpMix(float mix);
    void setBpMix(float mix);
    void setBpDelayDivisor(int divisor);
    void setLpDelayDivisor(int divisor);
    void reset();
    void process(float* samples, uint32_t frames);

private:
    WaveBuffer bpBuffer;
    WaveBuffer lpBuffer;
    struct ChannelFilters {
        NOrderBW_LH lowpass;
        NOrderBW_LH highpass;
        NOrderBW_BP bandpass;
        ChannelFilters() : lowpass(1), highpass(3), bandpass(3) {}
    } filters[2];

    float gain = 1.0f;
    uint32_t samplingRate = DEFAULT_SR;
    float lowCutHz = 120.0f;
    float highCutHz = 1200.0f;
    float hpMix = 1.2f;
    float bpMix = 1.0f;
    int bpDelayDivisor = 400;
    int lpDelayDivisor = 200;
};

class ClarityProcessor {
public:
    void setSamplingRate(uint32_t samplingRate);
    void setEnabled(bool enabled);
    bool isEnabled() const { return enabled; }
    void setMode(int mode);
    void setGainLinear(float linear);
    void setPostGainDb(float db);
    void setSafety(bool enabled, float thresholdDb, float releaseMs);
    void setNaturalLpfOffsetHz(int hz);
    void setOzoneFreqHz(int hz);
    void setXhifiParams(int lowCutHz, int highCutHz, float hpMix, float bpMix, int bpDelayDivisor, int lpDelayDivisor);
    void reset();
    void process(float* samples, uint32_t frames);

private:
    void applyMode(float* samples, uint32_t frames);
    void applyPostGainAndSafety(float* samples, uint32_t frames);
    void syncFilterGain();

    NoiseSharpening natural;
    HighShelf shelf[2];
    HiFi xhifi;

    bool enabled = false;
    Mode mode = Mode::NATURAL;
    uint32_t samplingRate = DEFAULT_SR;
    float gain = 0.0f;
    float postGainLinear = 1.0f;

    bool safetyEnabled = false;
    float safetyThreshold = 0.95f;
    float safetyEnv = 0.0f;
    float safetyReleaseCoef = 0.995f;
};

} // namespace clarity
