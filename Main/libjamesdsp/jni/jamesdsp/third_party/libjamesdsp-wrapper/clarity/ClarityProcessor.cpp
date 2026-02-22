#include "ClarityProcessor.h"

#include <algorithm>
#include <cmath>
#include <cstring>

namespace clarity {

static constexpr double PI = 3.14159265358979323846;

void IIR1::setLPF_BW(float frequency, uint32_t samplingRate) {
    const float omega2 = static_cast<float>(PI) * frequency / static_cast<float>(samplingRate);
    const float tanOmega2 = std::tanf(omega2);
    a1 = (1.0f - tanOmega2) / (1.0f + tanOmega2);
    b0 = tanOmega2 / (1.0f + tanOmega2);
    b1 = b0;
}

void IIR1::setHPF_BW(float frequency, uint32_t samplingRate) {
    const float omega2 = static_cast<float>(PI) * frequency / static_cast<float>(samplingRate);
    const float tanOmega2 = std::tanf(omega2);
    b0 = 1.0f / (1.0f + tanOmega2);
    b1 = -b0;
    a1 = (1.0f - tanOmega2) / (1.0f + tanOmega2);
}

float IIR1::process(float sample) {
    const float hist = sample * b1;
    sample = prevSample + sample * b0;
    prevSample = sample * a1 + hist;
    return sample;
}

NOrderBW_LH::NOrderBW_LH(uint32_t order) : filters(order) {
    mute();
}

void NOrderBW_LH::mute() {
    for (auto& f : filters) f.mute();
}

void NOrderBW_LH::setLPF(float frequency, uint32_t samplingRate) {
    for (auto& f : filters) f.setLPF_BW(frequency, samplingRate);
}

void NOrderBW_LH::setHPF(float frequency, uint32_t samplingRate) {
    for (auto& f : filters) f.setHPF_BW(frequency, samplingRate);
}

float NOrderBW_LH::process(float sample) {
    for (auto& f : filters) sample = f.process(sample);
    return sample;
}

NOrderBW_BP::NOrderBW_BP(uint32_t order) : lowpass(order), highpass(order) {
    mute();
}

void NOrderBW_BP::mute() {
    for (auto& f : lowpass) f.mute();
    for (auto& f : highpass) f.mute();
}

void NOrderBW_BP::setBPF(float lowCut, float highCut, uint32_t samplingRate) {
    for (auto& f : lowpass) f.setLPF_BW(highCut, samplingRate);
    for (auto& f : highpass) f.setHPF_BW(lowCut, samplingRate);
}

float NOrderBW_BP::process(float sample) {
    for (auto& f : lowpass) sample = f.process(sample);
    for (auto& f : highpass) sample = f.process(sample);
    return sample;
}

WaveBuffer::WaveBuffer(uint32_t channels_, uint32_t length) : buffer(length * channels_), channels(channels_) {}

void WaveBuffer::reset() { index = 0; }

float* WaveBuffer::pushZerosGetBuffer(uint32_t frames) {
    const uint32_t oldIndex = index;
    if (frames > 0) {
        const uint32_t required = channels * frames + index;
        if (required > buffer.size()) buffer.resize(required);
        std::memset(buffer.data() + index, 0, channels * frames * sizeof(float));
        index += channels * frames;
    }
    return buffer.data() + oldIndex;
}

void WaveBuffer::pushZeros(uint32_t frames) {
    if (frames == 0) return;
    const uint32_t required = channels * frames + index;
    if (required > buffer.size()) buffer.resize(required);
    std::memset(buffer.data() + index, 0, channels * frames * sizeof(float));
    index += channels * frames;
}

void WaveBuffer::popSamples(uint32_t frames) {
    const uint32_t amount = channels * frames;
    if (amount > index) return;
    index -= amount;
    std::memmove(buffer.data(), buffer.data() + amount, index * sizeof(float));
}

void NoiseSharpening::setSamplingRate(uint32_t sr) {
    samplingRate = sr;
    reset();
}

void NoiseSharpening::setGain(float g) { gain = g; }

void NoiseSharpening::setNyquistOffset(float hz) {
    nyquistOffsetHz = hz;
    reset();
}

void NoiseSharpening::reset() {
    const float cutoff = std::max(40.0f, static_cast<float>(samplingRate) * 0.5f - nyquistOffsetHz);
    for (int i = 0; i < 2; ++i) {
        filters[i].setLPF_BW(cutoff, samplingRate);
        filters[i].mute();
        prevIn[i] = 0.0f;
    }
}

void NoiseSharpening::process(float* buffer, uint32_t frames) {
    for (uint32_t i = 0; i < frames; ++i) {
        for (int ch = 0; ch < 2; ++ch) {
            const uint32_t idx = i * 2 + ch;
            const float sample = buffer[idx];
            const float prev = prevIn[ch];
            prevIn[ch] = sample;
            const float xIn = sample + (sample - prev) * gain;
            IIR1& filter = filters[ch];
            const float hist = xIn * filter.b1;
            const float out = filter.prevSample + xIn * filter.b0;
            filter.prevSample = xIn * filter.a1 + hist;
            buffer[idx] = out;
        }
    }
}

void HighShelf::setGainLinear(float gain) {
    gainDb = 20.0 * std::log10(std::max(0.000001f, gain));
}

void HighShelf::setSamplingRate(uint32_t samplingRate) {
    const double x = (2.0 * PI * frequency) / static_cast<double>(samplingRate);
    const double sinX = std::sin(x);
    const double cosX = std::cos(x);
    const double y = std::exp((gainDb * std::log(10.0)) / 40.0);

    x1 = x2 = y1 = y2 = 0.0;

    const double z = std::sqrt(y * 2.0) * sinX;
    const double a = (y - 1.0) * cosX;
    const double b = (y + 1.0) - a;
    const double c = z + b;
    const double d = (y + 1.0) * cosX;
    const double e = (y + 1.0) + a;
    const double f = (y - 1.0) - d;

    a0 = 1.0 / c;
    a1 = f * 2.0;
    a2 = b - z;
    b0 = (e + z) * y;
    b1 = -y * 2.0 * ((y - 1.0) + d);
    b2 = (e - z) * y;
}

float HighShelf::process(float sample) {
    const double out = (((x1 * b1 + sample * b0 + b2 * x2) - y1 * a1) - a2 * y2) * a0;
    y2 = y1;
    y1 = out;
    x2 = x1;
    x1 = sample;
    return static_cast<float>(out);
}

HiFi::HiFi() : bpBuffer(2, 0x800), lpBuffer(2, 0x800) {}

void HiFi::setSamplingRate(uint32_t sr) { samplingRate = sr; reset(); }
void HiFi::setGainLinear(float g) { gain = g; }
void HiFi::setLowCutHz(float hz) { lowCutHz = hz; reset(); }
void HiFi::setHighCutHz(float hz) { highCutHz = hz; reset(); }
void HiFi::setHpMix(float mix) { hpMix = mix; }
void HiFi::setBpMix(float mix) { bpMix = mix; }
void HiFi::setBpDelayDivisor(int divisor) { bpDelayDivisor = std::max(1, divisor); reset(); }
void HiFi::setLpDelayDivisor(int divisor) { lpDelayDivisor = std::max(1, divisor); reset(); }

void HiFi::reset() {
    for (int i = 0; i < 2; ++i) {
        filters[i].lowpass.setLPF(lowCutHz, samplingRate);
        filters[i].lowpass.mute();
        filters[i].highpass.setHPF(highCutHz, samplingRate);
        filters[i].highpass.mute();
        filters[i].bandpass.setBPF(lowCutHz, highCutHz, samplingRate);
        filters[i].bandpass.mute();
    }

    bpBuffer.reset();
    bpBuffer.pushZeros(static_cast<uint32_t>(samplingRate / bpDelayDivisor));
    lpBuffer.reset();
    lpBuffer.pushZeros(static_cast<uint32_t>(samplingRate / lpDelayDivisor));
}

void HiFi::process(float* samples, uint32_t frames) {
    if (frames == 0) return;

    float* bpWrite = bpBuffer.pushZerosGetBuffer(frames);
    float* lpWrite = lpBuffer.pushZerosGetBuffer(frames);

    for (uint32_t i = 0; i < frames * 2; ++i) {
        const int ch = static_cast<int>(i % 2);
        const float x = samples[i];
        const float lp = filters[ch].lowpass.process(x);
        const float hp = filters[ch].highpass.process(x);
        const float bp = filters[ch].bandpass.process(x);
        samples[i] = hp;
        lpWrite[i] = lp;
        bpWrite[i] = bp;
    }

    float* bpRead = bpBuffer.getBuffer();
    float* lpRead = lpBuffer.getBuffer();
    for (uint32_t i = 0; i < frames * 2; ++i) {
        const float hp = samples[i] * gain * hpMix;
        const float bp = bpRead[i] * gain * bpMix;
        samples[i] = hp + bp + lpRead[i];
    }

    bpBuffer.popSamples(frames);
    lpBuffer.popSamples(frames);
}

void ClarityProcessor::setSamplingRate(uint32_t sr) {
    samplingRate = sr;
    reset();
}

void ClarityProcessor::setEnabled(bool e) {
    if (e && !enabled) {
        reset();
    }
    enabled = e;
}

void ClarityProcessor::setMode(int m) {
    const auto newMode = static_cast<Mode>(std::clamp(m, 0, 2));
    if (newMode != mode) {
        mode = newMode;
        reset();
    }
}

void ClarityProcessor::setGainLinear(float linear) {
    const float nextGain = std::max(0.0f, linear);
    if (std::fabs(nextGain - gain) < 1e-7f) {
        return;
    }
    gain = nextGain;
    if (mode == Mode::OZONE) {
        reset();
    } else {
        syncFilterGain();
    }
}

void ClarityProcessor::setPostGainDb(float db) {
    postGainLinear = std::pow(10.0f, db / 20.0f);
}

void ClarityProcessor::setSafety(bool enabled_, float thresholdDb, float releaseMs) {
    safetyEnabled = enabled_;
    safetyThreshold = std::pow(10.0f, thresholdDb / 20.0f);
    const float t = std::max(0.001f, releaseMs / 1000.0f);
    safetyReleaseCoef = std::exp(-1.0f / (t * static_cast<float>(samplingRate)));
}

void ClarityProcessor::setNaturalLpfOffsetHz(int hz) {
    natural.setNyquistOffset(static_cast<float>(hz));
}

void ClarityProcessor::setOzoneFreqHz(int hz) {
    for (auto& s : shelf) s.setFrequency(static_cast<float>(hz));
    for (auto& s : shelf) s.setSamplingRate(samplingRate);
}

void ClarityProcessor::setXhifiParams(int lowCutHz, int highCutHz, float hpMix, float bpMix, int bpDelayDivisor, int lpDelayDivisor) {
    xhifi.setLowCutHz(static_cast<float>(lowCutHz));
    xhifi.setHighCutHz(static_cast<float>(highCutHz));
    xhifi.setHpMix(hpMix);
    xhifi.setBpMix(bpMix);
    xhifi.setBpDelayDivisor(bpDelayDivisor);
    xhifi.setLpDelayDivisor(lpDelayDivisor);
}

void ClarityProcessor::reset() {
    natural.setSamplingRate(samplingRate);
    natural.reset();
    syncFilterGain();
    for (auto& s : shelf) {
        s.setSamplingRate(samplingRate);
    }
    xhifi.setSamplingRate(samplingRate);
    xhifi.reset();
    safetyEnv = 0.0f;
}

void ClarityProcessor::syncFilterGain() {
    natural.setGain(gain);
    shelf[0].setGainLinear(gain + 1.0f);
    shelf[1].setGainLinear(gain + 1.0f);
    xhifi.setGainLinear(gain + 1.0f);
}

void ClarityProcessor::applyMode(float* samples, uint32_t frames) {
    switch (mode) {
        case Mode::NATURAL:
            natural.process(samples, frames);
            break;
        case Mode::OZONE:
            for (uint32_t i = 0; i < frames * 2; i += 2) {
                samples[i] = shelf[0].process(samples[i]);
                samples[i + 1] = shelf[1].process(samples[i + 1]);
            }
            break;
        case Mode::XHIFI:
            xhifi.process(samples, frames);
            break;
    }
}

void ClarityProcessor::applyPostGainAndSafety(float* samples, uint32_t frames) {
    const uint32_t count = frames * 2;
    for (uint32_t i = 0; i < count; ++i) {
        float x = samples[i] * postGainLinear;
        if (safetyEnabled) {
            const float absx = std::fabs(x);
            safetyEnv = std::max(absx, safetyEnv * safetyReleaseCoef);
            const float denom = std::max(safetyThreshold, 1e-6f);
            const float ratio = safetyEnv > denom ? (denom / safetyEnv) : 1.0f;
            x *= ratio;
        }
        samples[i] = x;
    }
}

void ClarityProcessor::process(float* samples, uint32_t frames) {
    if (!enabled || samples == nullptr || frames == 0) return;
    applyMode(samples, frames);
    applyPostGainAndSafety(samples, frames);
}

} // namespace clarity
