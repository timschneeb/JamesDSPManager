extern "C" {
#include "../jdsp_header.h"
}

#include <cstdint>
#include <cstdlib>
#include <new>

#include "ClarityProcessor.h"

namespace {

inline clarity::ClarityProcessor* getClarity(JamesDSPLib* jdsp) {
    return reinterpret_cast<clarity::ClarityProcessor*>(jdsp->clarityProcessor);
}

inline bool reserveInterleavedBuffer(JamesDSPLib* jdsp, size_t frames) {
    const size_t required = frames * 2;
    if (required == 0) {
        return true;
    }
    if (jdsp->clarityInterleavedBuffer != nullptr && jdsp->clarityInterleavedCapacity >= required) {
        return true;
    }

    float* resized = static_cast<float*>(realloc(jdsp->clarityInterleavedBuffer, required * sizeof(float)));
    if (resized == nullptr) {
        return false;
    }

    jdsp->clarityInterleavedBuffer = resized;
    jdsp->clarityInterleavedCapacity = required;
    return true;
}

inline bool ensureInterleavedBuffer(JamesDSPLib* jdsp, size_t frames) {
    const size_t required = frames * 2;
    return jdsp->clarityInterleavedBuffer != nullptr && jdsp->clarityInterleavedCapacity >= required;
}

} // namespace

extern "C" void ClarityReserveBuffer(JamesDSPLib* jdsp, size_t frames)
{
    if (jdsp == nullptr) {
        return;
    }

    jdsp_lock(jdsp);
    (void)reserveInterleavedBuffer(jdsp, frames);
    jdsp_unlock(jdsp);
}

extern "C" void ClarityConstructor(JamesDSPLib* jdsp)
{
    if (jdsp == nullptr) {
        return;
    }

    jdsp->clarityInterleavedBuffer = nullptr;
    jdsp->clarityInterleavedCapacity = 0;
    jdsp->clarityEnabled = 0;
    jdsp->clarityProcessor = new (std::nothrow) clarity::ClarityProcessor();
    if (jdsp->clarityProcessor == nullptr) {
        return;
    }

    auto* clarity = getClarity(jdsp);
    clarity->setSamplingRate(static_cast<uint32_t>(jdsp->fs));
    clarity->setEnabled(false);
    (void)reserveInterleavedBuffer(jdsp, jdsp->blockSizeMax);
}

extern "C" void ClarityDestructor(JamesDSPLib* jdsp)
{
    if (jdsp == nullptr) {
        return;
    }

    auto* clarity = getClarity(jdsp);
    delete clarity;
    jdsp->clarityProcessor = nullptr;

    if (jdsp->clarityInterleavedBuffer != nullptr) {
        free(jdsp->clarityInterleavedBuffer);
        jdsp->clarityInterleavedBuffer = nullptr;
    }
    jdsp->clarityInterleavedCapacity = 0;
}

extern "C" void ClaritySetSampleRate(JamesDSPLib* jdsp)
{
    if (jdsp == nullptr) {
        return;
    }

    jdsp_lock(jdsp);
    auto* clarity = getClarity(jdsp);
    if (clarity == nullptr) {
        jdsp_unlock(jdsp);
        return;
    }

    clarity->setSamplingRate(static_cast<uint32_t>(jdsp->fs));
    (void)reserveInterleavedBuffer(jdsp, jdsp->blockSizeMax);
    jdsp_unlock(jdsp);
}

extern "C" void ClaritySetParam(
    JamesDSPLib* jdsp,
    int mode,
    float gain,
    float postGainDb,
    char safetyEnabled,
    float safetyThresholdDb,
    float safetyReleaseMs,
    int naturalLpfOffsetHz,
    int ozoneFreqHz,
    int xhifiLowCutHz,
    int xhifiHighCutHz,
    float xhifiHpMix,
    float xhifiBpMix,
    int xhifiBpDelayDivisor,
    int xhifiLpDelayDivisor)
{
    if (jdsp == nullptr) {
        return;
    }

    jdsp_lock(jdsp);
    auto* clarity = getClarity(jdsp);
    if (clarity == nullptr) {
        jdsp_unlock(jdsp);
        return;
    }

    clarity->setMode(mode);
    clarity->setGainLinear(gain);
    clarity->setPostGainDb(postGainDb);
    clarity->setSafety(safetyEnabled != 0, safetyThresholdDb, safetyReleaseMs);
    clarity->setNaturalLpfOffsetHz(naturalLpfOffsetHz);
    clarity->setOzoneFreqHz(ozoneFreqHz);
    clarity->setXhifiParams(
        xhifiLowCutHz,
        xhifiHighCutHz,
        xhifiHpMix,
        xhifiBpMix,
        xhifiBpDelayDivisor,
        xhifiLpDelayDivisor
    );
    jdsp_unlock(jdsp);
}

extern "C" void ClarityEnable(JamesDSPLib* jdsp)
{
    if (jdsp == nullptr) {
        return;
    }

    jdsp_lock(jdsp);
    auto* clarity = getClarity(jdsp);
    if (clarity != nullptr) {
        clarity->setEnabled(true);
        jdsp->clarityEnabled = 1;
    } else {
        jdsp->clarityEnabled = 0;
    }
    jdsp_unlock(jdsp);
}

extern "C" void ClarityDisable(JamesDSPLib* jdsp)
{
    if (jdsp == nullptr) {
        return;
    }

    jdsp_lock(jdsp);
    auto* clarity = getClarity(jdsp);
    if (clarity != nullptr) {
        clarity->setEnabled(false);
    }
    jdsp->clarityEnabled = 0;
    jdsp_unlock(jdsp);
}

extern "C" void ClarityProcess(JamesDSPLib* jdsp, size_t n)
{
    if (jdsp == nullptr) {
        return;
    }

    auto* clarity = getClarity(jdsp);
    if (clarity == nullptr || n == 0) {
        return;
    }

    if (!ensureInterleavedBuffer(jdsp, n)) {
        return;
    }
    float* interleaved = jdsp->clarityInterleavedBuffer;

    for (size_t i = 0; i < n; ++i) {
        interleaved[i * 2] = jdsp->tmpBuffer[0][i];
        interleaved[i * 2 + 1] = jdsp->tmpBuffer[1][i];
    }

    clarity->process(interleaved, static_cast<uint32_t>(n));

    for (size_t i = 0; i < n; ++i) {
        jdsp->tmpBuffer[0][i] = interleaved[i * 2];
        jdsp->tmpBuffer[1][i] = interleaved[i * 2 + 1];
    }
}
