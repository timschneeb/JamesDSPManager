#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include "../jdsp_header.h"

#define ARB_EQ_CHANNEL_SPLIT "\n#CHANNEL_SPLIT#\n"
void ArbitraryResponseEqualizerConstructor(JamesDSPLib *jdsp)
{
	jdsp->arbMag.filterLen = InitArbitraryEq(&jdsp->arbMag.coeffGen[0], 0);
	InitArbitraryEq(&jdsp->arbMag.coeffGen[1], 0);
	FFTConvolver2x2Init(&jdsp->arbMag.convState);
	FFTConvolver2x2Init(&jdsp->arbMag.conv);
	float *kDelta = (float*)malloc(jdsp->arbMag.filterLen * sizeof(float));
	memset(kDelta, 0, jdsp->arbMag.filterLen * sizeof(float));
	kDelta[0] = 1.0f;
	FFTConvolver2x2LoadImpulseResponse(&jdsp->arbMag.convState, (unsigned int)jdsp->blockSize, kDelta, kDelta, jdsp->arbMag.filterLen);
	FFTConvolver2x2LoadImpulseResponse(&jdsp->arbMag.conv, (unsigned int)jdsp->blockSize, kDelta, kDelta, jdsp->arbMag.filterLen);
	free(kDelta);
}
void ArbitraryResponseEqualizerDestructor(JamesDSPLib *jdsp)
{
	EqNodesFree(&jdsp->arbMag.coeffGen[0]);
	EqNodesFree(&jdsp->arbMag.coeffGen[1]);
	FFTConvolver2x2Free(&jdsp->arbMag.convState);
	FFTConvolver2x2Free(&jdsp->arbMag.conv);
}
void ArbitraryResponseEqualizerStringParser(JamesDSPLib *jdsp, char *stringEq)
{
	char *leftEq = stringEq;
	char *rightEq = stringEq;
	char *leftCopy = 0;
	char *rightCopy = 0;
	char *split = strstr(stringEq, ARB_EQ_CHANNEL_SPLIT);
	if (split)
	{
		size_t leftLen = (size_t)(split - stringEq);
		const char *rightStart = split + strlen(ARB_EQ_CHANNEL_SPLIT);
		size_t rightLen = strlen(rightStart);
		leftCopy = (char*)malloc(leftLen + 1);
		rightCopy = (char*)malloc(rightLen + 1);
		memcpy(leftCopy, stringEq, leftLen);
		leftCopy[leftLen] = '\0';
		memcpy(rightCopy, rightStart, rightLen + 1);
		leftEq = leftCopy;
		rightEq = rightLen > 0 ? rightCopy : leftCopy;
	}
	ArbitraryEqString2SortedNodes(&jdsp->arbMag.coeffGen[0], leftEq);
	ArbitraryEqString2SortedNodes(&jdsp->arbMag.coeffGen[1], rightEq);
	float *eqFilL = jdsp->arbMag.coeffGen[0].GetFilter(&jdsp->arbMag.coeffGen[0], (float)jdsp->fs);
	float *eqFilR = jdsp->arbMag.coeffGen[1].GetFilter(&jdsp->arbMag.coeffGen[1], (float)jdsp->fs);
	FFTConvolver2x2RefreshImpulseResponse(&jdsp->arbMag.convState, &jdsp->arbMag.conv, eqFilL, eqFilR, jdsp->arbMag.filterLen);
	if (leftCopy)
		free(leftCopy);
	if (rightCopy)
		free(rightCopy);
}
void ArbitraryResponseEqualizerEnable(JamesDSPLib *jdsp, char enable)
{
	if (jdsp->arbMagForceRefresh)
	{
		float *eqFilL = jdsp->arbMag.coeffGen[0].GetFilter(&jdsp->arbMag.coeffGen[0], (float)jdsp->fs);
		float *eqFilR = jdsp->arbMag.coeffGen[1].GetFilter(&jdsp->arbMag.coeffGen[1], (float)jdsp->fs);
		jdsp_lock(jdsp);
		FFTConvolver2x2Free(&jdsp->arbMag.convState);
		FFTConvolver2x2Free(&jdsp->arbMag.conv);
		FFTConvolver2x2LoadImpulseResponse(&jdsp->arbMag.convState, (unsigned int)jdsp->blockSize, eqFilL, eqFilR, jdsp->arbMag.filterLen);
		FFTConvolver2x2LoadImpulseResponse(&jdsp->arbMag.conv, (unsigned int)jdsp->blockSize, eqFilL, eqFilR, jdsp->arbMag.filterLen);
		jdsp_unlock(jdsp);
		jdsp->arbMagForceRefresh = 0;
	}
	if (enable)
		jdsp->arbitraryMagEnabled = 1;
}
void ArbitraryResponseEqualizerDisable(JamesDSPLib *jdsp)
{
	jdsp->arbitraryMagEnabled = 0;
}
void ArbitraryResponseEqualizerProcess(JamesDSPLib *jdsp, size_t n)
{
	FFTConvolver2x2Process(&jdsp->arbMag.conv, jdsp->tmpBuffer[0], jdsp->tmpBuffer[1], jdsp->tmpBuffer[0], jdsp->tmpBuffer[1], (unsigned int)n);
}
