#include "colourUtil.hpp"
#include "spectrum.hpp"

#include <iostream>
#include <math.h>

using std::cout;
using std::endl;

using namespace spectrum;

static float X[61];
static float Y[61];
static float Z[61];

static const int nSamples = 60;

static void init()
{
    static bool called = false;

    if(!called)
    {
        for(int i = 0; i < nSamples; ++i)
        {
            float wl0 = lerp(float(i) / float(nSamples),
                    sampledLambdaStart, sampledLambdaEnd);
            float wl1 = lerp(float(i + 1) / float(nSamples),
                    sampledLambdaStart, sampledLambdaEnd);
            X[i] = AverageSpectrumSamples(CIE_lambda, CIE_X, nCIESamples, wl0,
                    wl1);
            Y[i] = AverageSpectrumSamples(CIE_lambda, CIE_Y, nCIESamples, wl0,
                    wl1);
            Z[i] = AverageSpectrumSamples(CIE_lambda, CIE_Z, nCIESamples, wl0,
                    wl1);
        }
    }

    called = true;
}

void ToXYZ(int n, float *vals, float xyz[3]) 
{
    xyz[0] = xyz[1] = xyz[2] = 0.f;
    
    for (int i = 0; i < n; ++i) {
        xyz[0] += X[i] * vals[i];
        xyz[1] += Y[i] * vals[i];
        xyz[2] += Z[i] * vals[i];
    }

    float scale = (float)(sampledLambdaEnd - sampledLambdaStart) / float(CIE_Y_integral * n);
    xyz[0] *= scale;
    xyz[1] *= scale;
    xyz[2] *= scale;
}

// Calls functions from PBRT
unsigned int getBlackbodyRGB(float temp)
{
    init();

    int n = nSamples;
    float wavelengths[n];

    for(int i = 0; i < n; ++i)
    {
        float wl0 = lerp(float(i) / float(nSamples),
                sampledLambdaStart, sampledLambdaEnd);
        float wl1 = lerp(float(i + 1) / float(nSamples),
                sampledLambdaStart, sampledLambdaEnd);
        wavelengths[i] = (wl0 + wl1)/2.0f;
    }

    float wl[301];
    for(int i = 0; i < 301; ++i)
    {
        wl[i] = i + sampledLambdaStart;
    }
    float blackbodyVals[301];
    BlackbodyNormalized(wl, 301, temp, blackbodyVals);

    float vals[n];
    for(int i = 0; i < n; ++i)
    {
        float wl0 = lerp(float(i) / float(n),
                sampledLambdaStart, sampledLambdaEnd);
        float wl1 = lerp(float(i + 1) / float(n),
                sampledLambdaStart, sampledLambdaEnd);
        vals[i] = AverageSpectrumSamples(wl, blackbodyVals, 301, wl0,
                wl1);
    }
    
    float xyz[3];
    ToXYZ(n, vals, xyz);

    float rgb[3];

    XYZToRGB(xyz, rgb);
    if(rgb[2] < 0)
    {
        rgb[2] = 0;
    }

    double gamma = 1.0/2.2;
    
    int red = std::min((uint)(255 << 16), ((uint)((pow(rgb[0], gamma) * 255)) << 16));
    int green = ((int)(pow(rgb[1], gamma) * 255) << 8);
    int blue = (int)(pow(rgb[2], gamma) * 255);

    return red + green + blue;
}
