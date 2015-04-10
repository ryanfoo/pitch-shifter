//
//  Shifter.cpp
//  fooHarmonizer
//
//  Created by Ryan Foo on 4/10/15.
//
//

#include "Shifter.h"

// Constructor
// Set Initial Parameters and sample rate (44100)
Shifter::Shifter() : currentSampleRate(INIT_SAMPLE_RATE)
{
    setParameters(Parameters());
    setSampleRate(INIT_SAMPLE_RATE);
}

// Deconstructor
Shifter::~Shifter()
{
}

// Set our new parameter settings
void Shifter::setParameters(const Parameters& newParams)
{
    parameters = newParams;
}

// Set new sample rate
void Shifter::setSampleRate(const double sampleRate)
{
    jassert(sampleRate > 0);
    currentSampleRate = sampleRate;
}

// Clear/initialize our buffers
void Shifter::prepareToPlay()
{
    
}

// Process Mono Data
void Shifter::processMono(float* const samples, const int numSamples) noexcept
{
    jassert (samples != nullptr);
    
    // Loop through sample buffer
    for (int i = 0; i < numSamples; i++)
    {
        // Interpolate data
        if (i % 2 == 0) samples[i] = processSampleL(samples[i]);
        else samples[i] = processSampleR(samples[i]);
    }
}

// Process Stereo Data
void Shifter::processStereo(float* const left, float* const right, const int numSamples) noexcept
{
    jassert (left != nullptr && right != nullptr);
    
    // Loop through sample buffers
    for (int i = 0; i < numSamples; i++)
    {
        left[i] = processSampleL(left[i]);
        right[i] = processSampleR(right[i]);
    }
}

inline float Shifter::processSampleL(float inSample)
{
    float y;
    
    return inSample*(1.0-parameters.mix) + y * parameters.mix;
}

inline float Shifter::processSampleR(float inSample)
{
    float y;
    
    return inSample*(1.0-parameters.mix) + y * parameters.mix;
}
