//
//  Shifter.h
//  fooHarmonizer
//
//  Created by Ryan Foo on 4/10/15.
//
//

#ifndef __fooHarmonizer__Shifter__
#define __fooHarmonizer__Shifter__

// #include <stdio.h>
#include "../JuceLibraryCode/JuceHeader.h"

#define INIT_SAMPLE_RATE        44100
#define WINDOW_SIZE             8192

class Shifter
{
public:
    // Constructor
    Shifter();
    // Deconstructor
    ~Shifter();
    
    // Struct holds pitch shifter's parameters
    struct Parameters
    {
        Parameters() noexcept
        : pitch(12),
        lpf(0.0f),
        hpf(0.0f),
        mix(0.50f)
        {}
        
        int pitch;
        float lpf;
        float hpf;
        float mix;
    };
    
    // Get parameters
    const Parameters& getParameters() const noexcept
    {
        return parameters;
    };
    
    // Set parameters
    void setParameters(const Parameters& newParams);
    
    // Set sample rate
    void setSampleRate(const double sampleRate);
    
    // Clear buffers
    void prepareToPlay();
    
    // Process Mono
    void processMono(float* const samples, const int numSamples) noexcept;
    
    // Process Stereo
    void processStereo(float* const left, float* const right, const int numSamples) noexcept;
    
    // Process Left Channel
    void processSampleL();
    
    // Process Right Channel
    void processSampleR();

    float inFifoL[WINDOW_SIZE], inFifoR[WINDOW_SIZE], outFifoL[WINDOW_SIZE], outFifoR[WINDOW_SIZE],
          fftData[WINDOW_SIZE*2], prevPhase[WINDOW_SIZE/2+1], sumPhase[WINDOW_SIZE/2+1], outData[WINDOW_SIZE*2],
          anaFreq[WINDOW_SIZE], anaMagn[WINDOW_SIZE], synFreq[WINDOW_SIZE], synMagn[WINDOW_SIZE];
    
    float magn, window, re, im, freqPerBin, expct;
    
    long gRover = 0, osamp, qpd, idx, inFifoLatency, stepSize, frameSize;
    
    bool monoStatus = false, stereoStatus = false;
    
protected:
    
private:
    // Init vars
    void initArrays();
    // STFT
    void stft(float* buf, float frameSize, float sign);
    // Pitch shifter's parameters
    Parameters parameters;
    // Pitch shifter's sample rate
    double currentSampleRate;
};


#endif /* defined(__fooHarmonizer__Shifter__) */
