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
// #include "fft.h"

#define INIT_SAMPLE_RATE        44100
#define WINDOW_SIZE             2048
#define HOP_SIZE                512

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
        : pitch(6),
        lpf(0),
        hpf(0),
        mix(0.50f)
        {}
        
        int pitch;
        float lpf, hpf, mix;
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
    float processSampleL(float inSample);
    
    // Process Right Channel
    float processSampleR(float inSample);
    
    float window[WINDOW_SIZE], prev_win[WINDOW_SIZE], cur_win[WINDOW_SIZE],
          omega[WINDOW_SIZE/2], phi0[WINDOW_SIZE/2], synMag[WINDOW_SIZE],
          synFreq[WINDOW_SIZE], sumPhase[WINDOW_SIZE/2], anaMagn[WINDOW_SIZE], anaFreq[WINDOW_SIZE],
          cur_magnitude[WINDOW_SIZE/2], cur_phase[WINDOW_SIZE/2],
          overlap, freqPerBin, fftData[WINDOW_SIZE],
          inFIFO[WINDOW_SIZE], outFIFO[WINDOW_SIZE], win, re, im, magn, phs, prev_phs[WINDOW_SIZE/2+1],
          outData[WINDOW_SIZE*2];
    int overlap_samples, osamp, gRover, inFifoLatency, stepSize;
    
protected:
    
private:
    // STFT
    void stft(float* buf, float frameSize, float sign);
    // Pitch shifter's parameters
    Parameters parameters;
    // Pitch shifter's sample rate
    double currentSampleRate;
};


#endif /* defined(__fooHarmonizer__Shifter__) */
