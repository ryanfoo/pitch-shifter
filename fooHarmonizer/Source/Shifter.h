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
#include "fft.h"

#define INIT_SAMPLE_RATE        44100
#define WINDOW_SIZE             256               // 8192 or 16384 for plugin
#define HOP_SIZE                (WINDOW_SIZE/4)

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
        : pitch(1.0f),
        lpf(20000.0f),
        hpf(0.0f),
        mix(0.50f),
        order(0),
        filter(1)
        {}
        
        float pitch;
        float lpf;
        float hpf;
        float mix;
        int order;
        int filter;
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
    
    // Process Left Channel
    void processSampleL(float* const samples, const int numSamples) noexcept;
    
    // Process Right Channel
    void processSampleR(float* const samples, const int numSamples) noexcept;
    
    // Update Lowpass Filter's Parameters
    void updateLPFilter(void);
    
    // Updates Highpass Filter's Parameters
    void updateHPFilter(void);
    
    // Process Filters
    void processFilters(float* const samples, const int numSamples);
        
    float cur_win[WINDOW_SIZE], pre_win[WINDOW_SIZE], om[WINDOW_SIZE/2], phi[WINDOW_SIZE/2], win[WINDOW_SIZE], cur_phs[WINDOW_SIZE/2], cur_mag[WINDOW_SIZE/2], prevPhase[WINDOW_SIZE/2+1], sumPhase[WINDOW_SIZE/2+1], outData[WINDOW_SIZE*2], anaFreq[WINDOW_SIZE], anaMagn[WINDOW_SIZE], synFreq[WINDOW_SIZE], synMagn[WINDOW_SIZE];
    
    float magn, freqPerBin, expct, overlap, overlap_samples;
    
    long osamp, qpd, stepSize, frameSize;
    
    bool monoStatus = false, stereoStatus = false;
    
protected:
    
private:
    // Init vars
    void initArrays();
    // Pitch shifter's parameters
    Parameters parameters;
    // Pitch shifter's sample rate
    double currentSampleRate;
    // Filters
    IIRFilter lpassFilter, hpassFilter;
};


#endif /* defined(__fooHarmonizer__Shifter__) */
