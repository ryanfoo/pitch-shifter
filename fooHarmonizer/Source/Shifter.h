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
#define WINDOW_SIZE             4096               // 256 is good (powers of two...) 4096
#define HOP_SIZE                (WINDOW_SIZE/4)

class Shifter
{
public:
    // Constructor
    Shifter();
    // Deconstructor
    ~Shifter();

    typedef struct
    {
        float win[WINDOW_SIZE];
        float cur_win[WINDOW_SIZE];
        float pre_win[WINDOW_SIZE];
        float om[WINDOW_SIZE/2];
        float phi[WINDOW_SIZE/2];
        float cur_phs[WINDOW_SIZE/2];
        float cur_mag[WINDOW_SIZE/2];
        float sumPhase[WINDOW_SIZE/2];
        float outData[WINDOW_SIZE*2];
        float anaFreq[WINDOW_SIZE];
        float anaMagn[WINDOW_SIZE];
        float synFreq[WINDOW_SIZE];
        float synMagn[WINDOW_SIZE];
    } data;
    
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
    
    // Init arrays
    void initArrays();
    void setBuffers();
    
    void processMono(float* const samples, const int numSamples);
    void processStereo(float* const left, float* const right, const int numSamples);
    
    // Process Mono
    void processChannel(float* const samples, const int numSamples) noexcept;
    
    // Process Left Channel
    void processLeftChannel(float* const samples, const int numSamples) noexcept;
    
    // Process Right Channel
    void processRightChannel(float* const samples, const int numSamples) noexcept;
    
    // Update Lowpass Filter's Parameters
    void updateLPFilter(void);
    
    // Updates Highpass Filter's Parameters
    void updateHPFilter(void);
    
    // Process Filters
    void processFilters(float* const samples, const int numSamples);
    
    // Arrays to process FFT, Magnitude, Phase
    float cur_win[WINDOW_SIZE], pre_win[WINDOW_SIZE], om[WINDOW_SIZE/2], phi[WINDOW_SIZE/2], win[WINDOW_SIZE], cur_phs[WINDOW_SIZE/2], cur_mag[WINDOW_SIZE/2], sumPhase[WINDOW_SIZE/2], outData[WINDOW_SIZE*2], anaFreq[WINDOW_SIZE], anaMagn[WINDOW_SIZE], synFreq[WINDOW_SIZE], synMagn[WINDOW_SIZE];
    
    // Variables for processing FFT windows
    float magn, freqPerBin, expct, overlap, overlap_samples;
    
    // Oversampling factor
    long osamp;
    
    // For initializing arrays upon start up
    bool monoStatus = false, stereoStatus = false;
    
    data monoData, leftData, rightData;
    
protected:
    
private:
    // Pitch shifter's parameters
    Parameters parameters;
    // Pitch shifter's sample rate
    double currentSampleRate;
    // Filters
    IIRFilter lpassFilter, hpassFilter;
};


#endif /* defined(__fooHarmonizer__Shifter__) */
