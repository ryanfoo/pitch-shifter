//
//  Shifter.cpp
//  fooHarmonizer
//
//  Created by Ryan Foo on 4/10/15.
//
//

#include "Shifter.h"
#include "fft.h"

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
    
    float tmp, freqPerBin;
    int j, idx;
    
    // Loop through sample buffer (STFT)
    for (int i = 0; i < numSamples; i += HOP_SIZE)
    {
#pragma mark - WINDOWING
        for (j = 0; j < WINDOW_SIZE; j++)
        {
            cur_win[i] = samples[i+j];
        }
        apply_window(cur_win, window, WINDOW_SIZE);
        
#pragma mark - FFT
        fftshift(cur_win, WINDOW_SIZE);
        rfft(cur_win, WINDOW_SIZE/2, FFT_FORWARD);
        // Convert to polar form:
        complex *cur_buf = (complex *)cur_win;
        for (j = 0; j < WINDOW_SIZE/2; j++)
        {
            cur_magnitude[j] = cmp_abs(cur_buf[j]);
            cur_phase[j] = atan2f(cur_buf[j].im, cur_buf[j].re);
        }
        for (j = 0; j < WINDOW_SIZE/2; j++)
        {
            tmp = cur_phase[j] - phi0[j];
            phi0[j] = cur_phase[j];
            
            tmp -= omega[j];
            tmp = princarg(tmp);
            tmp = osamp * tmp / (2. * M_PI);
            tmp = (float)j * freqPerBin + tmp * freqPerBin;
            anaFreq[j] = tmp;
        }
        
#pragma mark - SYNTHESIS
        // Initialize arrays
        memset(synMag, 0, WINDOW_SIZE*sizeof(float));
        memset(synFreq, 0, WINDOW_SIZE*sizeof(float));
        // Write new pitch
        for (j = 0; j < WINDOW_SIZE/2; j++)
        {
            idx = j * parameters.pitch;
            if (idx < WINDOW_SIZE/2)
            {
                synMag[idx] += cur_magnitude[j];
                synFreq[idx] = anaFreq[j] * parameters.pitch;
            }
        }
        
#pragma mark - PROCESSING
        for (j = 0; j < WINDOW_SIZE/2; j++)
        {
            cur_magnitude[j] = synMag[j];
            tmp = synFreq[j] - (float)j * freqPerBin;
            tmp /= freqPerBin;
            tmp = 2. * M_PI * tmp / (float)osamp;
            tmp += omega[j];
            sumPhase[j] += tmp;
            cur_phase[j] = sumPhase[j];
        }
        // Convert back to rect form:
        for (j = 0; j < WINDOW_SIZE/2; j++)
        {
            cur_buf[j].re = cur_magnitude[j] * cosf(cur_phase[j]);
            cur_buf[j].im = cur_magnitude[j] * sinf(cur_phase[j]);
        }
        // IFFT
        rfft((float*)cur_buf, WINDOW_SIZE/2, FFT_INVERSE);
        
#pragma mark - OUTPUT
        for (j = 0; j < WINDOW_SIZE; j++)
        {
            samples[i+j] = samples[i] * (1.0 - parameters.mix) + (prev_win[j+HOP_SIZE] + cur_win[j]) * parameters.mix;
        }
        
#pragma mark - WINDOWING
        for (j = 0; j < WINDOW_SIZE; j++)
        {
            prev_win[j] = (j < overlap_samples) ? prev_win[j + HOP_SIZE] : 0;
        }
        
        for (j = 0; j < WINDOW_SIZE; j++)
        {
            prev_win[j] += cur_win[j];
        }
        
        // Interpolate data
        // if (i % 2 == 0) samples[i] = processSampleL(samples[i]);
        // else samples[i] = processSampleR(samples[i]);
    }
}

// Process Stereo Data
void Shifter::processStereo(float* const left, float* const right, const int numSamples) noexcept
{
    jassert (left != nullptr && right != nullptr);
    
    // Loop through sample buffers (STFT)
    for (int i = 0; i < numSamples; i++)
    {
        left[i] = processSampleL(left[i]);
        right[i] = processSampleR(right[i]);
    }
}

inline float Shifter::processSampleL(float inSample)
{
    float y;
    // int i;
    

    // ANALYSIS and FFT transform:
    
    // PROCESSING (actual pitchshifting)
        
    // SYNTHESIS
    
    
    return inSample*(1.0-parameters.mix) + y * parameters.mix;
}

inline float Shifter::processSampleR(float inSample)
{
    float y;
    
    return inSample*(1.0-parameters.mix) + y * parameters.mix;
}
