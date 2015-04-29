//
//  Shifter.cpp
//  fooHarmonizer
//
//  Created by Ryan Foo on 4/10/15.
//
// Phase Vocoding Method is implemented in the fooHarmonizer
// Information obtained from DAFX book ch. 7
// Pitch Shifting Process:
//      1. Input Signal (Time Domain)
//      2. Window Data
//      3. STFT Signal (to FFT data)
//      4. Convert to Magnitude and Phase form
//      5. Multiply phases/frequencies with pitch parameter
//      6. Convert from Magnitude and Phase form
//      7. Overlap Add Signals
//      8. IFFT Signal (to Time Domain)


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

void Shifter::initArrays()
{
    osamp = WINDOW_SIZE/HOP_SIZE;
    frameSize   = WINDOW_SIZE/2;
    stepSize    = WINDOW_SIZE/osamp;
    freqPerBin  = currentSampleRate/WINDOW_SIZE;
    
    hanning(win, WINDOW_SIZE);
    memset(pre_win, 0, WINDOW_SIZE*sizeof(float));
    
    overlap = (WINDOW_SIZE - HOP_SIZE) / (float)WINDOW_SIZE;
    overlap_samples = overlap * WINDOW_SIZE;
    
    for (int i = 0; i < WINDOW_SIZE/2; i++)
    {
        om[i] = 2. * M_PI * i * osamp * (float)HOP_SIZE / (float)WINDOW_SIZE;
    }
    
    for (int i = 0; i < WINDOW_SIZE; i++)
    {
        win[i] *= 2. / osamp;
    }
    
    memset(phi, 0, WINDOW_SIZE/2*sizeof(float));
    memset(sumPhase, 0, WINDOW_SIZE/2*sizeof(float));
}

// Process Mono Data
void Shifter::processMono(float* const samples, const int numSamples) noexcept
{
    jassert (samples != nullptr);
    
    int i, j, index;
    float tmp;
    
    if (monoStatus == false)
    {
        initArrays();
        monoStatus = true;
    }
    
    // Bypass pitch shifter
    if (parameters.pitch == 0)
    {
        return;
    }
    
    for (i = 0; i < numSamples; i += HOP_SIZE)
    {
        // ANALYSIS
        /* Apply window to current frame */
        for (j = 0; j < WINDOW_SIZE; j++) {
            cur_win[j] = samples[i+j];
        }
        apply_window(cur_win, win, WINDOW_SIZE);
        
        // Obtain minimum phase by shifting time domain data before taking FFT
        fftshift(cur_win, WINDOW_SIZE);
        
        /* FFT */
        rfft(cur_win, WINDOW_SIZE/2, FFT_FORWARD);
        complex *cbuf = (complex *)cur_win;
        
        /* Get Magnitude and Phase (polar coordinates) */
        for (j = 0; j < WINDOW_SIZE/2; j++) {
            cur_mag[j] = cmp_abs(cbuf[j]);
            cur_phs[j] = atan2f(cbuf[j].im, cbuf[j].re);
        }
        
        for (j = 0; j < WINDOW_SIZE/2; j++) {
            // Get phase difference
            tmp = cur_phs[j] - phi[j];
            phi[j] = cur_phs[j];
            
            // Subtract expected phase difference
            tmp -= om[j];
            
            // Map to +/- Pi interval
            tmp = princarg(tmp);
            
            // get deviation from bin freq from the +/- pi interval
            tmp = osamp * tmp / (2. * M_PI);
            
            // compute the k-th partials' true frequency
            tmp = (float) j * freqPerBin + tmp * freqPerBin;
            
            // Store true frequency
            anaFreq[j] = tmp;
        }
        
        // PROCESSING
        memset(synMagn, 0, WINDOW_SIZE*sizeof(float));
        memset(synFreq, 0, WINDOW_SIZE*sizeof(float));
        for (j = 0; j < WINDOW_SIZE/2; j++) {
            index = j * parameters.pitch;
            if (index < WINDOW_SIZE/2) {
                synMagn[index] += cur_mag[j];
                synFreq[index] = anaFreq[j] * parameters.pitch;
                // This next line should commented, for for some reason it sounds better like this
                // data->gSynFreq[index] = data->gAnaFreq[j];
            }
        }
        
        // SYNTHESIS
        for (j = 0; j < WINDOW_SIZE/2; j++) {
            /* get magnitude and true frequency from synthesis arrays */
            cur_mag[j] = synMagn[j];
            
            /* subtract bin mid frequency */
            tmp = synFreq[j] - (float)j * freqPerBin;
            
            /* get bin deviation from freq deviation */
            tmp /= freqPerBin;
            
            /* take osamp into account */
            tmp = 2. * M_PI * tmp / (float)osamp;
            
            /* add the overlap phase advance back in */
            tmp += om[j];
            
            // accumulate delta phase to get bin phase
            sumPhase[j] += tmp;
            cur_phs[j] = sumPhase[j];
        }
        
        /* Back to Cartesian coordinates */
        for (j = 0; j < WINDOW_SIZE/2; j++) {
            cbuf[j].re = cur_mag[j] * cosf(cur_phs[j]);
            cbuf[j].im = cur_mag[j] * sinf(cur_phs[j]);
        }
        
        /* Back to Time Domain */
        rfft((float*)cbuf, WINDOW_SIZE/2, FFT_INVERSE);
        
        /* Assign to the output */
        for (j = 0; j < HOP_SIZE; j++) {
            outData[i+j] = pre_win[j + HOP_SIZE] + cur_win[j];
        }
        
        /* Move previous window */
        for (j = 0; j < WINDOW_SIZE; j++) {
            pre_win[j] = (j < overlap_samples) ?
            pre_win[j + HOP_SIZE] : 0;
        }
        
        /* Update previous window */
        for (j = 0; j < WINDOW_SIZE; j++) {
            pre_win[j] += cur_win[j];
        }
        
        for (j = 0; j < HOP_SIZE; j++)
        {
            samples[i+j] = samples[i+j] * (1.0 - parameters.mix) + outData[i+j] * parameters.mix;
        }
    }
}

void Shifter::processStereo(float* const left, float* const right, const int numSamples) noexcept
{
    jassert (left != nullptr && right != nullptr);
}

// Process Left Channel Stereo Data
void Shifter::processSampleL(float* const buf, const int numSamples) noexcept
{
    jassert (buf != nullptr);
    
    float tmp;
    int i, j, index;
    
    if (stereoStatus == false)
    {
        initArrays();
        stereoStatus = true;
    }
    
    // Bypass pitch shifter
    if (parameters.pitch == 0)
    {
        // return;
    }
    else
    {
        for (i = 0; i < numSamples; i += HOP_SIZE)
        {
            // ANALYSIS
            /* Apply window to current frame */
            for (j = 0; j < WINDOW_SIZE; j++) {
                cur_win[j] = buf[i+j];
            }
            apply_window(cur_win, win, WINDOW_SIZE);
            
            // Obtain minimum phase by shifting time domain data before taking FFT
            fftshift(cur_win, WINDOW_SIZE);
            
            /* FFT */
            rfft(cur_win, WINDOW_SIZE/2, FFT_FORWARD);
            complex *cbuf = (complex *)cur_win;
            
            /* Get Magnitude and Phase (polar coordinates) */
            for (j = 0; j < WINDOW_SIZE/2; j++) {
                cur_mag[j] = cmp_abs(cbuf[j]);
                cur_phs[j] = atan2f(cbuf[j].im, cbuf[j].re);
            }
            
            for (j = 0; j < WINDOW_SIZE/2; j++) {
                // Get phase difference
                tmp = cur_phs[j] - phi[j];
                phi[j] = cur_phs[j];
                
                // Subtract expected phase difference
                tmp -= om[j];
                
                // Map to +/- Pi interval
                tmp = princarg(tmp);
                
                // get deviation from bin freq from the +/- pi interval
                tmp = osamp * tmp / (2. * M_PI);
                
                // compute the k-th partials' true frequency
                tmp = (float) j * freqPerBin + tmp * freqPerBin;
                
                // Store true frequency
                anaFreq[j] = tmp;
            }
            
            // PROCESSING
            memset(synMagn, 0, WINDOW_SIZE*sizeof(float));
            memset(synFreq, 0, WINDOW_SIZE*sizeof(float));
            for (j = 0; j < WINDOW_SIZE/2; j++) {
                index = j * parameters.pitch;
                if (index < WINDOW_SIZE/2) {
                    synMagn[index] += cur_mag[j];
                    synFreq[index] = anaFreq[j] * parameters.pitch;
                    // This next line should commented, for for some reason it sounds better like this
                    // data->gSynFreq[index] = data->gAnaFreq[j];
                }
            }
            
            // SYNTHESIS
            for (j = 0; j < WINDOW_SIZE/2; j++) {
                /* get magnitude and true frequency from synthesis arrays */
                cur_mag[j] = synMagn[j];
                
                /* subtract bin mid frequency */
                tmp = synFreq[j] - (float)j * freqPerBin;
                
                /* get bin deviation from freq deviation */
                tmp /= freqPerBin;
                
                /* take osamp into account */
                tmp = 2. * M_PI * tmp / (float)osamp;
                
                /* add the overlap phase advance back in */
                tmp += om[j];
                
                // accumulate delta phase to get bin phase
                sumPhase[j] += tmp;
                cur_phs[j] = sumPhase[j];
            }
            
            /* Back to Cartesian coordinates */
            for (j = 0; j < WINDOW_SIZE/2; j++) {
                cbuf[j].re = cur_mag[j] * cosf(cur_phs[j]);
                cbuf[j].im = cur_mag[j] * sinf(cur_phs[j]);
            }
            
            /* Back to Time Domain */
            rfft((float*)cbuf, WINDOW_SIZE/2, FFT_INVERSE);
            
            /* Assign to the output */
            for (j = 0; j < HOP_SIZE; j++) {
                outData[i+j] = pre_win[j + HOP_SIZE] + cur_win[j];
            }
            
            /* Move previous window */
            for (j = 0; j < WINDOW_SIZE; j++) {
                pre_win[j] = (j < overlap_samples) ?
                pre_win[j + HOP_SIZE] : 0;
            }
            
            /* Update previous window */
            for (j = 0; j < WINDOW_SIZE; j++) {
                pre_win[j] += cur_win[j];
            }
            
            for (j = 0; j < HOP_SIZE; j++)
            {
                buf[i+j] = buf[i+j] * (1.0 - parameters.mix) + outData[i+j] * parameters.mix;
            }
        }
    }
}

// Process Right Channel Stereo
void Shifter::processSampleR(float* const buf, const int numSamples) noexcept
{
    jassert (buf != nullptr);

    float tmp;
    int i, j, index;
    
    if (stereoStatus == false)
    {
        initArrays();
        stereoStatus = true;
    }

    // Bypass pitch shifter
    if (parameters.pitch == 0)
    {
        // return;
    }
    else
    {
        for (i = 0; i < numSamples; i += HOP_SIZE)
        {
            // ANALYSIS
            /* Apply window to current frame */
            for (j = 0; j < WINDOW_SIZE; j++) {
                cur_win[j] = buf[i+j];
            }
            apply_window(cur_win, win, WINDOW_SIZE);
            
            // Obtain minimum phase by shifting time domain data before taking FFT
            fftshift(cur_win, WINDOW_SIZE);
            
            /* FFT */
            rfft(cur_win, WINDOW_SIZE/2, FFT_FORWARD );
            complex *cbuf = (complex *)cur_win;
            
            /* Get Magnitude and Phase (polar coordinates) */
            for (j = 0; j < WINDOW_SIZE/2; j++) {
                cur_mag[j] = cmp_abs(cbuf[j]);
                cur_phs[j] = atan2f(cbuf[j].im, cbuf[j].re);
            }
            
            for (j = 0; j < WINDOW_SIZE/2; j++) {
                // Get phase difference
                tmp = cur_phs[j] - phi[j];
                phi[j] = cur_phs[j];
                
                // Subtract expected phase difference
                tmp -= om[j];
                
                // Map to +/- Pi interval
                tmp = princarg(tmp);
                
                // get deviation from bin freq from the +/- pi interval
                tmp = osamp * tmp / (2. * M_PI);
                
                // compute the k-th partials' true frequency
                tmp = (float) j * freqPerBin + tmp * freqPerBin;
                
                // Store true frequency
                anaFreq[j] = tmp;
            }
            
            // PROCESSING
            memset(synMagn, 0, WINDOW_SIZE*sizeof(float));
            memset(synFreq, 0, WINDOW_SIZE*sizeof(float));
            for (j = 0; j < WINDOW_SIZE/2; j++) {
                index = j * parameters.pitch;
                if (index < WINDOW_SIZE/2) {
                    synMagn[index] += cur_mag[j];
                    synFreq[index] = anaFreq[j] * parameters.pitch;
                    // This next line should commented, for for some reason it sounds better like this
                    // data->gSynFreq[index] = data->gAnaFreq[j];
                }
            }
            
            // SYNTHESIS
            for (j = 0; j < WINDOW_SIZE/2; j++) {
                /* get magnitude and true frequency from synthesis arrays */
                cur_mag[j] = synMagn[j];
                
                /* subtract bin mid frequency */
                tmp = synFreq[j] - (float)j * freqPerBin;
                
                /* get bin deviation from freq deviation */
                tmp /= freqPerBin;
                
                /* take osamp into account */
                tmp = 2. * M_PI * tmp / (float)osamp;
                
                /* add the overlap phase advance back in */
                tmp += om[j];
                
                // accumulate delta phase to get bin phase
                sumPhase[j] += tmp;
                cur_phs[j] = sumPhase[j];
            }
            
            /* Back to Cartesian coordinates */
            for (j = 0; j < WINDOW_SIZE/2; j++) {
                cbuf[j].re = cur_mag[j] * cosf(cur_phs[j]);
                cbuf[j].im = cur_mag[j] * sinf(cur_phs[j]);
            }
            
            /* Back to Time Domain */
            rfft((float*)cbuf, WINDOW_SIZE/2, FFT_INVERSE);
            
            /* Assign to the output */
            for (j = 0; j < HOP_SIZE; j++) {
                outData[i+j] = pre_win[j + HOP_SIZE] + cur_win[j];
            }
            
            /* Move previous window */
            for (j = 0; j < WINDOW_SIZE; j++) {
                pre_win[j] = (j < overlap_samples) ?
                pre_win[j + HOP_SIZE] : 0;
            }
            
            /* Update previous window */
            for (j = 0; j < WINDOW_SIZE; j++) {
                pre_win[j] += cur_win[j];
            }
            
            for (j = 0; j < HOP_SIZE; j++)
            {
                buf[i+j] = buf[i+j] * (1.0 - parameters.mix) + outData[i+j] * parameters.mix;
            }
        }
    }
}

void Shifter::processChannel(float* const ch, const int numSamples)
{
    int i, j;
    long tmp, index;
    
    for (i = 0; i < numSamples; i += HOP_SIZE)
    {
        // ANALYSIS
        /* Apply window to current frame */
        for (j = 0; j < WINDOW_SIZE; j++) {
            cur_win[j] =  ch[i+j];
        }
        apply_window(cur_win, win, WINDOW_SIZE);
        
        // Obtain minimum phase by shifting time domain data before taking FFT
        fftshift(cur_win, WINDOW_SIZE);
        
        /* FFT */
        rfft(cur_win, WINDOW_SIZE/2, FFT_FORWARD );
        complex *cbuf = (complex *)cur_win;
        
        /* Get Magnitude and Phase (polar coordinates) */
        for (j = 0; j < WINDOW_SIZE/2; j++) {
            cur_mag[j] = cmp_abs(cbuf[j]);
            cur_phs[j] = atan2f(cbuf[j].im, cbuf[j].re);
        }
        
        for (j = 0; j < WINDOW_SIZE/2; j++) {
            // Get phase difference
            tmp = cur_phs[j] - phi[j];
            phi[j] = cur_phs[j];
            
            // Subtract expected phase difference
            tmp -= om[j];
            
            // Map to +/- Pi interval
            tmp = princarg(tmp);
            
            // get deviation from bin freq from the +/- pi interval
            tmp = osamp * tmp / (2. * M_PI);
            
            // compute the k-th partials' true frequency
            tmp = (float) j * freqPerBin + tmp * freqPerBin;
            
            // Store true frequency
            anaFreq[j] = tmp;
        }
        
        // PROCESSING
        memset(synMagn, 0, WINDOW_SIZE*sizeof(float));
        memset(synFreq, 0, WINDOW_SIZE*sizeof(float));
        for (j = 0; j < WINDOW_SIZE/2; j++) {
            index = j * parameters.pitch;
            if (index < WINDOW_SIZE/2) {
                synMagn[index] += cur_mag[j];
                synFreq[index] = anaFreq[j] * parameters.pitch;
                // This next line should commented, for for some reason it sounds better like this
                // data->gSynFreq[index] = data->gAnaFreq[j];
            }
        }
        
        // SYNTHESIS
        for (j = 0; j < WINDOW_SIZE/2; j++) {
            /* get magnitude and true frequency from synthesis arrays */
            cur_mag[j] = synMagn[j];
            
            /* subtract bin mid frequency */
            tmp = synFreq[j] - (float)j * freqPerBin;
            
            /* get bin deviation from freq deviation */
            tmp /= freqPerBin;
            
            /* take osamp into account */
            tmp = 2. * M_PI * tmp / (float)osamp;
            
            /* add the overlap phase advance back in */
            tmp += om[j];
            
            // accumulate delta phase to get bin phase
            sumPhase[j] += tmp;
            cur_phs[j] = sumPhase[j];
        }
        
        /* Back to Cartesian coordinates */
        for (j = 0; j < WINDOW_SIZE/2; j++) {
            cbuf[j].re = cur_mag[j] * cosf(cur_phs[j]);
            cbuf[j].im = cur_mag[j] * sinf(cur_phs[j]);
        }
        
        /* Back to Time Domain */
        rfft((float*)cbuf, WINDOW_SIZE/2, FFT_INVERSE);
        
        /* Assign to the output */
        for (j = 0; j < HOP_SIZE; j++) {
            outData[i+j] = pre_win[j + HOP_SIZE] + cur_win[j];
        }
        
        /* Move previous window */
        for (j = 0; j < WINDOW_SIZE; j++) {
            pre_win[j] = (j < overlap_samples) ?
            pre_win[j + HOP_SIZE] : 0;
        }
        
        /* Update previous window */
        for (j = 0; j < WINDOW_SIZE; j++) {
            pre_win[j] += cur_win[j];
        }
        
        for (j = 0; j < HOP_SIZE; j++)
        {
            ch[i+j] = ch[i+j] * (1.0 - parameters.mix) + outData[i+j] * parameters.mix;
        }
    }
}
