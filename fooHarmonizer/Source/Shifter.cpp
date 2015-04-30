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
//      7. IFFT Signal (to Time Domain)
//      8. Overlap Add Signals
//
// I use the fft.c/h written by:
//          Ge Wang (gewang@cs.princeton.edu)
//          Perry R. Cook (prc@cs.princeton.edu)
// 
// fft.c/h makes a window for the input signal and converts signals to and from FFT with rfft function
//

#include "Shifter.h"

# pragma mark - Initialization and Constructors -

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

// Initialize our Arrays upon startup
void Shifter::initArrays()
{
    // Set FFT Oversampling Factor - determines the overlap between adjacent STFT frames
    osamp = WINDOW_SIZE/HOP_SIZE;
    // Set Frequencies Per Bin - # of frequencies in each bin to be analyzed = SR/WINDOW_SIZE
    freqPerBin  = currentSampleRate/WINDOW_SIZE;
    // Apply hanning window to our main window
    hanning(win, WINDOW_SIZE);
    // Zero out previous window
    memset(pre_win, 0, WINDOW_SIZE*sizeof(float));
    // Set Overlap Percentage - # of samples that will overlap
    overlap = (WINDOW_SIZE - HOP_SIZE) / (float)WINDOW_SIZE;
    // Set Overlap Samples - how much overlap there will be between frames
    overlap_samples = overlap * WINDOW_SIZE;
    
    // Set expected omega frequency values
    for (int i = 0; i < WINDOW_SIZE/2; i++)
    {
        om[i] = 2. * M_PI * i * osamp * (float)HOP_SIZE / (float)WINDOW_SIZE;
    }
    // Scale window for overlap add
    for (int i = 0; i < WINDOW_SIZE; i++)
    {
        win[i] *= 2. / osamp;
    }
    
    // Zero out buffers
    memset(phi, 0, WINDOW_SIZE/2*sizeof(float));
    memset(sumPhase, 0, WINDOW_SIZE/2*sizeof(float));
}

# pragma mark - Mono Channel Processing -
// Process Mono Data
void Shifter::processMono(float* const samples, const int numSamples) noexcept
{
    // Assert that the samples are not null
    jassert (samples != nullptr);
    
    // Init vars
    int i, j, index;
    float tmp;
    
    // Init our arrays upon start-up
    if (monoStatus == false)
    {
        initArrays();
        monoStatus = true;
    }
    
//    // Bypass pitch shifter
//    if (parameters.pitch == 1.0)
//    {
//        samples[i] = samples[i] * 0.5;
//        return;
//    }
    
    // Process our samples
    for (i = 0; i < numSamples; i += HOP_SIZE)
    {
# pragma mark - Analysis
        // Set our incoming samples to the current stft window
        for (j = 0; j < WINDOW_SIZE; j++) cur_win[j] = samples[i+j];
        // Applies a hanning window to data
        apply_window(cur_win, win, WINDOW_SIZE);
        
        // Obtain minimum phase by shifting time domain data before taking FFT
        fftshift(cur_win, WINDOW_SIZE);
        
# pragma mark - FFT/Convert to Magnitudes + Phases
        // FFT real values (Convert to frequency domain)
        rfft(cur_win, WINDOW_SIZE/2, FFT_FORWARD);
        // Get real and imaginary #s of the FFT'd window
        complex *cbuf = (complex *)cur_win;
        
        // Get Magnitude and Phase (polar coordinates)
        for (j = 0; j < WINDOW_SIZE/2; j++) {
            cur_mag[j] = cmp_abs(cbuf[j]);
            cur_phs[j] = atan2f(cbuf[j].im, cbuf[j].re);
        }
        // Get frequencies of FFT'd signal (analysis stage)
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
        
# pragma mark - Processing
        // Zero our processing buffers
        memset(synMagn, 0, WINDOW_SIZE*sizeof(float));
        memset(synFreq, 0, WINDOW_SIZE*sizeof(float));
        // Set new frequencies according to our pitch value
        for (j = 0; j < WINDOW_SIZE/2; j++) {
            index = j * parameters.pitch;
            if (index < WINDOW_SIZE/2) {
                synMagn[index] += cur_mag[j];
                synFreq[index] = anaFreq[j] * parameters.pitch;
            }
        }
        
# pragma mark - Synthesis
        // Write our new magnitudes and phases
        for (j = 0; j < WINDOW_SIZE/2; j++) {
            // get magnitude and true frequency from synthesis arrays
            cur_mag[j] = synMagn[j];
            
            // subtract bin mid frequency
            tmp = synFreq[j] - (float)j * freqPerBin;
            
            // get bin deviation from freq deviation
            tmp /= freqPerBin;
            
            // Factor in overlap factor
            tmp = 2. * M_PI * tmp / (float)osamp;
            
            // add the overlap phase advance back in
            tmp += om[j];
            
            // accumulate delta phase to get bin phase
            sumPhase[j] += tmp;
            cur_phs[j] = sumPhase[j];
        }
        
        // Back to Cartesian coordinates
        for (j = 0; j < WINDOW_SIZE/2; j++) {
            cbuf[j].re = cur_mag[j] * cosf(cur_phs[j]);
            cbuf[j].im = cur_mag[j] * sinf(cur_phs[j]);
        }
        
        // FFT back to time domain signal
        rfft((float*)cbuf, WINDOW_SIZE/2, FFT_INVERSE);
 
# pragma mark - Output
        // Write to output
        for (j = 0; j < HOP_SIZE; j++) {
            outData[i+j] = pre_win[j + HOP_SIZE] + cur_win[j];
        }
        
        // Move previous window
        for (j = 0; j < WINDOW_SIZE; j++) {
            pre_win[j] = (j < overlap_samples) ?
            pre_win[j + HOP_SIZE] : 0;
        }
        
        // Update previous window
        for (j = 0; j < WINDOW_SIZE; j++) {
            pre_win[j] += cur_win[j];
        }
        
        // Combine input data with output data
        for (j = 0; j < HOP_SIZE; j++)
        {
            samples[i+j] = samples[i+j] * (1.0 - parameters.mix) + outData[i+j] * parameters.mix;
        }
    }
    
    // Filter data if filter button is on
    if (parameters.filter) processFilters(samples, numSamples);
}

# pragma mark - Stereo Channel Processing -

// Process Left Channel Stereo Data
void Shifter::processSampleL(float* const samples, const int numSamples) noexcept
{
    // Assert that the samples are not null
    jassert (samples != nullptr);
    
    // Init vars
    int i, j, index;
    float tmp;
    
    // Init our arrays upon start-up
    if (monoStatus == false)
    {
        initArrays();
        monoStatus = true;
    }
    
    //    // Bypass pitch shifter
    //    if (parameters.pitch == 1.0)
    //    {
    //        samples[i] = samples[i] * 0.5;
    //        return;
    //    }
    
    // Process our samples
    for (i = 0; i < numSamples; i += HOP_SIZE)
    {
# pragma mark - Analysis
        // Set our incoming samples to the current stft window
        for (j = 0; j < WINDOW_SIZE; j++) cur_win[j] = samples[i+j];
            // Applies a hanning window to data
            apply_window(cur_win, win, WINDOW_SIZE);
            
            // Obtain minimum phase by shifting time domain data before taking FFT
            fftshift(cur_win, WINDOW_SIZE);
            
# pragma mark - FFT/Convert to Magnitudes + Phases
            // FFT real values (Convert to frequency domain)
            rfft(cur_win, WINDOW_SIZE/2, FFT_FORWARD);
            // Get real and imaginary #s of the FFT'd window
            complex *cbuf = (complex *)cur_win;
            
            // Get Magnitude and Phase (polar coordinates)
            for (j = 0; j < WINDOW_SIZE/2; j++) {
                cur_mag[j] = cmp_abs(cbuf[j]);
                cur_phs[j] = atan2f(cbuf[j].im, cbuf[j].re);
            }
        // Get frequencies of FFT'd signal (analysis stage)
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
        
# pragma mark - Processing
        // Zero our processing buffers
        memset(synMagn, 0, WINDOW_SIZE*sizeof(float));
        memset(synFreq, 0, WINDOW_SIZE*sizeof(float));
        // Set new frequencies according to our pitch value
        for (j = 0; j < WINDOW_SIZE/2; j++) {
            index = j * parameters.pitch;
            if (index < WINDOW_SIZE/2) {
                synMagn[index] += cur_mag[j];
                synFreq[index] = anaFreq[j] * parameters.pitch;
            }
        }
        
# pragma mark - Synthesis
        // Write our new magnitudes and phases
        for (j = 0; j < WINDOW_SIZE/2; j++) {
            // get magnitude and true frequency from synthesis arrays
            cur_mag[j] = synMagn[j];
            
            // subtract bin mid frequency
            tmp = synFreq[j] - (float)j * freqPerBin;
            
            // get bin deviation from freq deviation
            tmp /= freqPerBin;
            
            // Factor in overlap factor
            tmp = 2. * M_PI * tmp / (float)osamp;
            
            // add the overlap phase advance back in
            tmp += om[j];
            
            // accumulate delta phase to get bin phase
            sumPhase[j] += tmp;
            cur_phs[j] = sumPhase[j];
        }
        
        // Back to Cartesian coordinates
        for (j = 0; j < WINDOW_SIZE/2; j++) {
            cbuf[j].re = cur_mag[j] * cosf(cur_phs[j]);
            cbuf[j].im = cur_mag[j] * sinf(cur_phs[j]);
        }
        
        // FFT back to time domain signal
        rfft((float*)cbuf, WINDOW_SIZE/2, FFT_INVERSE);
        
# pragma mark - Output
        // Write to output
        for (j = 0; j < HOP_SIZE; j++) {
            outData[i+j] = pre_win[j + HOP_SIZE] + cur_win[j];
        }
        
        // Move previous window
        for (j = 0; j < WINDOW_SIZE; j++) {
            pre_win[j] = (j < overlap_samples) ?
            pre_win[j + HOP_SIZE] : 0;
        }
        
        // Update previous window
        for (j = 0; j < WINDOW_SIZE; j++) {
            pre_win[j] += cur_win[j];
        }
        
        // Combine input data with output data
        for (j = 0; j < HOP_SIZE; j++)
        {
            samples[i+j] = samples[i+j] * (1.0 - parameters.mix) + outData[i+j] * parameters.mix;
        }
    }
    
    // Filter data if filter button is on
    if (parameters.filter) processFilters(samples, numSamples);
}

// Process Right Channel Stereo
void Shifter::processSampleR(float* const samples, const int numSamples) noexcept
{
    // Assert that the samples are not null
    jassert (samples != nullptr);
    
    // Init vars
    int i, j, index;
    float tmp;
    
    // Init our arrays upon start-up
    if (monoStatus == false)
    {
        initArrays();
        monoStatus = true;
    }
    
    //    // Bypass pitch shifter
    //    if (parameters.pitch == 1.0)
    //    {
    //        samples[i] = samples[i] * 0.5;
    //        return;
    //    }
    
    // Process our samples
    for (i = 0; i < numSamples; i += HOP_SIZE)
    {
# pragma mark - Analysis
        // Set our incoming samples to the current stft window
        for (j = 0; j < WINDOW_SIZE; j++) cur_win[j] = samples[i+j];
            // Applies a hanning window to data
            apply_window(cur_win, win, WINDOW_SIZE);
            
            // Obtain minimum phase by shifting time domain data before taking FFT
            fftshift(cur_win, WINDOW_SIZE);
            
# pragma mark - FFT/Convert to Magnitudes + Phases
            // FFT real values (Convert to frequency domain)
            rfft(cur_win, WINDOW_SIZE/2, FFT_FORWARD);
            // Get real and imaginary #s of the FFT'd window
            complex *cbuf = (complex *)cur_win;
            
            // Get Magnitude and Phase (polar coordinates)
            for (j = 0; j < WINDOW_SIZE/2; j++) {
                cur_mag[j] = cmp_abs(cbuf[j]);
                cur_phs[j] = atan2f(cbuf[j].im, cbuf[j].re);
            }
        // Get frequencies of FFT'd signal (analysis stage)
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
        
# pragma mark - Processing
        // Zero our processing buffers
        memset(synMagn, 0, WINDOW_SIZE*sizeof(float));
        memset(synFreq, 0, WINDOW_SIZE*sizeof(float));
        // Set new frequencies according to our pitch value
        for (j = 0; j < WINDOW_SIZE/2; j++) {
            index = j * parameters.pitch;
            if (index < WINDOW_SIZE/2) {
                synMagn[index] += cur_mag[j];
                synFreq[index] = anaFreq[j] * parameters.pitch;
            }
        }
        
# pragma mark - Synthesis
        // Write our new magnitudes and phases
        for (j = 0; j < WINDOW_SIZE/2; j++) {
            // get magnitude and true frequency from synthesis arrays
            cur_mag[j] = synMagn[j];
            
            // subtract bin mid frequency
            tmp = synFreq[j] - (float)j * freqPerBin;
            
            // get bin deviation from freq deviation
            tmp /= freqPerBin;
            
            // Factor in overlap factor
            tmp = 2. * M_PI * tmp / (float)osamp;
            
            // add the overlap phase advance back in
            tmp += om[j];
            
            // accumulate delta phase to get bin phase
            sumPhase[j] += tmp;
            cur_phs[j] = sumPhase[j];
        }
        
        // Back to Cartesian coordinates
        for (j = 0; j < WINDOW_SIZE/2; j++) {
            cbuf[j].re = cur_mag[j] * cosf(cur_phs[j]);
            cbuf[j].im = cur_mag[j] * sinf(cur_phs[j]);
        }
        
        // FFT back to time domain signal
        rfft((float*)cbuf, WINDOW_SIZE/2, FFT_INVERSE);
        
# pragma mark - Output
        // Write to output
        for (j = 0; j < HOP_SIZE; j++) {
            outData[i+j] = pre_win[j + HOP_SIZE] + cur_win[j];
        }
        
        // Move previous window
        for (j = 0; j < WINDOW_SIZE; j++) {
            pre_win[j] = (j < overlap_samples) ?
            pre_win[j + HOP_SIZE] : 0;
        }
        
        // Update previous window
        for (j = 0; j < WINDOW_SIZE; j++) {
            pre_win[j] += cur_win[j];
        }
        
        // Combine input data with output data
        for (j = 0; j < HOP_SIZE; j++)
        {
            samples[i+j] = samples[i+j] * (1.0 - parameters.mix) + outData[i+j] * parameters.mix;
        }
    }
    
    // Filter data if filter button is on
    if (parameters.filter) processFilters(samples, numSamples);
}

# pragma mark - Filter Processing

// Updates Lowpass Filter's Parameters
void Shifter::updateLPFilter(void)
{
    IIRCoefficients low_coef = IIRCoefficients::makeLowPass(currentSampleRate, parameters.lpf);
    lpassFilter.setCoefficients(low_coef);
}

// Updates Highpass Filter's Parameters
void Shifter::updateHPFilter(void)
{
    IIRCoefficients high_coef = IIRCoefficients::makeLowPass(currentSampleRate, parameters.hpf);
    hpassFilter.setCoefficients(high_coef);
}

// Apply filtering to our input data
void Shifter::processFilters(float* const samples, const int numSamples)
{
    // If LPF->HPF order selected
    if (!parameters.order)
    {
        lpassFilter.processSamples(samples, numSamples);
        hpassFilter.processSamples(samples, numSamples);
    }
    // If LPF<-HPF order selected
    else
    {
        hpassFilter.processSamples(samples, numSamples);
        lpassFilter.processSamples(samples, numSamples);
    }
}