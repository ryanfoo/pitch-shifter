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
    // Set Overlap Percentage - # of samples that will overlap
    overlap = (WINDOW_SIZE - HOP_SIZE) / (float)WINDOW_SIZE;
    // Set Overlap Samples - how much overlap there will be between frames
    overlap_samples = overlap * WINDOW_SIZE;
    
    // Apply blackman window to our main window
    blackman(monoData.win, WINDOW_SIZE);
    blackman(leftData.win, WINDOW_SIZE);
    blackman(rightData.win, WINDOW_SIZE);
    // Zero out previous window
    memset(monoData.pre_win, 0, WINDOW_SIZE*sizeof(float));
    memset(leftData.pre_win, 0, WINDOW_SIZE*sizeof(float));
    memset(rightData.pre_win, 0, WINDOW_SIZE*sizeof(float));
    
    // Set expected omega frequency values
    for (int i = 0; i < WINDOW_SIZE/2; i++)
    {
        monoData.om[i] = 2. * M_PI * i * osamp * (float)HOP_SIZE / (float)WINDOW_SIZE;
        leftData.om[i] = 2. * M_PI * i * osamp * (float)HOP_SIZE / (float)WINDOW_SIZE;
        rightData.om[i] = 2. * M_PI * i * osamp * (float)HOP_SIZE / (float)WINDOW_SIZE;
    }
    // Scale window for overlap add
    for (int i = 0; i < WINDOW_SIZE; i++)
    {
        monoData.win[i] *= 2. / osamp;
        leftData.win[i] *= 2. / osamp;
        rightData.win[i] *= 2. / osamp;
    }
    
    setBuffers();
}

void Shifter::setBuffers()
{
    // Zero out buffers
    memset(monoData.phi, 0, WINDOW_SIZE/2*sizeof(float));
    memset(monoData.sumPhase, 0, WINDOW_SIZE/2*sizeof(float));
    memset(leftData.phi, 0, WINDOW_SIZE/2*sizeof(float));
    memset(leftData.sumPhase, 0, WINDOW_SIZE/2*sizeof(float));
    memset(rightData.phi, 0, WINDOW_SIZE/2*sizeof(float));
    memset(rightData.sumPhase, 0, WINDOW_SIZE/2*sizeof(float));
}

# pragma mark - Mono Channel Processing -
void Shifter::processMono(float* const samples, const int numSamples)
{
    processChannel(samples, numSamples);
}

// Process Mono Data
inline void Shifter::processChannel(float* const samples, const int numSamples) noexcept
{
    // Assert that the samples are not null
    jassert (samples != nullptr);
    
    // Init vars
    long i, j, index;
    float tmp;
    
    // Set Frequencies Per Bin - # of frequencies in each bin to be analyzed = SR/WINDOW_SIZE
    freqPerBin  = currentSampleRate/(float)WINDOW_SIZE;
    
    // Init our arrays upon start-up
    if (monoStatus == false)
    {
        setBuffers();
        monoStatus = true;
    }
    
    // Process our samples
    for (i = 0; i < numSamples; i += HOP_SIZE)
    {
# pragma mark - Analysis
        // Set our incoming samples to the current stft window
        for (j = 0; j < WINDOW_SIZE; j++) monoData.cur_win[j] = samples[i+j];
        // Applies a hanning window to data
        apply_window(monoData.cur_win, monoData.win, WINDOW_SIZE);
        
        // Obtain minimum phase by shifting time domain data before taking FFT
        fftshift(monoData.cur_win, WINDOW_SIZE);
        
# pragma mark - FFT/Convert to Magnitudes + Phases
        // FFT real values (Convert to frequency domain)
        rfft(monoData.cur_win, WINDOW_SIZE/2, FFT_FORWARD);
        // Get real and imaginary #s of the FFT'd window
        complex *cbuf = (complex *)monoData.cur_win;
        
        // Get Magnitude and Phase (polar coordinates)
        for (j = 0; j < WINDOW_SIZE/2; j++) {
            monoData.cur_mag[j] = cmp_abs(cbuf[j]);
            monoData.cur_phs[j] = atan2f(cbuf[j].im, cbuf[j].re);
        }
        // Get frequencies of FFT'd signal (analysis stage)
        for (j = 0; j < WINDOW_SIZE/2; j++) {
            // Get phase difference
            tmp = monoData.cur_phs[j] - monoData.phi[j];
            monoData.phi[j] = monoData.cur_phs[j];
            
            // Subtract expected phase difference
            tmp -= monoData.om[j];
            
            // Map to +/- Pi interval
            tmp = princarg(tmp);
            
            // get deviation from bin freq from the +/- pi interval
            tmp = osamp * tmp / (2. * M_PI);
            
            // compute the k-th partials' true frequency
            tmp = (float) j * freqPerBin + tmp * freqPerBin;
            
            // Store true frequency
            monoData.anaFreq[j] = tmp;
        }
        
# pragma mark - Processing
        // Zero our processing buffers
        memset(monoData.synMagn, 0, WINDOW_SIZE*sizeof(float));
        memset(monoData.synFreq, 0, WINDOW_SIZE*sizeof(float));
        // Set new frequencies according to our pitch value
        for (j = 0; j < WINDOW_SIZE/2; j++) {
            index = j * parameters.pitch;
            if (index < WINDOW_SIZE/2) {
                monoData.synMagn[index] += monoData.cur_mag[j];
                monoData.synFreq[index] = monoData.anaFreq[j] * parameters.pitch;
            }
        }
        
# pragma mark - Synthesis
        // Write our new magnitudes and phases
        for (j = 0; j < WINDOW_SIZE/2; j++) {
            // get magnitude and true frequency from synthesis arrays
            monoData.cur_mag[j] = monoData.synMagn[j];
            
            // subtract bin mid frequency
            tmp = monoData.synFreq[j] - (float)j * freqPerBin;
            
            // get bin deviation from freq deviation
            tmp /= freqPerBin;
            
            // Factor in overlap factor
            tmp = 2. * M_PI * tmp / (float)osamp;
            
            // add the overlap phase advance back in
            tmp += monoData.om[j];
            
            // accumulate delta phase to get bin phase
            monoData.sumPhase[j] += tmp;
            monoData.cur_phs[j] = monoData.sumPhase[j];
        }
        
        // Back to Cartesian coordinates
        for (j = 0; j < WINDOW_SIZE/2; j++) {
            cbuf[j].re = monoData.cur_mag[j] * cosf(monoData.cur_phs[j]);
            cbuf[j].im = monoData.cur_mag[j] * sinf(monoData.cur_phs[j]);
        }
        
        // FFT back to time domain signal
        rfft((float*)cbuf, WINDOW_SIZE/2, FFT_INVERSE);
 
# pragma mark - Output
        // Write to output
        for (j = 0; j < HOP_SIZE; j++) {
            monoData.outData[i+j] = monoData.pre_win[j + HOP_SIZE] + monoData.cur_win[j];
        }
        
        // Filter data if filter button is on
        if (parameters.filter) processFilters(monoData.outData, HOP_SIZE);
        
        // Move previous window
        for (j = 0; j < WINDOW_SIZE; j++) {
            monoData.pre_win[j] = (j < overlap_samples) ?
            monoData.pre_win[j + HOP_SIZE] : 0;
        }
        
        // Update previous window
        for (j = 0; j < WINDOW_SIZE; j++) {
            monoData.pre_win[j] += monoData.cur_win[j];
        }
        
        // Combine input data with output data
        for (j = 0; j < HOP_SIZE; j++)
        {
            samples[i+j] = samples[i+j] * (1.0 - parameters.mix) + monoData.outData[i+j] * parameters.mix;
        }
    }
    
    // Filter data if filter button is on
    // if (parameters.filter) processFilters(samples, numSamples);
}

# pragma mark - Stereo Channel Processing -

void Shifter::processStereo(float* const left, float* const right, const int numSamples)
{
    processLeftChannel(left, numSamples);
    processRightChannel(right, numSamples);
}

// Process Mono Data
inline void Shifter::processLeftChannel(float* const samples, const int numSamples) noexcept
{
    // Assert that the samples are not null
    jassert (samples != nullptr);
    
    // Init vars
    long i, j, index;
    float tmp;
    
    // Set Frequencies Per Bin - # of frequencies in each bin to be analyzed = SR/WINDOW_SIZE
    freqPerBin  = currentSampleRate/(float)WINDOW_SIZE;
    
    // Init our arrays upon start-up
    if (stereoStatus == false)
    {
        setBuffers();
        stereoStatus = true;
    }
    
    // Process our samples
    for (i = 0; i < numSamples; i += HOP_SIZE)
    {
# pragma mark - Analysis
        // Set our incoming samples to the current stft window
        for (j = 0; j < WINDOW_SIZE; j++) leftData.cur_win[j] = samples[i+j];
            // Applies a hanning window to data
            apply_window(leftData.cur_win, leftData.win, WINDOW_SIZE);
            
            // Obtain minimum phase by shifting time domain data before taking FFT
            fftshift(leftData.cur_win, WINDOW_SIZE);
            
# pragma mark - FFT/Convert to Magnitudes + Phases
            // FFT real values (Convert to frequency domain)
            rfft(leftData.cur_win, WINDOW_SIZE/2, FFT_FORWARD);
            // Get real and imaginary #s of the FFT'd window
            complex *cbuf = (complex *)leftData.cur_win;
            
            // Get Magnitude and Phase (polar coordinates)
            for (j = 0; j < WINDOW_SIZE/2; j++) {
                leftData.cur_mag[j] = cmp_abs(cbuf[j]);
                leftData.cur_phs[j] = atan2f(cbuf[j].im, cbuf[j].re);
            }
        // Get frequencies of FFT'd signal (analysis stage)
        for (j = 0; j < WINDOW_SIZE/2; j++) {
            // Get phase difference
            tmp = leftData.cur_phs[j] - leftData.phi[j];
            leftData.phi[j] = leftData.cur_phs[j];
            
            // Subtract expected phase difference
            tmp -= leftData.om[j];
            
            // Map to +/- Pi interval
            tmp = princarg(tmp);
            
            // get deviation from bin freq from the +/- pi interval
            tmp = osamp * tmp / (2. * M_PI);
            
            // compute the k-th partials' true frequency
            tmp = (float) j * freqPerBin + tmp * freqPerBin;
            
            // Store true frequency
            leftData.anaFreq[j] = tmp;
        }
        
# pragma mark - Processing
        // Zero our processing buffers
        memset(leftData.synMagn, 0, WINDOW_SIZE*sizeof(float));
        memset(leftData.synFreq, 0, WINDOW_SIZE*sizeof(float));
        // Set new frequencies according to our pitch value
        for (j = 0; j < WINDOW_SIZE/2; j++) {
            index = j * parameters.pitch;
            if (index < WINDOW_SIZE/2) {
                leftData.synMagn[index] += leftData.cur_mag[j];
                leftData.synFreq[index] = leftData.anaFreq[j] * parameters.pitch;
            }
        }
        
# pragma mark - Synthesis
        // Write our new magnitudes and phases
        for (j = 0; j < WINDOW_SIZE/2; j++) {
            // get magnitude and true frequency from synthesis arrays
            leftData.cur_mag[j] = leftData.synMagn[j];
            
            // subtract bin mid frequency
            tmp = leftData.synFreq[j] - (float)j * freqPerBin;
            
            // get bin deviation from freq deviation
            tmp /= freqPerBin;
            
            // Factor in overlap factor
            tmp = 2. * M_PI * tmp / (float)osamp;
            
            // add the overlap phase advance back in
            tmp += om[j];
            
            // accumulate delta phase to get bin phase
            leftData.sumPhase[j] += tmp;
            leftData.cur_phs[j] = leftData.sumPhase[j];
        }
        
        // Back to Cartesian coordinates
        for (j = 0; j < WINDOW_SIZE/2; j++) {
            cbuf[j].re = leftData.cur_mag[j] * cosf(leftData.cur_phs[j]);
            cbuf[j].im = leftData.cur_mag[j] * sinf(leftData.cur_phs[j]);
        }
        
        // FFT back to time domain signal
        rfft((float*)cbuf, WINDOW_SIZE/2, FFT_INVERSE);
        
# pragma mark - Output
        // Write to output
        for (j = 0; j < HOP_SIZE; j++) {
            leftData.outData[i+j] = leftData.pre_win[j + HOP_SIZE] + leftData.cur_win[j];
        }
        
        // Filter data if filter button is on
        if (parameters.filter) processFilters(leftData.outData, HOP_SIZE);
        
        // Move previous window
        for (j = 0; j < WINDOW_SIZE; j++) {
            leftData.pre_win[j] = (j < overlap_samples) ?
            leftData.pre_win[j + HOP_SIZE] : 0;
        }
        
        // Update previous window
        for (j = 0; j < WINDOW_SIZE; j++) {
            leftData.pre_win[j] += leftData.cur_win[j];
        }
        
        // Combine input data with output data
        for (j = 0; j < HOP_SIZE; j++)
        {
            samples[i+j] = samples[i+j] * (1.0 - parameters.mix) + leftData.outData[i+j] * parameters.mix;
        }
    }
    
    // Filter data if filter button is on
    // if (parameters.filter) processFilters(samples, numSamples);
}

# pragma mark - Mono Channel Processing -
// Process Mono Data
inline void Shifter::processRightChannel(float* const samples, const int numSamples) noexcept
{
    // Assert that the samples are not null
    jassert (samples != nullptr);
    
    // Init vars
    long i, j, index;
    float tmp;
    
    // Set Frequencies Per Bin - # of frequencies in each bin to be analyzed = SR/WINDOW_SIZE
    freqPerBin  = currentSampleRate/(float)WINDOW_SIZE;
    
    // Init our arrays upon start-up
    if (stereoStatus == false)
    {
        setBuffers();
        stereoStatus = true;
    }

    // Process our samples
    for (i = 0; i < numSamples; i += HOP_SIZE)
    {
# pragma mark - Analysis
        // Set our incoming samples to the current stft window
        for (j = 0; j < WINDOW_SIZE; j++) rightData.cur_win[j] = samples[i+j];
            // Applies a hanning window to data
            apply_window(rightData.cur_win, rightData.win, WINDOW_SIZE);
            
            // Obtain minimum phase by shifting time domain data before taking FFT
            fftshift(rightData.cur_win, WINDOW_SIZE);
            
# pragma mark - FFT/Convert to Magnitudes + Phases
            // FFT real values (Convert to frequency domain)
            rfft(rightData.cur_win, WINDOW_SIZE/2, FFT_FORWARD);
            // Get real and imaginary #s of the FFT'd window
            complex *cbuf = (complex *)rightData.cur_win;
            
            // Get Magnitude and Phase (polar coordinates)
            for (j = 0; j < WINDOW_SIZE/2; j++) {
                rightData.cur_mag[j] = cmp_abs(cbuf[j]);
                rightData.cur_phs[j] = atan2f(cbuf[j].im, cbuf[j].re);
            }
        // Get frequencies of FFT'd signal (analysis stage)
        for (j = 0; j < WINDOW_SIZE/2; j++) {
            // Get phase difference
            tmp = rightData.cur_phs[j] - rightData.phi[j];
            rightData.phi[j] = rightData.cur_phs[j];
            
            // Subtract expected phase difference
            tmp -= rightData.om[j];
            
            // Map to +/- Pi interval
            tmp = princarg(tmp);
            
            // get deviation from bin freq from the +/- pi interval
            tmp = osamp * tmp / (2. * M_PI);
            
            // compute the k-th partials' true frequency
            tmp = (float) j * freqPerBin + tmp * freqPerBin;
            
            // Store true frequency
            rightData.anaFreq[j] = tmp;
        }
        
# pragma mark - Processing
        // Zero our processing buffers
        memset(rightData.synMagn, 0, WINDOW_SIZE*sizeof(float));
        memset(rightData.synFreq, 0, WINDOW_SIZE*sizeof(float));
        // Set new frequencies according to our pitch value
        for (j = 0; j < WINDOW_SIZE/2; j++) {
            index = j * parameters.pitch;
            if (index < WINDOW_SIZE/2) {
                rightData.synMagn[index] += rightData.cur_mag[j];
                rightData.synFreq[index] = rightData.anaFreq[j] * parameters.pitch;
            }
        }
        
# pragma mark - Synthesis
        // Write our new magnitudes and phases
        for (j = 0; j < WINDOW_SIZE/2; j++) {
            // get magnitude and true frequency from synthesis arrays
            rightData.cur_mag[j] = rightData.synMagn[j];
            
            // subtract bin mid frequency
            tmp = rightData.synFreq[j] - (float)j * freqPerBin;
            
            // get bin deviation from freq deviation
            tmp /= freqPerBin;
            
            // Factor in overlap factor
            tmp = 2. * M_PI * tmp / (float)osamp;
            
            // add the overlap phase advance back in
            tmp += om[j];
            
            // accumulate delta phase to get bin phase
            rightData.sumPhase[j] += tmp;
            rightData.cur_phs[j] = rightData.sumPhase[j];
        }
        
        // Back to Cartesian coordinates
        for (j = 0; j < WINDOW_SIZE/2; j++) {
            cbuf[j].re = rightData.cur_mag[j] * cosf(rightData.cur_phs[j]);
            cbuf[j].im = rightData.cur_mag[j] * sinf(rightData.cur_phs[j]);
        }
        
        // FFT back to time domain signal
        rfft((float*)cbuf, WINDOW_SIZE/2, FFT_INVERSE);
        
# pragma mark - Output
        // Write to output
        for (j = 0; j < HOP_SIZE; j++) {
            rightData.outData[i+j] = rightData.pre_win[j + HOP_SIZE] + rightData.cur_win[j];
        }
        
        // Filter data if filter button is on
        if (parameters.filter) processFilters(rightData.outData, HOP_SIZE);
        
        // Move previous window
        for (j = 0; j < WINDOW_SIZE; j++) {
            rightData.pre_win[j] = (j < overlap_samples) ?
            rightData.pre_win[j + HOP_SIZE] : 0;
        }
        
        // Update previous window
        for (j = 0; j < WINDOW_SIZE; j++) {
            rightData.pre_win[j] += rightData.cur_win[j];
        }
        
        // Combine input data with output data
        for (j = 0; j < HOP_SIZE; j++)
        {
            samples[i+j] = samples[i+j] * (1.0 - parameters.mix) + rightData.outData[i+j] * parameters.mix;
        }
    }
    
    // Filter data if filter button is on
    // if (parameters.filter) processFilters(samples, numSamples);
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
    else
    {
        hpassFilter.processSamples(samples, numSamples);
        lpassFilter.processSamples(samples, numSamples);
    }
}