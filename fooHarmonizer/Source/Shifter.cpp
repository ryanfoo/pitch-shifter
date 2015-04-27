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
    osamp = win_size/hop_size;
    frameSize   = win_size/2;
    stepSize    = win_size/osamp;
    freqPerBin  = currentSampleRate/win_size;
    
    hanning(win, win_size);
    memset(pre_win, 0, win_size*sizeof(float));
    
    overlap = (win_size - hop_size) / (float)win_size;
    overlap_samples = overlap * win_size;
    
    for (int i = 0; i < win_size/2; i++)
    {
        om[i] = 2. * M_PI * i * osamp * (float)hop_size / (float)win_size;
    }
    
    for (int i = 0; i < win_size; i++)
    {
        win[i] *= 2. / osamp;
    }
    
    memset(phi, 0, win_size/2*sizeof(float));
    memset(sumPhase, 0, win_size/2*sizeof(float));
}

/*
 * STFT - short time fourier transform
 * Variables -
 *      buf: Incoming Signal (either time-domain or frequency-domain)
 *      frameSize: Window Size of FFT analysis
 *      sign: Indicates FFT (-1) or IFFT (+1)
 */
void Shifter::stft(float* buf, float frameSize, float sign)
{
    float wr, wi, arg, *p1, *p2, temp;
    float tr, ti, ur, ui, *p1r, *p1i, *p2r, *p2i;
    long i, bitm, j, le, le2, k;
    
    for (i = 2; i < 2*WINDOW_SIZE-2; i += 2) {
        for (bitm = 2, j = 0; bitm < 2*WINDOW_SIZE; bitm <<= 1) {
            if (i & bitm) j++;
            j <<= 1;
        }
        if (i < j) {
            p1 = buf+i; p2 = buf+j;
            temp = *p1; *(p1++) = *p2;
            *(p2++) = temp; temp = *p1;
            *p1 = *p2; *p2 = temp;
        }
    }
    for (k = 0, le = 2; k < (long)(log(WINDOW_SIZE)/log(2.)+.5); k++) {
        le <<= 1;
        le2 = le>>1;
        ur = 1.0;
        ui = 0.0;
        arg = M_PI / (le2>>1);
        wr = cos(arg);
        wi = sign*sin(arg);
        for (j = 0; j < le2; j += 2) {
            p1r = buf+j; p1i = p1r+1;
            p2r = p1r+le2; p2i = p2r+1;
            for (i = j; i < 2*WINDOW_SIZE; i += le) {
                tr = *p2r * ur - *p2i * ui;
                ti = *p2r * ui + *p2i * ur;
                *p2r = *p1r - tr; *p2i = *p1i - ti;
                *p1r += tr; *p1i += ti;
                p1r += le; p1i += le;
                p2r += le; p2i += le;
            }
            tr = ur*wr - ui*wi;
            ui = ur*wi + ui*wr;
            ur = tr;
        }
    }
}

// Process Mono Data
void Shifter::processMono(float* const samples, const int numSamples) noexcept
{
    jassert (samples != nullptr);

    win_size = (float)numSamples/4;
    hop_size = win_size/4;
    
    int i, j, index;
    float tmp;
    
    if (monoStatus == false)
    {
        initArrays();
        monoStatus = true;
    }
    
    for (i = 0; i < numSamples; i += hop_size)
    {
        // ANALYSIS
        /* Apply window to current frame */
        for (j = 0; j < win_size; j++) {
            cur_win[j] = samples[i+j];
        }
        apply_window(cur_win, win, win_size);
        
        // Obtain minimum phase by shifting time domain data before taking FFT
        fftshift(cur_win, win_size);
        
        /* FFT */
        rfft(cur_win, win_size/2, FFT_FORWARD);
        complex *cbuf = (complex *)cur_win;
        
        /* Get Magnitude and Phase (polar coordinates) */
        for (j = 0; j < win_size/2; j++) {
            cur_mag[j] = cmp_abs(cbuf[j]);
            cur_phs[j] = atan2f(cbuf[j].im, cbuf[j].re);
        }
        
        for (j = 0; j < win_size/2; j++) {
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
        memset(synMagn, 0, win_size*sizeof(float));
        memset(synFreq, 0, win_size*sizeof(float));
        for (j = 0; j < win_size/2; j++) {
            index = j * parameters.pitch;
            if (index < win_size/2) {
                synMagn[index] += cur_mag[j];
                synFreq[index] = anaFreq[j] * parameters.pitch;
                // This next line should commented, for for some reason it sounds better like this
                // data->gSynFreq[index] = data->gAnaFreq[j];
            }
        }
        
        // SYNTHESIS
        for (j = 0; j < win_size/2; j++) {
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
        for (j = 0; j < win_size/2; j++) {
            cbuf[j].re = cur_mag[j] * cosf(cur_phs[j]);
            cbuf[j].im = cur_mag[j] * sinf(cur_phs[j]);
        }
        
        /* Back to Time Domain */
        rfft((float*)cbuf, win_size/2, FFT_INVERSE);
        
        /* Assign to the output */
        for (j = 0; j < hop_size; j++) {
            outData[i+j] = pre_win[j + hop_size] + cur_win[j];
        }
        
        /* Move previous window */
        for (j = 0; j < win_size; j++) {
            pre_win[j] = (j < overlap_samples) ?
            pre_win[j + hop_size] : 0;
        }
        
        /* Update previous window */
        for (j = 0; j < win_size; j++) {
            pre_win[j] += cur_win[j];
        }
        
        for (j = 0; j < hop_size; j++)
        {
            samples[i+j] = samples[i+j] * (1.0 - parameters.mix) + outData[i+j] * parameters.mix;
        }
    }
}

// Process Stereo Data
void Shifter::processStereo(float* const left, float* const right, const int numSamples) noexcept
{
    jassert (left != nullptr && right != nullptr);
    
    win_size = (float)numSamples/4;
    hop_size = win_size/4;
    
    // float tmpL, tmpR;
    
    if (stereoStatus == false)
    {
        initArrays();
        stereoStatus = true;
    }
}

void Shifter::processSampleL(float* const buf, const int numSamples) noexcept
{
    jassert (buf != nullptr);
    
    win_size = numSamples/4;
    hop_size = win_size/4;
    
    float tmp;
    int i, j, index;
    
    if (stereoStatus == false)
    {
        initArrays();
        stereoStatus = true;
    }
    
    for (i = 0; i < numSamples; i += hop_size)
    {
        // ANALYSIS
        /* Apply window to current frame */
        for (j = 0; j < win_size; j++) {
            cur_win[j] = buf[i+j];
        }
        apply_window(cur_win, win, win_size);
        
        // Obtain minimum phase by shifting time domain data before taking FFT
        fftshift(cur_win, win_size);
        
        /* FFT */
        rfft(cur_win, win_size/2, FFT_FORWARD);
        complex *cbuf = (complex *)cur_win;
        
        /* Get Magnitude and Phase (polar coordinates) */
        for (j = 0; j < win_size/2; j++) {
            cur_mag[j] = cmp_abs(cbuf[j]);
            cur_phs[j] = atan2f(cbuf[j].im, cbuf[j].re);
        }
        
        for (j = 0; j < win_size/2; j++) {
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
        memset(synMagn, 0, win_size*sizeof(float));
        memset(synFreq, 0, win_size*sizeof(float));
        for (j = 0; j < win_size/2; j++) {
            index = j * parameters.pitch;
            if (index < win_size/2) {
                synMagn[index] += cur_mag[j];
                synFreq[index] = anaFreq[j] * parameters.pitch;
                // This next line should commented, for for some reason it sounds better like this
                // data->gSynFreq[index] = data->gAnaFreq[j];
            }
        }
        
        // SYNTHESIS
        for (j = 0; j < win_size/2; j++) {
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
        for (j = 0; j < win_size/2; j++) {
            cbuf[j].re = cur_mag[j] * cosf(cur_phs[j]);
            cbuf[j].im = cur_mag[j] * sinf(cur_phs[j]);
        }
        
        /* Back to Time Domain */
        rfft((float*)cbuf, win_size/2, FFT_INVERSE);
        
        /* Assign to the output */
        for (j = 0; j < hop_size; j++) {
            outData[i+j] = pre_win[j + hop_size] + cur_win[j];
        }
        
        /* Move previous window */
        for (j = 0; j < win_size; j++) {
            pre_win[j] = (j < overlap_samples) ?
            pre_win[j + hop_size] : 0;
        }
        
        /* Update previous window */
        for (j = 0; j < win_size; j++) {
            pre_win[j] += cur_win[j];
        }
        
        for (j = 0; j < hop_size; j++)
        {
            buf[i+j] = buf[i+j] * (1.0 - parameters.mix) + outData[i+j] * parameters.mix;
        }
    }
}

void Shifter::processSampleR(float* const buf, const int numSamples) noexcept
{
    jassert (buf != nullptr);
    
    win_size = numSamples/4;
    hop_size = win_size/4;
    overlap = (win_size - hop_size) / (float)win_size;
    overlap_samples = overlap * win_size;
    
    float tmp;
    int i, j, index;
    
    if (stereoStatus == false)
    {
        initArrays();
        stereoStatus = true;
    }
    
    for (i = 0; i < numSamples; i += hop_size)
    {
        // ANALYSIS
        /* Apply window to current frame */
        for (j = 0; j < win_size; j++) {
            cur_win[j] = buf[i+j];
        }
        apply_window(cur_win, win, win_size);
        
        // Obtain minimum phase by shifting time domain data before taking FFT
        fftshift(cur_win, win_size);
        
        /* FFT */
        rfft(cur_win, win_size/2, FFT_FORWARD );
        complex *cbuf = (complex *)cur_win;
        
        /* Get Magnitude and Phase (polar coordinates) */
        for (j = 0; j < win_size/2; j++) {
            cur_mag[j] = cmp_abs(cbuf[j]);
            cur_phs[j] = atan2f(cbuf[j].im, cbuf[j].re);
        }
        
        for (j = 0; j < win_size/2; j++) {
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
        memset(synMagn, 0, win_size*sizeof(float));
        memset(synFreq, 0, win_size*sizeof(float));
        for (j = 0; j < win_size/2; j++) {
            index = j * parameters.pitch;
            if (index < win_size/2) {
                synMagn[index] += cur_mag[j];
                synFreq[index] = anaFreq[j] * parameters.pitch;
                // This next line should commented, for for some reason it sounds better like this
                // data->gSynFreq[index] = data->gAnaFreq[j];
            }
        }
        
        // SYNTHESIS
        for (j = 0; j < win_size/2; j++) {
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
        for (j = 0; j < win_size/2; j++) {
            cbuf[j].re = cur_mag[j] * cosf(cur_phs[j]);
            cbuf[j].im = cur_mag[j] * sinf(cur_phs[j]);
        }
        
        /* Back to Time Domain */
        rfft((float*)cbuf, win_size/2, FFT_INVERSE);
        
        /* Assign to the output */
        for (j = 0; j < hop_size; j++) {
            outData[i+j] = pre_win[j + hop_size] + cur_win[j];
        }
        
        /* Move previous window */
        for (j = 0; j < win_size; j++) {
            pre_win[j] = (j < overlap_samples) ?
            pre_win[j + hop_size] : 0;
        }
        
        /* Update previous window */
        for (j = 0; j < win_size; j++) {
            pre_win[j] += cur_win[j];
        }
        
        for (j = 0; j < hop_size; j++)
        {
            buf[i+j] = buf[i+j] * (1.0 - parameters.mix) + outData[i+j] * parameters.mix;
        }
    }
}
