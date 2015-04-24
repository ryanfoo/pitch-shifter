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
    osamp = 32;
    frameSize   = WINDOW_SIZE/2;
    stepSize    = WINDOW_SIZE/osamp;
    freqPerBin  = currentSampleRate/WINDOW_SIZE;
    expct = 2. * M_PI * stepSize / WINDOW_SIZE;
    inFifoLatency = WINDOW_SIZE-stepSize;
    if (gRover == 0) gRover = inFifoLatency;
    
    memset(inFifoL,  0, WINDOW_SIZE*sizeof(float));
    memset(inFifoR,  0, WINDOW_SIZE*sizeof(float));
    memset(outFifoL, 0, WINDOW_SIZE*sizeof(float));
    memset(outFifoR, 0, WINDOW_SIZE*sizeof(float));
    memset(fftData, 0, WINDOW_SIZE*2*sizeof(float));
    memset(prevPhase, 0, (WINDOW_SIZE/2+1)*sizeof(float));
    memset(sumPhase, 0, (WINDOW_SIZE/2+1)*sizeof(float));
    memset(anaFreq, 0, WINDOW_SIZE*sizeof(float));
    memset(anaMagn, 0, WINDOW_SIZE*sizeof(float));
    memset(outData, 0, WINDOW_SIZE*2*sizeof(float));
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

    if (monoStatus == false)
    {
        initArrays();
        monoStatus = true;
    }
    
    for (int i = 0; i < numSamples; i++)
    {
        inFifoL[gRover] = samples[i*2];
        inFifoR[gRover] = samples[i*2+1];
        samples[i*2] = samples[i*2] * (1.0 - parameters.mix) + outFifoL[gRover-inFifoLatency] * parameters.mix;
        samples[i*2+1] = samples[i*2+1] * (1.0 - parameters.mix) + outFifoR[gRover-inFifoLatency] * parameters.mix;
        gRover++;
        
        if (gRover >= WINDOW_SIZE)
        {
            gRover = inFifoLatency;
            processSampleL();
            processSampleR();
        }
    }
}

// Process Stereo Data
void Shifter::processStereo(float* const left, float* const right, const int numSamples) noexcept
{
    jassert (left != nullptr && right != nullptr);
    
    if (stereoStatus == false)
    {
        initArrays();
        stereoStatus = true;
    }

    for (int i = 0; i < numSamples; i++)
    {
        inFifoL[gRover] = left[i];
        inFifoR[gRover] = right[i];
        left[i] = inFifoL[gRover] * (1.0 - parameters.mix) + outFifoL[gRover-inFifoLatency] * parameters.mix;
        right[i] = inFifoR[gRover] * (1.0 - parameters.mix) + outFifoR[gRover-inFifoLatency] * parameters.mix;
        gRover++;
        
        if (gRover >= WINDOW_SIZE)
        {
            gRover = inFifoLatency;
            processSampleL();
            processSampleR();
        }
    }
}

inline void Shifter::processSampleL()
{
    float tmp, phase;
    
    for (int i = 0; i < WINDOW_SIZE; i++)
    {
        window = -0.5f * cosf(2. * M_PI * (float)i/(float)WINDOW_SIZE) + 0.5f;
        fftData[i*2] = inFifoL[i] * window;
        fftData[i*2+1] = 0.;
    }
    
    stft(fftData, WINDOW_SIZE, -1);
    
    for (int i = 0; i <= frameSize; i++)
    {
        re = fftData[i*2];
        im = fftData[i*2+1];
        
        magn = 2. * sqrt(pow(re,2) + pow(im,2));
        phase = atan2f(im, re);
        
        tmp = phase - prevPhase[i];
        prevPhase[i] = phase;
        
        tmp -= (float)i*expct;
        
        qpd = tmp/M_PI;
        if (qpd >= 0) qpd += (qpd & 1);
        else qpd -= (qpd & 1);
        tmp -= M_PI*(float)qpd;
        
        tmp = osamp*tmp/(2*M_PI);
        tmp = (float)i*freqPerBin + tmp*freqPerBin;
        
        anaMagn[i] = magn;
        anaFreq[i] = tmp;
    }
    
    memset(synMagn, 0, WINDOW_SIZE*sizeof(float));
    memset(synFreq, 0, WINDOW_SIZE*sizeof(float));
    
    for (int i = 0; i <= frameSize; i++)
    {
        idx = i*parameters.pitch;
        if (idx <= frameSize)
        {
            synMagn[idx] += anaMagn[i];
            synFreq[idx] = anaFreq[i] * parameters.pitch;
        }
    }
    
    for (int i = 0; i <= frameSize; i++)
    {
        magn = synMagn[i];
        tmp = synFreq[i];
        
        tmp -= (float)i*freqPerBin;
        
        tmp /= freqPerBin;
        
        tmp = 2. * M_PI * tmp / osamp;
        
        tmp += (float)i*expct;
        
        sumPhase[i] += tmp;
        phase = sumPhase[i];
        
        fftData[i*2]    = magn * cosf(phase);
        fftData[i*2+1]  = magn * sinf(phase);
    }
    
    for (int i = WINDOW_SIZE+2; i < WINDOW_SIZE*2; i++) fftData[i] = 0.;
    
    stft(fftData, WINDOW_SIZE, 1);
    
    for (int i = 0; i < WINDOW_SIZE; i++)
    {
        window = -0.5f * cosf(2.*M_PI*(float)i/WINDOW_SIZE) + 0.5f;
        outData[i] += 2.*window*fftData[i*2]/(frameSize*osamp);
        if (isnan(outData[i])) outData[i] = 0;
    }
    
    for (int i = 0; i < stepSize; i++) outFifoL[i] = outData[i];
    
    memmove(outData, outData+stepSize, WINDOW_SIZE*sizeof(float));
    
    for (int i = 0; i < inFifoLatency; i++) inFifoL[i] = inFifoL[i+stepSize];
}

inline void Shifter::processSampleR()
{
    float tmp, phase;
    
    for (int i = 0; i < WINDOW_SIZE; i++)
    {
        window = -0.5f * cosf(2. * M_PI * (float)i/(float)WINDOW_SIZE) + 0.5f;
        fftData[i*2] = inFifoR[i] * window;
        fftData[i*2+1] = 0.;
    }
    
    stft(fftData, WINDOW_SIZE, -1);
    
    for (int i = 0; i <= frameSize; i++)
    {
        re = fftData[i*2];
        im = fftData[i*2+1];
        
        magn = 2. * sqrt(pow(re,2) + pow(im,2));
        phase = atan2f(im, re);
        
        tmp = phase - prevPhase[i];
        prevPhase[i] = phase;
        
        tmp -= (float)i*expct;
        
        qpd = tmp/M_PI;
        if (qpd >= 0) qpd += (qpd & 1);
        else qpd -= (qpd & 1);
        tmp -= M_PI*(float)qpd;
        
        tmp = osamp*tmp/(2*M_PI);
        tmp = (float)i*freqPerBin + tmp*freqPerBin;
        
        anaMagn[i] = magn;
        anaFreq[i] = tmp;
    }
    
    memset(synMagn, 0, WINDOW_SIZE*sizeof(float));
    memset(synFreq, 0, WINDOW_SIZE*sizeof(float));
    
    for (int i = 0; i <= frameSize; i++)
    {
        idx = i*parameters.pitch;
        if (idx <= frameSize)
        {
            synMagn[idx] += anaMagn[i];
            synFreq[idx] = anaFreq[i] * parameters.pitch;
        }
    }
    
    for (int i = 0; i <= frameSize; i++)
    {
        magn = synMagn[i];
        tmp = synFreq[i];
        
        tmp -= (float)i*freqPerBin;
        
        tmp /= freqPerBin;
        
        tmp = 2. * M_PI * tmp / osamp;
        
        tmp += (float)i*expct;
        
        sumPhase[i] += tmp;
        phase = sumPhase[i];
        
        fftData[i*2]    = magn * cosf(phase);
        fftData[i*2+1]  = magn * sinf(phase);
    }
    
    for (int i = WINDOW_SIZE+2; i < WINDOW_SIZE*2; i++) fftData[i] = 0.;
    
    stft(fftData, WINDOW_SIZE, 1);
    
    for (int i = 0; i < WINDOW_SIZE; i++)
    {
        window = -0.5f * cosf(2.*M_PI*(float)i/WINDOW_SIZE) + 0.5f;
        outData[i] += 2.*window*fftData[i*2]/(frameSize*osamp);
        if (isnan(outData[i])) outData[i] = 0;
    }
    
    for (int i = 0; i < stepSize; i++) outFifoR[i] = outData[i];
    
    memmove(outData, outData+stepSize, WINDOW_SIZE*sizeof(float));
    
    for (int i = 0; i < inFifoLatency; i++) inFifoR[i] = inFifoR[i+stepSize];
}
