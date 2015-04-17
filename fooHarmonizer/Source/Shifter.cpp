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

void Shifter::stft(float* buf, float frameSize, float sign)
{
    float wr, wi, arg, tmp, tr, ti, ur, ui;
    int i, j, bit, le, le2, k;
    
    for (i = 2; i < frameSize*2-2; i += 2)
    {
        for (bit = 2, j = 0; bit < frameSize*2; bit <<= 1)
        {
            if ((i & bit) != 0) j++;
            j <<= 1;
        }
        if (i < j)
        {
            tmp = buf[i];
            buf[i] = buf[j];
            buf[j] = tmp;
            tmp = buf[i+1];
            buf[i+1] = buf[j+1];
            buf[j+1] = tmp;
        }
    }
    
    float max = (float)(logl(frameSize)/logl(2.0) + 0.5);
    for (k = 0, le = 2; k < max; k++)
    {
        le <<= 1;
        le2 = le >> 1;
        ur = 1.0f;
        ui = 0.0f;
        arg = (float)M_PI / (le2 >> 1);
        wr = (float)cosf(arg);
        wi = (float)(sign * sinf(arg));
        for (j = 0; j < le2; j+= 2)
        {
            for (i = j; i < frameSize*2; i += le)
            {
                tr = buf[i+le2] * ur - buf[i+le2+1] * ui;
                ti = buf[i+le2] * ui + buf[i+le2+1] * ur;
                buf[i+le2] = buf[i] - tr;
                buf[i+le2+1] = buf[i+1] - ti;
                buf[i] += tr;
                buf[i+1] += ti;
            }
            tr = ur * wr - ui * wi;
            ui = ur * wi + ui * wr;
            ur = tr;
        }
    }
}

// Process Mono Data
void Shifter::processMono(float* const samples, const int numSamples) noexcept
{
    jassert (samples != nullptr);

    // Loop through sample buffer (STFT)
    for (int i = 0; i < numSamples; i++)
    {
        // Interpolate data
        if (i % 2 == 0) samples[i] = processSampleL(samples[i]);
        else samples[i] = processSampleR(samples[i]);
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
    float y, tmp, expct;
    long j;
    int i;

    inFIFO[gRover] = inSample;
    y = outFIFO[gRover - inFifoLatency];
    expct = 2.0f * M_PI * HOP_SIZE / WINDOW_SIZE;
    gRover++;
    
    if (gRover >= WINDOW_SIZE)
    {
        gRover = inFifoLatency;
        
        // WINDOWING
        for (i = 0; i < WINDOW_SIZE; i++)
        {
            win = -0.5 *cosf(2.0f * M_PI * (float)i / (float)WINDOW_SIZE) + 0.5;
            fftData[i*2] = (float)(inFIFO[i] * win);
            fftData[i*2+1] = 0.0f;
        }
        
        // ANALYSIS
        stft(fftData, WINDOW_SIZE, -1);
        for (i = 0; i <= WINDOW_SIZE/2; i++)
        {
            re = fftData[i*2];
            im = fftData[i*2+1];
            
            magn = 2.0f * sqrt(pow(re, 2) + pow(im, 2));
            phs = atan2f(im, re);
            
            tmp = phs - prev_phs[i];
            prev_phs[i] = phs;
            
            tmp -= (float)i * expct;
            
            int qpd = (float)(tmp/M_PI);
            if (qpd >= 0) qpd += qpd & 1;
            else qpd -= qpd & 1;
            tmp -= M_PI*(float)qpd;
            
            tmp = osamp * tmp / (2. * M_PI);
            
            tmp = (float)i * freqPerBin + tmp * freqPerBin;
            anaMagn[i] = (float)magn;
            anaFreq[i] = (float)tmp;
        }
        
        // PROCESSING
        memset(synMag, 0, WINDOW_SIZE*sizeof(float));
        memset(synFreq, 0, WINDOW_SIZE*sizeof(float));
        for (i = 0; i <= WINDOW_SIZE/2; i++)
        {
            j = (float)(i * parameters.pitch);
            if (j <= WINDOW_SIZE/2)
            {
                synMag[j] += anaMagn[i];
                synFreq[j] = anaFreq[i] * parameters.pitch;
            }
        }
        
        // SYNTHESIS
        for (i = 0; i <= WINDOW_SIZE/2; i++)
        {
            magn = synMag[i];
            tmp = synFreq[i];
            
            tmp -= (float)i * freqPerBin;
            tmp /= freqPerBin;
            tmp = 2.0f * M_PI * tmp / osamp;
            tmp += (float)i * expct;
            
            sumPhase[i] += (float)tmp;
            phs = sumPhase[i];
            
            fftData[i*2] = (float)(magn * cosf(phs));
            fftData[i*2+1] = (float(magn * sinf(phs)));
        }
        
        for (i = WINDOW_SIZE+2; i < WINDOW_SIZE*2; i++) fftData[i] = 0.0f;
        stft(fftData, WINDOW_SIZE, 1);
        
        for (i = 0; i < WINDOW_SIZE; i++)
        {
            win = -0.5 * cosf(2. * M_PI * (float)i / (float)WINDOW_SIZE) + 0.5;
            outData[i] += (float)(2. * win * fftData[i*2] / (WINDOW_SIZE*osamp));
        }
        
        for (i = 0; i < WINDOW_SIZE; i++)
        {
            outData[i] = outData[i+HOP_SIZE];
        }
        
        for (i = 0; i < inFifoLatency; i++) inFIFO[i] = inFIFO[i+HOP_SIZE];
    }
    
    return inSample*(1.0-parameters.mix) + y * parameters.mix;
}

inline float Shifter::processSampleR(float inSample)
{
    float y, tmp, expct;
    long j;
    int i;
    
    inFIFO[gRover] = inSample;
    y = outFIFO[gRover - inFifoLatency];
    expct = 2.0f * M_PI * HOP_SIZE / WINDOW_SIZE;
    gRover++;
    
    if (gRover >= WINDOW_SIZE)
    {
        gRover = inFifoLatency;
        
        // WINDOWING
        for (i = 0; i < WINDOW_SIZE; i++)
        {
            win = -0.5 *cosf(2.0f * M_PI * (float)i / (float)WINDOW_SIZE) + 0.5;
            fftData[i*2] = (float)(inFIFO[i] * win);
            fftData[i*2+1] = 0.0f;
        }
        
        // ANALYSIS
        stft(fftData, WINDOW_SIZE, -1);
        for (i = 0; i <= WINDOW_SIZE/2; i++)
        {
            re = fftData[i*2];
            im = fftData[i*2+1];
            
            magn = 2.0f * sqrt(pow(re, 2) + pow(im, 2));
            phs = atan2f(im, re);
            
            tmp = phs - prev_phs[i];
            prev_phs[i] = phs;
            
            tmp -= (float)i * expct;
            
            int qpd = (float)(tmp/M_PI);
            if (qpd >= 0) qpd += qpd & 1;
            else qpd -= qpd & 1;
            tmp -= M_PI*(float)qpd;
            
            tmp = osamp * tmp / (2. * M_PI);
            
            tmp = (float)i * freqPerBin + tmp * freqPerBin;
            anaMagn[i] = (float)magn;
            anaFreq[i] = (float)tmp;
        }
        
        // PROCESSING
        memset(synMag, 0, WINDOW_SIZE*sizeof(float));
        memset(synFreq, 0, WINDOW_SIZE*sizeof(float));
        for (i = 0; i <= WINDOW_SIZE/2; i++)
        {
            j = (float)(i * parameters.pitch);
            if (j <= WINDOW_SIZE/2)
            {
                synMag[j] += anaMagn[i];
                synFreq[j] = anaFreq[i] * parameters.pitch;
            }
        }
        
        // SYNTHESIS
        for (i = 0; i <= WINDOW_SIZE/2; i++)
        {
            magn = synMag[i];
            tmp = synFreq[i];
            
            tmp -= (float)i * freqPerBin;
            tmp /= freqPerBin;
            tmp = 2.0f * M_PI * tmp / osamp;
            tmp += (float)i * expct;
            
            sumPhase[i] += (float)tmp;
            phs = sumPhase[i];
            
            fftData[i*2] = (float)(magn * cosf(phs));
            fftData[i*2+1] = (float(magn * sinf(phs)));
        }
        
        for (i = WINDOW_SIZE+2; i < WINDOW_SIZE*2; i++) fftData[i] = 0.0f;
        stft(fftData, WINDOW_SIZE, 1);
        
        for (i = 0; i < WINDOW_SIZE; i++)
        {
            win = -0.5 * cosf(2. * M_PI * (float)i / (float)WINDOW_SIZE) + 0.5;
            outData[i] += (float)(2. * win * fftData[i*2] / (WINDOW_SIZE*osamp));
        }
        
        for (i = 0; i < WINDOW_SIZE; i++)
        {
            outData[i] = outData[i+HOP_SIZE];
        }
        
        for (i = 0; i < inFifoLatency; i++) inFIFO[i] = inFIFO[i+HOP_SIZE];
    }

    
    return inSample*(1.0-parameters.mix) + y * parameters.mix;
}
