/*
  ==============================================================================

    This file was auto-generated by the Introjucer!

    It contains the basic framework code for a JUCE plugin processor.

  ==============================================================================
*/

#include "PluginProcessor.h"
#include "PluginEditor.h"


//==============================================================================
FooHarmonizerAudioProcessor::FooHarmonizerAudioProcessor()
{
    shifter.monoData = new Shifter::data();
    shifter.leftData = new Shifter::data();
    shifter.rightData = new Shifter::data();
    
    shifter.initArrays(shifter.monoData);
    shifter.initArrays(shifter.leftData);
    shifter.initArrays(shifter.rightData);
}

FooHarmonizerAudioProcessor::~FooHarmonizerAudioProcessor()
{
}

//==============================================================================
const String FooHarmonizerAudioProcessor::getName() const
{
    return JucePlugin_Name;
}

int FooHarmonizerAudioProcessor::getNumParameters()
{
    return 0;
}

float FooHarmonizerAudioProcessor::getParameter (int index)
{
    return 0.0f;
}

void FooHarmonizerAudioProcessor::setParameter (int index, float newValue)
{
}

const String FooHarmonizerAudioProcessor::getParameterName (int index)
{
    return String();
}

const String FooHarmonizerAudioProcessor::getParameterText (int index)
{
    return String();
}

const String FooHarmonizerAudioProcessor::getInputChannelName (int channelIndex) const
{
    return String (channelIndex + 1);
}

const String FooHarmonizerAudioProcessor::getOutputChannelName (int channelIndex) const
{
    return String (channelIndex + 1);
}

bool FooHarmonizerAudioProcessor::isInputChannelStereoPair (int index) const
{
    return true;
}

bool FooHarmonizerAudioProcessor::isOutputChannelStereoPair (int index) const
{
    return true;
}

bool FooHarmonizerAudioProcessor::acceptsMidi() const
{
   #if JucePlugin_WantsMidiInput
    return true;
   #else
    return false;
   #endif
}

bool FooHarmonizerAudioProcessor::producesMidi() const
{
   #if JucePlugin_ProducesMidiOutput
    return true;
   #else
    return false;
   #endif
}

bool FooHarmonizerAudioProcessor::silenceInProducesSilenceOut() const
{
    return false;
}

double FooHarmonizerAudioProcessor::getTailLengthSeconds() const
{
    return 0.0;
}

int FooHarmonizerAudioProcessor::getNumPrograms()
{
    return 1;   // NB: some hosts don't cope very well if you tell them there are 0 programs,
                // so this should be at least 1, even if you're not really implementing programs.
}

int FooHarmonizerAudioProcessor::getCurrentProgram()
{
    return 0;
}

void FooHarmonizerAudioProcessor::setCurrentProgram (int index)
{
}

const String FooHarmonizerAudioProcessor::getProgramName (int index)
{
    return String();
}

void FooHarmonizerAudioProcessor::changeProgramName (int index, const String& newName)
{
}

// Updates the Pitch Shifter's Parameters every time a knob is turned
void FooHarmonizerAudioProcessor::updateShifter(void)
{
    // Get Shifter Parameters
    Shifter::Parameters shifterParams = shifter.getParameters();
    
    // Set Shifter Parameters
    shifterParams.pitch = pitchVal;
    shifterParams.mix = mixVal;
    shifterParams.lpf = lpVal;
    shifterParams.hpf = hpVal;
    shifterParams.order = order;
    shifterParams.filter = filter;
    
    // Actually set them in parameters struct
    shifter.setParameters(shifterParams);
    
    // If filter is on, update the filter values
    if (shifterParams.filter)
    {
        shifter.updateLPFilter();
        shifter.updateHPFilter();
    }
}

//==============================================================================
void FooHarmonizerAudioProcessor::prepareToPlay (double sampleRate, int samplesPerBlock)
{
    // Use this method as the place to do any pre-playback
    // initialisation that you need..
}

void FooHarmonizerAudioProcessor::releaseResources()
{
    // When playback stops, you can use this as an opportunity to free up any
    // spare memory, etc.
}

void FooHarmonizerAudioProcessor::processBlock (AudioSampleBuffer& buffer, MidiBuffer& midiMessages)
{
    // In case we have more outputs than inputs, this code clears any output
    // channels that didn't contain input data, (because these aren't
    // guaranteed to be empty - they may contain garbage).
    // I've added this to avoid people getting screaming feedback
    // when they first compile the plugin, but obviously you don't need to
    // this code if your algorithm already fills all the output channels.
    for (int i = getNumInputChannels(); i < getNumOutputChannels(); ++i)
        buffer.clear (i, 0, buffer.getNumSamples());
    
    // Process Audio
    if (getNumInputChannels() == 1)
    {
        // Obtain channel 1 data
        float *monoChannel = buffer.getWritePointer(0);
        
        // Pitch Shifting processing
//        shifter.processChannel(monoChannel, buffer.getNumSamples());
        shifter.processMono(monoChannel, buffer.getNumSamples());
    }
    else if (getNumInputChannels() == 2)
    {
        // Obtain channel 1 and 2 data
        float *leftChannel = buffer.getWritePointer(0), *rightChannel = buffer.getWritePointer(1);
        
        // Pitch Shifting processing
        shifter.processStereo(leftChannel, rightChannel, buffer.getNumSamples());
//        shifter.processMono(leftChannel, buffer.getNumSamples());
//        shifter.processMono(rightChannel, buffer.getNumSamples());
    }
}

//==============================================================================
bool FooHarmonizerAudioProcessor::hasEditor() const
{
    return true; // (change this to false if you choose to not supply an editor)
}

AudioProcessorEditor* FooHarmonizerAudioProcessor::createEditor()
{
    return new FooHarmonizerAudioProcessorEditor (*this);
}

//==============================================================================
void FooHarmonizerAudioProcessor::getStateInformation (MemoryBlock& destData)
{
    // You should use this method to store your parameters in the memory block.
    // You could do that either as raw data, or use the XML or ValueTree classes
    // as intermediaries to make it easy to save and load complex data.
}

void FooHarmonizerAudioProcessor::setStateInformation (const void* data, int sizeInBytes)
{
    // You should use this method to restore your parameters from this memory block,
    // whose contents will have been created by the getStateInformation() call.
}

//==============================================================================
// This creates new instances of the plugin..
AudioProcessor* JUCE_CALLTYPE createPluginFilter()
{
    return new FooHarmonizerAudioProcessor();
}
