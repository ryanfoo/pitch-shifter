/*
 ==============================================================================
 
 This file was auto-generated by the Introjucer!
 
 It contains the basic framework code for a JUCE plugin editor.
 
 ==============================================================================
 */

#include "PluginProcessor.h"
#include "PluginEditor.h"

// Creates a slider with customized params
void FooHarmonizerAudioProcessorEditor::createSlider(Slider &slider, Slider::SliderStyle style, double defaultVal,
                                                     double min, double max, double incr, std::string name)
{
    // Define Slider parameters
    slider.setSliderStyle(style);
    slider.setRange(min,max,incr);
    slider.setTextBoxStyle(Slider::NoTextBox, false, 90, 0);
    slider.setPopupDisplayEnabled(true, this);
    slider.setTextValueSuffix(" " + name);
    slider.setValue(defaultVal);
    
    // Add a listener (user interaction)
    slider.addListener(this);
    // Add slider to editor
    addAndMakeVisible(&slider);
    // set name
    slider.setComponentID(name);
}

// Creates a label with customized params
void FooHarmonizerAudioProcessorEditor::createLabel(Label &label, std::string name)
{
    label.setSize(50,20);
    label.setEnabled(true);
    label.setText(name, dontSendNotification);
    addAndMakeVisible(label);
    label.setVisible(true);
    label.isAlwaysOnTop();
}

//==============================================================================
FooHarmonizerAudioProcessorEditor::FooHarmonizerAudioProcessorEditor (FooHarmonizerAudioProcessor& p)
: AudioProcessorEditor (&p), processor (p)
{
    // Make sure that before the constructor has finished, you've set the
    // editor's size to whatever you need it to be.
    setSize (500, 500);
    
    // Create Sliders
    createSlider(mixSlider, Slider::Rotary, processor.mixVal, 0.0, 1.0f, 0.01f, "Mix");
    createSlider(pitchSlider, Slider::Rotary, processor.pitchVal, -12, 12, 1, "Semitone");
    createSlider(lowpassSlider, Slider::Rotary, processor.lpVal, 0.0, 20000, 10, "Lowpass Filter");
    createSlider(highpassSlider, Slider::Rotary, processor.hpVal, 0.0, 20000, 10, "Highpass Filter");
    // Create Labels
    createLabel(mixText, "Mix");
    createLabel(pitchText, "Semitone");
    createLabel(lowpassText, "LP Filter");
    createLabel(highpassText, "HP Filter");
}

FooHarmonizerAudioProcessorEditor::~FooHarmonizerAudioProcessorEditor()
{
}

//==============================================================================
void FooHarmonizerAudioProcessorEditor::paint (Graphics& g)
{
    // Background
    g.fillAll (Colours::white);
    
    // Texts and Line Graph
    g.setColour (Colours::black);
    g.setFont (40.0f);
    g.drawFittedText ("Foo Harmonizer", getLocalBounds(), Justification::centredTop, 1);
    g.setFont(25.0f);
    g.drawFittedText("by Ryan Foo", getLocalBounds(), Justification::centredBottom, 2);
    
    // Backgrounds
    g.setColour(Colours::grey);
}

void FooHarmonizerAudioProcessorEditor::resized()
{
    // This is generally where you'll want to lay out the positions of any
    // subcomponents in your editor..
    
    // Set bounds for sliders
    mixSlider.setBounds(100, 125, 50, 50);
    pitchSlider.setBounds(350, 125, 50, 50);
    lowpassSlider.setBounds(100, 375, 50, 50);
    highpassSlider.setBounds(350, 375, 50, 50);
    
    // Set bounds for texts
    mixText.setBounds(100, 180, 10, 10);
    pitchText.setBounds(350, 180, 10, 10);
    lowpassText.setBounds(100, 430, 50, 10);
    highpassText.setBounds(350, 430, 50, 10);
}

// Callback from listener
void FooHarmonizerAudioProcessorEditor::sliderValueChanged(Slider* slider)
{
    if (slider->getComponentID().compare("Mix") == 0) {
        processor.mixVal = slider->getValue();
    }
    else if (slider->getComponentID().compare("Semitone") == 0) {
        processor.pitchVal = slider->getValue(); 
    }
    else if (slider->getComponentID().compare("Lowpass Filter") == 0) {
        processor.lpVal = slider->getValue();
    }
    else if (slider->getComponentID().compare("Highpass Filter") == 0) {
        processor.hpVal = slider->getValue();
    }
    
    processor.updateShifter();
}