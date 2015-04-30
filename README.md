# FooHarmonizer
Pitch Shifter Plug-In

Author: Ryan Foo

Course: NYU's Fundamental of Digital Signal Theory Final Project

AU and VST formats supported. 

Follow the directions at http://www.insaneinthemainframe.at/2014/09/22/compile-the-juce-audio-plugin-demo-on-mac-osx-xcode-6-0-1/ to add the necessary AU and VST SDKs.



Resources: DAFX book Ch.7 for phase vocoding information


# Understanding the Plugin
1. Semitone Range is from 0.1 to 4.0. 1.0 is original pitch, 2.0 is one octave up, 0.5 is one octave down
2. Lowpass and Highpass filter options. You can choose the signal chain order of the two.
3. Uses STFT to pitch shift. 
