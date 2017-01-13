"""
MODULE:
    SignalMixer.py
DESCRIPTION:
    Module to mix real and complex signals with a sine wave, 
    cosine wave or complex exponential of a given frequency.
ALGORITHM:
    Code adapted from 
    http://www.haystack.mit.edu/cgi-bin/millstone_viewcvs.cgi/prototypes/SignalProcessing/c/mixer.c?rev=1.8&content-type=text/vnd.viewcvs-markup
EXCEPTIONS:
    FrequencyTooHigh
DEPENDENCIES:
    math, numpy
NOTE: 
    test() requires pylab and SignalGenerator.py
    
$Id: SignalMixer.py 2415 2009-12-21 14:57:45Z brideout $
"""

import math
import numpy

class FrequencyTooHigh(Exception):
    """
    DESCRIPTION:
        Exception class if the sampling frequency is too low compared to the signal frequency
    USAGE:
        raise FrequencyTooHigh(mixFreq, mixSampleFreq)
    PUBLIC ATTRIBUTES:
        self.freq (float)       = The frequency to be mixed with the signal in Hz
        self.sampleFreq (float) = The frequency at which to sample the signal to be mixed 
                                  (should be same as signal sample rate)
    AFFECTS:
        Creates an exception object
    EXCEPTIONS:
        None
    DEPENDENCIES:
        None
    WRITTEN BY:
        'Harendra Guturu':mailto:hguturu@haystack.mit.edu July 17, 2007
    """
    def __init__(self, freq, sampleFreq):
        """
        DESCRIPTION:
            Intializes the object variables of the class. 
        INPUTS:
            freq (float)              = The frequency to be mixed with the signal in Hz
            sampleFreq (float)        = The frequency at which to sample the signal to be mixed 
                                        (should be same as signal sample rate)
        RETURNS:
            Param1 (FrequencyTooHigh) = The pointer to the newly intialized object
        AFFECTS:
            Creates an exception object and intializes public attributes
        EXCEPTIONS:
            None
        DEPENDENCIES:
            None
        """
        self.freq = freq
        self.sampleFreq = sampleFreq

    def __str__(self):
        """
        DESCRIPTION:
            Helper function to print the exception.
        USAGE:
            None...automatically called when the object is printed.
        INPUTS:
            None
        RETURNS:
            Param1 (string) = The object as represented by a string.
        AFFECTS:
            None
        EXCEPTIONS:
            None
        DEPENDENCIES:
            None
        """
        return "Signal frequency too high compared to sample frequency, Rate: %s, Sample Rate: %s" \
                                                                        % (self.freq,self.sampleFreq)

def mix_complex_exponential(mixFreq, mixSampleFreq, oscilPhase, signal):
    """
    DESCRIPTION:
        Mixes a signal with a complex exponential of a given frequency 
    INPUTS:
        mixFreq (float)              = The frequency to be mixed with the signal in Hz
        mixSampleFreq (float)        = The frequency at which to sample the signal to be mixed in Hz 
                                       (should be same as signal sample rate)
        oscilPhase (float)           = Oscillator phase offset, radians
        signal (list, numpy.ndarray) = The array with all the signal values
    RETURNS:
        tuple of
            Param1 (numpy.ndarray)       = The array with the mixed signal values
            oscilPhase                   = phase of oscillator for next sample when called again
    AFFECTS:
        None
    EXCEPTIONS:
        FrequencyTooHigh             = If the mixFreq/mixSampleFreq >= pi/4.0
    DEPENDENCIES:
        numpy
    """
    mixedSignal = []
    if(abs(mixFreq/mixSampleFreq) >= numpy.pi/4.0):
        raise FrequencyTooHigh(mixFreq, mixSampleFreq)
    else:
        oscilCoef = 2.0*numpy.pi*mixFreq/mixSampleFreq
        t = numpy.arange(len(signal))
        expVal = numpy.cos(oscilCoef*t + oscilPhase)+numpy.sin(oscilCoef*t + oscilPhase)*1j
        signal = numpy.array(signal)
        newPhase = oscilCoef*len(signal) + oscilPhase
        while newPhase > numpy.pi/2.0:
            newPhase -= numpy.pi
        return((signal*expVal, newPhase))

def mix_sine(mixFreq, mixSampleFreq, signal):
    """
    DESCRIPTION:
        Mixes a signal with a sine wave of a given frequency
    INPUTS:
        mixFreq (float)              = The frequency to be mixed with the signal in Hz
        mixSampleFreq (float)        = The frequency at which to sample the signal to be mixed 
                                       (should be same as signal sample rate)
        signal (list, numpy.ndarray) = The array with all the signal values
    RETURNS:
        Param1 (numpy.ndarray)       = The array with the mixed signal values
    AFFECTS:
        None
    EXCEPTIONS:
        FrequencyTooHigh             = If the mixFreq/mixSampleFreq >= pi/4.0
    DEPENDENCIES:
        numpy
    """
    mixedSignal = []
    if(mixFreq/mixSampleFreq >= numpy.pi/4.0):
        raise FrequencyTooHigh(mixFreq, mixSampleFreq)
    else:
        oscilCoef = 2.0*numpy.pi*mixFreq/mixSampleFreq
        t = numpy.arange(len(signal))
        sinVal = numpy.sin(oscilCoef*t)
        signal = numpy.array(signal)
        return signal*sinVal

    
def mix_cosine(mixFreq, mixSampleFreq, signal):
    """
    DESCRIPTION:
        Mixes a signal with a cosine wave of a given frequency
    INPUTS:
        mixFreq (float)              = The frequency to be mixed with the signal in Hz
        mixSampleFreq (float)        = The frequency at which to sample the signal to be mixed 
                                       (should be same as signal sample rate)
        signal (list, numpy.ndarray) = The array with all the signal values
    RETURNS:
        Param1 (numpy.ndarray)       = The array with the mixed signal values
    AFFECTS:
        None
    EXCEPTIONS:
        FrequencyTooHigh             = If the mixFreq/mixSampleFreq >= pi/4.0
    DEPENDENCIES:
        numpy
    """
    mixedSignal = []
    if(mixFreq/mixSampleFreq >= numpy.pi/4.0):
        raise FrequencyTooHigh(mixFreq, mixSampleFreq)
    else:
        oscilCoef = 2.0*numpy.pi*mixFreq/mixSampleFreq
        t = numpy.arange(len(signal))
        cosVal = numpy.cos(oscilCoef*t)
        signal = numpy.array(signal)
        return signal*cosVal

def test():
    """
    DESCRIPTION:
        Tests the Module SignalMixer.py, also provides an example of how to use it
    USAGE:
        test()
    INPUTS:
        None
    RETURNS:
        None
    AFFECTS:
        Tests the various methods in this module
    EXCEPTIONS:
        None
    DEPENDENCIES:
        pylab, SignalGenerator.py
    """
    import pylab
    import SignalGenerator
    seed = -500
    toneFrequency = 5 # Hz
    toneStrength = 0 # dBm
    toneSNR = 40
    samplesPerPeriod = 10
    sampleRate = SignalGenerator.compute_sampling_rate(toneFrequency,samplesPerPeriod)
    timePeriod = 1 # second
    totalPeriods = SignalGenerator.compute_total_periods(toneFrequency,1)
    resistance = 50.0 # ohms

    signalGen = SignalGenerator.SignalGenerator(toneFrequency, toneStrength, resistance, 
                                                toneSNR, samplesPerPeriod, seed)
    
    del resistance, samplesPerPeriod, toneSNR, toneStrength, toneFrequency, seed
    
    totalSamples = int(signalGen.samplesPerPeriod*totalPeriods) # samples to get full periods
    
    rArray = []
    iArray = []
    cArray = []  
    for s in range(0,totalSamples):
        sigVal = signalGen.toneSignal.next()
        rArray.append(sigVal[0])
        iArray.append(sigVal[1])
        cArray.append(sigVal[0]+sigVal[1]*1j)
        del sigVal

    del totalSamples, totalPeriods

    print signalGen
    print "Examine plots to check if mixing changes the frequency spectrum appropriately."

    timeArray = numpy.arange(0,signalGen.numSamples)
    timeArray = timeArray*((1.0*signalGen.numPeriods)/(signalGen.numSamples*signalGen.signalFreq))
    
    (freqArray,fftArray) = SignalGenerator.compute_ctft(cArray,
                                                        signalGen.signalSampleRate, 
                                                        noShift=False)
    
    pylab.figure(1)
    pylab.plot(freqArray,numpy.abs(fftArray))
    pylab.axis([freqArray[0],freqArray[-1],-0.25,0.25])
    pylab.title('Original Spectra')
    pylab.xlabel('Frequency (Hz)')
    pylab.ylabel('Magnitude of Frequency Response [dBm] (Frequency domain)')
    pylab.figure(2)
    pylab.plot(freqArray,numpy.angle(fftArray))
    pylab.title('Original Spectra')
    pylab.xlabel('Frequency (Hz)')
    pylab.ylabel('Phase of Frequency Response [rads] (Frequency domain)')
        
    try:
        complexMixArray = mix_complex_exponential(-signalGen.signalFreq,signalGen.signalSampleRate,cArray)
    except FrequencyTooHigh, e:
        print e
    else:
        (freqCMixArray, fftCMixArray) = SignalGenerator.compute_ctft(complexMixArray,
                                                                     signalGen.signalSampleRate,
                                                                     noShift=False)
    pylab.figure(3)
    pylab.plot(freqArray,numpy.abs(fftCMixArray))
    pylab.axis([freqArray[0],freqArray[-1],-0.25,0.25])
    pylab.title('Spectra mixed with complex exponential')
    pylab.xlabel('Frequency (Hz)')
    pylab.ylabel('Magnitude of Frequency Response [dBm] (Frequency domain)')
    pylab.figure(4)
    pylab.plot(freqArray,numpy.angle(fftCMixArray))
    pylab.title('Spectra mixed with complex exponential')
    pylab.xlabel('Frequency (Hz)')
    pylab.ylabel('Phase of Frequency Response [rads] (Frequency domain)')
    
    try:
        sineMixArray = mix_sine(signalGen.signalFreq,signalGen.signalSampleRate,cArray)
    except FrequencyTooHigh, e:
        print e
    else:
        (freqSineArray, fftSineArray) = SignalGenerator.compute_ctft(sineMixArray,
                                                                     signalGen.signalSampleRate,
                                                                     noShift=False)
    pylab.figure(5)
    pylab.plot(freqArray,numpy.abs(fftSineArray))
    pylab.axis([freqArray[0],freqArray[-1],-0.25,0.25])
    pylab.title('Spectra mixed with sine wave')
    pylab.xlabel('Frequency (Hz)')
    pylab.ylabel('Magnitude of Frequency Response [dBm] (Frequency domain)')
    pylab.figure(6)
    pylab.plot(freqArray,numpy.angle(fftSineArray))
    pylab.title('Spectra mixed with sine wave')
    pylab.xlabel('Frequency (Hz)')
    pylab.ylabel('Phase of Frequency Response [rads] (Frequency domain)')
    
    try:
        cosMixArray = mix_cosine(signalGen.signalFreq,signalGen.signalSampleRate,cArray)
    except FrequencyTooHigh, e:
        print e
    else:
        (freqCosArray, fftCosArray) = SignalGenerator.compute_ctft(cosMixArray,
                                                                   signalGen.signalSampleRate,
                                                                   noShift=False)
    pylab.figure(7)
    pylab.plot(freqArray,numpy.abs(fftCosArray))
    pylab.axis([freqArray[0],freqArray[-1],-0.25,0.25])
    pylab.title('Spectra mixed with cosine wave')
    pylab.xlabel('Frequency (Hz)')
    pylab.ylabel('Magnitude of Frequency Response [dBm] (Frequency domain)')
    pylab.figure(8)
    pylab.plot(freqArray,numpy.angle(fftCosArray))
    pylab.title('Spectra mixed with cosine wave')
    pylab.xlabel('Frequency (Hz)')
    pylab.ylabel('Phase of Frequency Response [rads] (Frequency domain)')
    
    try:
        exceptArray = mix_cosine(signalGen.signalSampleRate,signalGen.signalSampleRate,cArray)
    except FrequencyTooHigh, e:
        print e
        print "Exception working."
    else:
        print "Exception not raised...may not be working."
    
    pylab.show()
    del freqCosArray, fftCosArray, freqSineArray, fftSineArray, freqCMixArray, fftCMixArray
    del freqArray, fftArray, timeArray
    del pylab, SignalGenerator
    del cArray, iArray, rArray
    del signalGen

if __name__ == "__main__":
    test()
