"""
MODULE:
    SignalGenerator.py
DESCRIPTION:
    Module to generate a tone signals
    The tone is generated in time domain in terms of voltages.
    This Module uses the tone frequency in hz, number of periods
    to look over and the number of samples to compute per period, the resistance of
    the generator, the strength of the signal and the signal to noise ratio and
    a seeded random generator for the noise
EXCEPTIONS:
    None
DEPENDENCIES:
    math, random, numpy
NOTE: 
    test() requires pylab
    
$Id: SignalGenerator.py 113 2007-07-25 21:01:11Z hguturu $
"""

import math
import random
import numpy

def compute_sampling_rate(signalFreq, samplesPerPeriod):
    """
    DESCRIPTION:
        Helper function to compute values needed to pass onto the SignalGenerator class.
        This function computes the sampling rate of a given frequency when the number of 
        samples per period is given. 
    INPUTS:
        signalFreq (float)       = The frequency of the signal in Hz 
        samplesPerPeriod (float) = The number of samples for a given period 
    RETURNS:
        Param1 (float)           = The sampling rate in Hz 
    AFFECTS:
        None
    EXCEPTIONS:
        None
    DEPENDENCIES:
        None
    """
    numPeriodsPerSecond = compute_total_periods(signalFreq, 1.0)
    return samplesPerPeriod*numPeriodsPerSecond*1.0


def compute_total_periods(signalFreq, samplingTimePeriod):
    """
    DESCRIPTION:
        Helper function to compute values need to pass only the SignalGenerator class.
        This function computes the number of periods of a given signal frequency when 
        a samplingTimePeriod (window) is defined.   
    INPUTS:
        signalFreq (float)         = The frequency of the signal in Hz 
        samplingTimePeriod (float) = The window over which the signal is being looked at in seconds 
    RETURNS:
        Param1 (float)             = The total number of periods of the signal being looked 
    AFFECTS:
        None
    EXCEPTIONS:
        None
    DEPENDENCIES:
        None
    """
    return samplingTimePeriod*signalFreq*1.0

def compute_samples_per_period(signalFreq, samplesPerSecond):
    """
    DESCRIPTION:
        Helper function to compute values need to pass only the SignalGenerator class.
        This function computes the number of samples taken in for each period of the signal
        given the signal frequency and the total samples per second  
    INPUTS:
        signalFreq (float)        = The frequency of the signal in Hz 
        samplingPerSecond (float) = The total number of samples in a second 
    RETURNS:
        Param1 (float)            = The numbers of samples for a single period of the signal 
    AFFECTS:
        None
    EXCEPTIONS:
        None
    DEPENDENCIES:
        None
    """
    return (samplesPerSecond*1.0)/signalFreq

def compute_power(strength):
    """
    DESCRIPTION: 
        Returns the power of the signal when a dBm value is passed
        The equation PdBm = 10 Log (P/1e-3) is solved for P when PdBm is given
    INPUTS: 
        strength (float) = Strength of the signal in dBm 
    RETURNS: 
        Param1 (float)   = Power in Watts. 
    AFFECTS:
        None
    EXCEPTIONS:
        None
    DEPENDENCIES:
        math.pow(...)
    """
    return math.pow(10.0,(strength/10.0))*1e-3 # solve the eqn: tonedBm = 10log(tonePower/1mW)

def compute_voltage(power,resistance):
    """
    DESCRIPTION:
        Returns the voltage by solving the equation P = V^2/R
    INPUTS:
        power (float)      = The power in Watts 
        resistance (float) = The resistance in Ohms 
    RETURNS:
        Param1 (float)     = The voltage in Volts 
    AFFECTS:
        None
    EXCEPTIONS:
        None
    DEPENDENCIES:
        math.sqrt(...)
    """
    return math.sqrt(power*resistance*1.0)

            
def compute_dtft(signalVals,noShift=True):
    """
    DESCRIPTION:
        compute_dtft calculate DTFT of a signal
    INPUTS:
        signalVals (list or numpy.ndarray) = finite-length input vector, whose length is L 
        noShift (true, false)              = if the value is True then the fft is evaluated from [0,2pi)
                                             else the fft is valued from [-pi,pi) 
    RETURNS:
        Param1 (numpy.ndarray)             = Vector of freqs where DTFT is computed in radians
        Param2 (numpy.ndarray)             = DTFT values as complex numbers    
    AFFECTS:
        None
    EXCEPTIONS:
        None
    DEPENDENCIES:
        numpy.array(...), numpy.fft.fft(...), numpy.fft.fftshift(...)
    ALGORITHM:
        http://www.eedsp.gatech.edu/Information/MATLAB_User_Guide/node96.html#figDTFTplotexample
    REFERENCE:
        http://en.wikipedia.org/wiki/Fourier_transform
    """
    L = len(signalVals)
    mid = int(math.ceil(L/2.0))
    transformVals = numpy.fft.fft( signalVals)
    if(noShift): # from [0,2pi)
        frequencies = (1.0/L) * numpy.array(range(0,L))
        #NOTE: alternate way using numpy.fft.fftfreq
        #frequencies = numpy.fft.fftfreq(L)
        #frequencies[mid:L+1] = frequencies[mid:L+1] + 1.0   # <--- move [pi,2pi) to [-pi,0)
    else: # from [-pi,pi)
        transformVals = numpy.fft.fftshift(transformVals)
        #frequencies = (1.0/L) * numpy.array(range(0,L))
        #frequencies[mid:L+1] = frequencies[mid:L+1] - 1.0   # <--- move [pi,2pi) to [-pi,0)
        #frequencies = numpy.fft.fftshift(frequencies)
        #NOTE: Alternate way using numpy.fft.fftfreq
        frequencies = numpy.fft.fftshift(numpy.fft.fftfreq(L))

    frequencies = frequencies*2.0*math.pi
    return (frequencies, transformVals)

def compute_ctft(signalVals, samplingRate, noShift=True):
    """
    DESCRIPTION:
        compute_ctft calculates CTFT of a signal
    INPUTS:
        signalVals (list or numpy.ndarray) = finite-length input vector, whose length is L
                                             in time domain
        samplingRate (float)               = the rate at which the continuous signal was sampled   
        noShift (true, false)              = if the value is True then the fft is evaluated from [0,2pi)
                                             else the fft is valued from [-pi,pi) 
    RETURN:
        Param1 (numpy.ndarray)             = Vector of freqs where CTFT is computed in frequency
        Param2 (numpy.ndarray)             = CTFT values as complex numbers
    AFFECTS:
        None
    EXCEPTIONS:
        None
    DEPENDENCIES:
        None
    ALGORITHM:
        Basically just converts from DTFT axis to CTFT axis
    REFERENCE:
        http://en.wikipedia.org/wiki/Fourier_transform
    """
    L = len(signalVals)
    (freqArray,fftArray) = compute_dtft(signalVals,noShift)
    
    #fftArray = fftArray/samplingRate # scaling by T since dtft has a 1/T scale, samplingRate = 1/T
    #numPeriods = L/(samplingRate/signalFreq)
    #fftArray = (fftArray/numPeriods)*signalFreq # to compensate for multiple periods and frequencies
    #NOTE: Above simplifies to fftArray = fftArray/L
    fftArray = fftArray/(1.0*L)
    freqArray = (freqArray*samplingRate)/(2.0*math.pi) # to get freq instead of angular freq
    
    return (freqArray,fftArray)

def compute_inv_dtft(transformVals, notShifted=True):
    """
    DESCRIPTION:
        compute_dtft calculate DTFT of a signal
    INPUTS:
        transformVals (list or numpy.ndarray) = finite-length input vector, whose length is L 
        notShifted (true, false)              = if the value is True then the fft was evaluated 
                                                from [0,2pi) else the fft was valued from [-pi,pi) 
    RETURNS:
        Param1 (numpy.ndarray)                = Vector of samples were DTFT is computed
        Param2 (numpy.ndarray)                = inverse DTFT values as complex numbers    
    AFFECTS:
        None
    EXCEPTIONS:
        None
    DEPENDENCIES:
        numpy.array(...), numpy.fft.fft(...), numpy.fft.fftshift(...)
    ALGORITHM:
        http://www.eedsp.gatech.edu/Information/MATLAB_User_Guide/node96.html#figDTFTplotexample
    REFERENCE:
        http://en.wikipedia.org/wiki/Fourier_transform
    NOTE:
        The right value of notShifted will ensure the spectrum is shifted properly before the ifft 
    """
    L = len(transformVals)
    mid = int(math.ceil(L/2.0))
    if(notShifted): # from [0,2pi)
        signalVals = numpy.fft.ifft(transformVals)
    else: # from [-pi,pi)
        transformVals = numpy.fft.ifftshift(transformVals)
        signalVals = numpy.fft.ifft(transformVals)

    samples = numpy.arange(0,L)
    return (samples, signalVals)


def compute_inv_ctft(transformVals, samplingRate, notShifted=True):
    """
    DESCRIPTION:
        compute_inv_ctft calculates inverse CTFT of a signal
    INPUTS:
        transformVals (list or numpy.ndarray) = finite-length input vector, whose length is L
                                                evaluation of [0,2pi) (in frequency domain)
        samplingRate (float)                  = the rate at which the continuous signal was sampled   
        signalFreq (float)                    = the frequency of the tone
        notShifted (true, false)              = if the value is True then the fft was evaluated 
                                                from [0,2pi) else the fft was valued from [-pi,pi) 
    RETURN:
        Param1 (numpy.ndarray)                = Vector of time where inverse CTFT is computed 
        Param2 (numpy.ndarray)                = inverse CFTT values as complex numbers
    AFFECTS:
        None
    EXCEPTIONS:
        None
    DEPENDENCIES:
        None
    ALGORITHM:
        Basically just converts the scaling due to sampling, also converts samples to time
    REFERENCE:
        http://en.wikipedia.org/wiki/Fourier_transform
    """
    L = len(transformVals)
    transformVals = transformVals*L
    (timeArray,signalArray) = compute_inv_dtft(transformVals,notShifted)
    timeArray = timeArray/(1.0*samplingRate)
   
    return (timeArray,signalArray)


class SignalGenerator:
    """
    DESCRIPTION:
        A container for a Tone signal generator with random white noise,
        The generated tone is in Time domain, the Frequency domain results 
        can be obtained by passing the samples to the CTFT function
    USAGE:
        toneGen = SignalGenerator(frequency, strength, snr, resitance, samplesperperiod, seed)
        sigVal1 = toneGen.toneSignal.next()
        sigVal2 = toneGen.toneSignal.next()
        ... so on 
    PUBLIC ATTRIBUTES:
        self.signalFreq (float)          = The frequency of the signal to be generated in Hz
        self.signalStrength (float)      = The strength of the signal in dBm
        self.signalSNR (float or None)   = The signal to noise ratio in dB
        self.signalPower (float)         = The power of the tone in Watts (computed from self.signalStrength)
        self.generatorResistance (float) = The resistances of the generator
        self.signalVoltage (float)       = The voltage of the signal in Volts 
                                           (computed using self.signalPower and self.signalResistance)
        self.randSeed (float or int)     = The seed to use for the pseudo random generator used for noise
        self.samplesPerPeriod (float)    = The number of samples to take in a given period of the tone
        self.signalSampleRate (float)    = The rate at which the signal is sampled 
                                           (computed from self.samplesPerPeriod and self.signalFreq)
        self.numSamples (int)            = The number of samples taken 
                                           (updated everytime self.signalSignal is accessed for a sample)
        self.numPeriods (float)          = The number of periods taken so far 
                                           (updated using self.numSamples and self.samplesPerPeriod after 
                                           every sampling is taken)
        self.toneSignal (generator)      = The tone generator object that will return the tone
        self.signalPhase (float)         = The phase of the signal in radians
    AFFECTS:
        Creates a generator object
    EXCEPTIONS:
        None
    DEPENDENCIES:
        random, math, numpy
    WRITTEN BY:
        'Harendra Guturu':mailto:hguturu@haystack.mit.edu July 17, 2007
    """
    def __init__(self, frequency, strength, resistance=50.0,
                 signalToNoise=None, samplesPerPeriod=1000.0, seed=None, phase=0.0):
        """
        DESCRIPTION:
            Intializes the object variables of the class.  
        INPUTS:
            signalFrequency (float)   = The frequency of the tone in Hz
            signalStrength (float)    = The strength of the tone in dBm    
            resitance (float)         = The resistance of the generator in Ohms (DEFAULT=50.0)    
            signalSNR (float or None) = The signal to noise ratio in dB 
                                        (DEFAULT=None, The signal will have no noise)    
            samplesPerPeriod (float)  = Number of samples to be taken per period (DEFAULT=1000.0)    
            seed (float or int)       = The seed for the noise 
                                        (DEFAULT=None, random.Random() default value will be 
                                        used, check Random documentation for version specific default, 
                                        may be system time or os random generator)
            phase (float)             = The phase of the signal (DEFAULT=0.0)
        RETURNS:
            Param1 (SignalGenerator)  = The pointer to the newly intialized object
        AFFECTS:
            Creates a generator object and intializes public attributes
        EXCEPTIONS:
            None
        DEPENDENCIES:
            None
        """
        self.signalFreq = frequency*1.0
        self.signalStrength = strength*1.0
        self.signalSNR = signalToNoise
        self.signalPower = compute_power(self.signalStrength)
        self.generatorResistance = resistance*1.0
        self.signalVoltage = compute_voltage(self.signalPower, self.generatorResistance)
        self.randNoiseSeed = seed
        self.samplesPerPeriod = samplesPerPeriod*1.0
        self.signalSampleRate = compute_sampling_rate(self.signalFreq, self.samplesPerPeriod)
        self.numSamples = 0
        self.numPeriods = 0.0
        self.signalPhase = phase*1.0
        self.toneSignal = self.toneGenerator(self.randNoiseSeed, self.signalFreq,
                                             self.signalVoltage, self.signalSampleRate,
                                             self.signalSNR, self.signalPhase)

    def toneGenerator(self,seed, toneFreq, toneVoltage=50e-3, 
                      samplingRate=100.0, signalToNoise=None,tonePhase=0.0):
        """
        DESCRIPTION:
            Generating a Tone of frequency fb
            x(t) = e^i*a*t = (cos(at) + isin(at)) (in time domain) (w/ noise if required)
        USAGE:
            a = toneGenerator(randSeed, toneFreq, toneVoltage, signalSampleRate, toneSNR)
            signalVal1 = a.next()
            signalVal2 = a.next()
            ... so on
        INPUTS:
            seed (float or int)           = The seed used for the pseudorandom number generator 
                                            for the noise simulation
            toneFreq (float)              = The frequency of the tone in hz    
            signalVoltage (float)         = The max voltage of the signal, DEFAULT = 50mV    
            samplingRate (float)          = The number of samples taken for a period, DEFAULT=100.0    
            signalToNoise (float or None) = the ratio between signal and noise in dB, 
                                            DEFAULT=None which results in no noise being added
            tonePhase (float)             = The phase of the tone in radians, DEFAULT=0
        RETURNS:
            Param1 (generator)            = The function it self is returned and use a.next() to get each 
                                            subsequenct signal value which is a tuple (real,imag).
        AFFECTS:
            Creates a generator object
        EXCEPTIONS:
            None
        DEPENDENCIES:
            math.sqrt(...), math.cos(...), math.sin(...), random.Random(...)
        """
        if(seed == None):
            rGen = random.Random() # seed is current time at run time
        else:
            rGen = random.Random(seed)

        samplingInterval = (2.0*math.pi)/samplingRate
        N = 0

        while(1):
            real = toneVoltage*math.cos(toneFreq*samplingInterval*N+tonePhase)
            imag = toneVoltage*math.sin(toneFreq*samplingInterval*N+tonePhase)
            
            toneMag = toneVoltage #math.sqrt(math.pow(real,2.0)+math.pow(imag,2.0)) #(alternate way to solve)

            if(signalToNoise != None):
                # noise with snr, extra term for signalMag to account for sampling  
                noiseMag = (toneMag*math.sqrt(samplingRate*1.0/toneFreq))/(math.pow(10.0,signalToNoise/20.0))       
                rVal = rGen.random() * 2.0 * math.pi
                real = real + noiseMag*math.cos(rVal)
                imag = imag + noiseMag*math.sin(rVal)
                
            self.numSamples = self.numSamples + 1
            self.numPeriods = self.numSamples*1.0/self.samplesPerPeriod
            yield (real,imag)
            N = N + 1
    
    def __repr__(self):
        """
        DESCRIPTION:
            Helper function to print all the attributes of the generator.    
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
        stringInfo = "Signal Frequency: %s [Hz] \n" % self.signalFreq
        stringInfo = stringInfo + "Signal Phase: %s [rads] \n" % self.signalPhase
        stringInfo = stringInfo + "Signal Strength: %s [dBm] \n" % self.signalStrength
        stringInfo = stringInfo + "Signal SNR: %s [dB] \n" % self.signalSNR
        stringInfo = stringInfo + "Signal Power: %s [Watts] \n" % self.signalPower
        stringInfo = stringInfo + "Generator Resistance: %s [Ohms] \n" % self.generatorResistance
        stringInfo = stringInfo + "Signal Voltage: %s [Volts] \n" % self.signalVoltage
        stringInfo = stringInfo + "Random Noise Seed: %s \n" % self.randNoiseSeed
        stringInfo = stringInfo + "Samples Per Period: %s [samples/period] \n" % self.samplesPerPeriod
        stringInfo = stringInfo + "Sample Rate: %s [samples/second] \n" % self.signalSampleRate
        stringInfo = stringInfo + "Number of Samples: %s [samples] \n" % self.numSamples
        stringInfo = stringInfo + "Number of Periods: %s [periods] " % self.numPeriods
        return stringInfo

def test():
    """
    DESCRIPTION:
        Tests the Module SignalGenerator.py, also provides an example of how to use it
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
        pylab
    """
    import pylab
    seed = -500
    toneFrequency = 10 # Hz
    tonePhase = math.pi*2.0 # rads
    toneStrength = 0 # dBm
    toneSNR = 40 # dB
    samplesPerPeriod = 15
    sampleRate = compute_sampling_rate(toneFrequency,samplesPerPeriod)
    if(sampleRate == 150.0):
        print "compute_sampling_rate(...) working."
    else:
        print "compute_sampling_rate(...) NOT working."
    if(samplesPerPeriod == compute_samples_per_period(toneFrequency, sampleRate)):
        print "compute_samples_per_period(...) working."
    else:
        print "compute_samples_per_period(...) NOT working."
    timePeriod = 1 # second
    totalPeriods = compute_total_periods(toneFrequency,1)
    if(totalPeriods == 10.0):
        print "compute_total_periods(...) working."
    else:
        print "compute_total_periods(...) NOT working."
    
    resistance = 50.0 # ohms

    gen = SignalGenerator(toneFrequency, toneStrength, resistance, toneSNR, 
                          samplesPerPeriod, seed, tonePhase)

    del resistance, samplesPerPeriod, toneSNR, toneStrength, toneFrequency, seed
    
    totalSamples = int(gen.samplesPerPeriod*totalPeriods) # samples to get full periods
    
    rArray = []
    iArray = []
    cArray = []

    for s in range(0,totalSamples):
        sigVal = gen.toneSignal.next()
        rArray.append(sigVal[0])
        iArray.append(sigVal[1])
        cArray.append(sigVal[0]+sigVal[1]*1j)
        del sigVal

    del totalSamples, totalPeriods

    print gen
    print "Examine plots to check if generator attributes match the plots."
    
    timeArray = numpy.arange(0,gen.numSamples)*((1.0*gen.numPeriods)/(gen.numSamples*gen.signalFreq))
    
    pylab.figure()
    pylab.plot(timeArray,rArray)
    pylab.plot(timeArray,iArray)
    pylab.title('Real (Blue) & Imaginary (Green) Voltage Component of tone')
    pylab.xlabel('Time (s)')
    pylab.ylabel('Signal [V] (Time domain)')


    (freqArray1,fftArray1) = compute_ctft(cArray,gen.signalSampleRate, noShift=True)

    pylab.figure()
    pylab.plot(freqArray1,numpy.log10(numpy.abs(fftArray1)/(math.sqrt(1.0e-3*gen.generatorResistance)))*20.0)
    pylab.title('Frequency spectra of tone with axis going from [0,2pi)')
    pylab.xlabel('Frequency (Hz)')
    pylab.ylabel('Magnitude of Frequency Response [dBm] (Frequency domain)')
    
#    pylab.figure()
#    pylab.plot(freqArray1,numpy.angle(fftArray1))
#    pylab.title('Frequency spectra of tone with axis going from [0,2pi)')
#    pylab.xlabel('Frequency (Hz)')
#    pylab.ylabel('Phase of Frequency Response [rads] (Frequency domain)')
    
    (timeArray1,cArray1) = compute_inv_ctft(fftArray1,gen.signalSampleRate, notShifted=True)
    pylab.figure()
    pylab.plot(timeArray1,numpy.real(cArray1))
    pylab.plot(timeArray1,numpy.imag(cArray1))
    pylab.title('Real (Blue) & Imaginary (Green) Voltage Component of tone')
    pylab.xlabel('Time (s)')
    pylab.ylabel('Signal [V] (Time domain)')
    
    (freqArray2,fftArray2) = compute_ctft(cArray,gen.signalSampleRate, noShift=False)

    pylab.figure()
    pylab.plot(freqArray2,numpy.log10(numpy.abs(fftArray2)/(math.sqrt(1.0e-3*gen.generatorResistance)))*20.0)
    pylab.title('Frequency spectra of tone with axis going from [-pi,pi)')
    pylab.xlabel('Frequency (Hz)')
    pylab.ylabel('Magnitude of Frequency Response [dBm] (Frequency domain)')
    
#    pylab.figure()
#    pylab.plot(freqArray2,numpy.angle(fftArray2))
#    pylab.title('Frequency spectra of tone with axis going from [0,2pi)')
#    pylab.xlabel('Frequency (Hz)')
#    pylab.ylabel('Phase of Frequency Response [rads] (Frequency domain)')

    (timeArray2,cArray2) = compute_inv_ctft(fftArray2,gen.signalSampleRate, notShifted=False)
    pylab.figure()
    pylab.plot(timeArray2,numpy.real(cArray2))
    pylab.plot(timeArray2,numpy.imag(cArray2))
    pylab.title('Real (Blue) & Imaginary (Green) Voltage Component of tone')
    pylab.xlabel('Time (s)')
    pylab.ylabel('Signal [V] (Time domain)')



    pylab.show()
    del timeArray2, cArray2, timeArray1, cArray1
    del freqArray2, fftArray2, freqArray1, fftArray1
    del pylab
    del timeArray, cArray, iArray, rArray
    del gen
 
if __name__ == "__main__":
    test()
   
