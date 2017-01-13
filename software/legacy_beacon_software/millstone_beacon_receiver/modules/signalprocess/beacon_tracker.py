"""
MODULE:
    SignalTracker.py
DESCRIPTION:
    Module to find the frequency of a given signal and to return a copy of the signal
    after frequency shifting it to baseband. An ACF is used to find the frequency.
ALGORITHM:
    ACF code found in 
    http://www.haystack.mit.edu/cgi-bin/millstone_viewcvs.cgi/prototypes/SignalProcessing/python/lpm.py?rev=1.19&content-type=text/vnd.viewcvs-markup
EXCEPTIONS:
    None
DEPENDENCIES:
    math, numpy, warnings, lpmpkg.Lpm
NOTE: 
    test() requires time, pylab, SignalGenerator.py
    
$Id: SignalTracker.py 2416 2009-12-21 15:00:43Z brideout $
"""

import math
import numpy
import matplotlib.mlab
import matplotlib.pylab
import scipy.optimize

import sys

import warnings
# sys.path.append('../lpmpkg')
import Lpm as lpm
warnings.filterwarnings(action='ignore',category=DeprecationWarning,module='Lpm')

def __weighted_leastsq_min_fn__(slope, lag, phase, weight):
    return sum(weight*(phase - slope*lag)**2)

def __detect_shift_from_baseband_acfmethod__(samplingRate,dataArray,integerLags=5):
    """
    DESCRIPTION
        Determines how far from baseband the average weight of the frequency spectra is shifted
        Function coded with the idea of using ACF to determine the shift
    INPUTS:
        samplingRate (float)            = The sampling rate of the data in Hz
        dataArray (list, numpy.ndarray) = The voltage data array
        integerLags (int)               = The number of integers to lag by when computing
                                          the LagProfileMatrix, DEFAULT=5
    RETURNS:
        Param1 (float)                  =  The frequency that it is shifted from baseband
                                           in Hz. It will be from [-BW/2,BW/2)
    AFFECTS:
        None
    EXCEPTIONS:
        None
    DEPENDENCIES:
        numpy.linspace(...)
    Note:
        * The returned frequency will the aliased frequency if the signal is sampled below nyquist.
          __convert_unaliased_freq__ should change it to unaliased value
        * The acf inverse trapezoidal (this method) is apparently better for 
          beacon signal acquisitionn, while acf trapezoidal is better for radar signals 
    """
    nZeros = len(dataArray)
    fractionality = 1
    if integerLags > nZeros:
        # limit lag calculations
        intLags = nZeros - 1
    else:
        intLags = integerLags

    lpm  = lpm.LagProfileMatrix(nZeros, fractionality, intLags)
    lpm.AcfMultiplyAccumulate(dataArray)
    
    maxGate = lpm.getMaximumGatingAcfMatrixInverseTrapezoidal()
    acfMat = lpm.AcfMatrixInverseTrapezoidal(maxGate)
    
    phase = map(__find_neg_phase__, numpy.real(acfMat)[0], numpy.imag(acfMat)[0])
    phase = __fix_wrapping__(phase)
    
    lag = numpy.linspace(0,intLags-1, intLags)

    """
    # This code does equal-weight linear fit 
    #[slope, intercept] = numpy.polyfit(lag, phase, 1) # linear fit
    [slope, intercept] = matplotlib.mlab.polyfit(lag, phase, 1) # linear fit
    """
    
    # Weighted least squares using real part as weight and forced zero intercept
    initialguess = phase[-1]/lag[-1]
    (slope, fopt, iter, funcalls, warnflag) = scipy.optimize.fmin(__weighted_leastsq_min_fn__, initialguess, args=(lag, phase, numpy.real(acfMat)[0]), disp=0, full_output=True)
    if warnflag == 2 or abs(slope) > 1e20:
        # try a plain old linear (equal weight) leastsquares fit
        equalweight = numpy.ones(len(phase), numpy.float)
        (slope, fopt, iter, funcalls, warnflag) = scipy.optimize.fmin(__weighted_leastsq_min_fn__, initialguess, args=(lag, phase, equalweight), disp=0, full_output=True)
    if warnflag == 2 or abs(slope) > 1e20:
        matplotlib.pylab.figure()
        matplotlib.pylab.subplot(311)
        matplotlib.pylab.plot(lag, phase)
        matplotlib.pylab.subplot(312)
        matplotlib.pylab.plot(lag, acfMat.real[0])
        matplotlib.pylab.subplot(313)
        matplotlib.pylab.plot(lag, acfMat.imag[0])
        matplotlib.pylab.xlabel('Lag increment = %.0f usec' % (1.0e-6 / samplingRate))
        
        matplotlib.pylab.show()
        
        raise ValueError, 'Weighted least squares fit failed'

    """
    matplotlib.pylab.figure()
    matplotlib.pylab.subplot(311)
    matplotlib.pylab.plot(lag, phase)
    matplotlib.pylab.subplot(312)
    matplotlib.pylab.plot(lag, acfMat.real[0])
    matplotlib.pylab.subplot(313)
    matplotlib.pylab.plot(lag, acfMat.imag[0])
    matplotlib.pylab.xlabel('Lag increment = %.0f usec' % (1.0e-6 / samplingRate))
    
    matplotlib.pylab.show()
    """
    
    return (slope*samplingRate)/(2.0*math.pi)


def __detect_shift_from_baseband_fftmethod__(samplingRate,dataArray,nFFT=None):
    """
    DESCRIPTION
        Determines how far from baseband the average weight of the frequency spectra is shifted,
        using a spectral method.
        
        Not implemented yet.

    INPUTS:
 
        samplingRate (float)            = The sampling rate of the data in Hz
        dataArray (list, numpy.ndarray) = The voltage data array
        nFFT (int)                      = Number of frequency bins to use when computing
                                          FFT for frequency detection.  Default = 
                                          length of dataArray.
    RETURNS:
        Param1 (float)                  =  The frequency that it is shifted from baseband
                                           in Hz. It will be from [-BW/2,BW/2)
    AFFECTS:
        None
    EXCEPTIONS:
        None
    DEPENDENCIES:
        numpy.fft.fft(...)
    Note:
        * The returned frequency will be the aliased frequency if the signal is sampled below nyquist.
          __convert_unaliased_freq__ should change it to unaliased value
    """
    
    raise ValueError, 'Method not yet implemented'
    
def __fix_wrapping__(phaseArray):
    """
    DESCRIPTION:
        Fixes the wrapping of phase by keeping track of each previous wrap
        EX: phaseArray = [0, pi/2, -pi/2] = [0, pi/2, 3pi/2] (it will keep it going in one dir)
    INPUTS:
        phaseArray (list, numpy.ndarray) = The array with phases between -pi to pi
    RETURNS:
        Param1 (numpy.ndarray)           =  An array with phases wrapped accordingly with 
                                            no restriction -pi to pi.
                                            The wrapping will be calculated based on previous term.
    AFFECTS:
        None
    EXCEPTIONS:
        None
    DEPENDENCIES:
        numpy.array(...), math.pi
    """
    length = len(phaseArray)
    if(length<3):
        return phaseArray
    else:
        origArray = numpy.array(phaseArray)
        delta = 0
        for i in range(2,len(phaseArray)):
            if(abs(origArray[i]-origArray[i-1]) > math.pi): # increments the delta everytime there is a wrap
                if(origArray[1] < 0):
                    delta -= 2.0*math.pi
                else:
                    delta += 2.0*math.pi
            phaseArray[i] += delta
        return phaseArray

def __find_neg_phase__(re, im):
    """
    DESCRIPTION:
        Returns the negative phase of a complex value
        The phase DOES take into account the sign of re, im. It is restricted to -pi to pi
    INPUTS:
        re (float)     = The real component
        im (float)     = The imaginary component
    RETURNS:
        Param1 (float) = -math.atan2(im/re) in radians
    AFFECTS:
        None
    EXCEPTIONS:
        None
    DEPENDENCIES:
        math.atan2(...)
    """
    phase = math.atan2(im,re)
    return -phase

def convert_to_unaliased_freq(sateFreq, samplingRate, aliasedFreq):
    """
    DESCRIPTION:
        Converts the aliasedFreq to the unaliased frequency.
        It checks for spectral inversion by making sure the number of aliases are even or odd.
    INPUTS:
        sateFreq (float)      = The central frequency of the signal ie. the beacon frequencies in Hz
                                without being doppler shifted...this value should be constant
                                for a whole pass and most likely be the CERTO beacon value
        samplingRate (float)  = The sampling rate of the signal in Hz
        aliasedFreq (float)   = The frequency that is detected due to aliasing from [-BW/2,BW/2)
    RETURNS:
        Param1 (float)        = The unaliased frequency in Hz
    AFFECTS:
        None
    EXCEPTIONS:
        None
    DEPENDENCIES:
        numpy.floor(...)
    """
    fullAliases = numpy.floor((sateFreq*1.0)/samplingRate)
    partAlias = (sateFreq*1.0)/samplingRate - fullAliases

    modFreq = sateFreq % samplingRate
    if((samplingRate/2.0 - modFreq >= 0) and aliasedFreq >= 0): # center and alias in [0,-BW/2)
        #print 1
        return (sateFreq - modFreq) + aliasedFreq
    elif((samplingRate/2.0 - modFreq >=0) and aliasedFreq < 0): # center in [0,-BW/2) and alias in [-BW/2,0)
        if((modFreq-aliasedFreq)*2.0 >= samplingRate): # the alias has wrapped, >= instead of > since
                                                       # right boundary is missing [0,BW/2)
            #print '2a'
            return (sateFreq - modFreq) + aliasedFreq + samplingRate
        else: # has not wrapped
            #print '2b'
            return (sateFreq - modFreq) + aliasedFreq
    elif((samplingRate/2.0 - modFreq < 0) and aliasedFreq >= 0): # center in [-BW/2,0) and alias in [0,BW/2)
        if((modFreq-aliasedFreq)*2.0 < samplingRate): # the alias has wrapped, < instead of <= since
                                                       # right boundary is missing [0,BW/2)
            #print '3a'
            return (sateFreq - modFreq) + aliasedFreq
        else: # has not wrapped
            #print '3b'
            return (sateFreq - modFreq) + aliasedFreq + samplingRate
    else: # center and alias in [-BW/2,0)
        #print 4
        return (sateFreq - modFreq) + aliasedFreq + samplingRate
    
def convert_to_aliased_freq(unaliasedFreq, samplingRate):
    """
    DESCRIPTION:
        Converts the sateFreq to the aliased frequency.
        It checks for spectral inversion by making sure the number of aliases are even or odd.
    INPUTS:
        unaliasedFreq (float) = The true frequency of the signal 
        samplingRate (float)  = The sampling rate of the signal in Hz
        shifted (True, False) = If shifted=True then the Param1 is in range [-BW/2,BW/2)
                                if it is False, then it is in the range [0,-BW)
                                DEFAULT=True
    RETURNS:
        Param1 (float)        = The aliased frequency in Hz in the range of [-BW/2,BW/2) or [0,-BW) 
                                depending on shifted, DEFAULT=[-BW/2,BW/2)
    AFFECTS:
        None
    EXCEPTIONS:
        None
    DEPENDENCIES:
        numpy.floor(...)
    """
    fullAliases = numpy.floor((unaliasedFreq*1.0)/samplingRate)
    partAlias = (unaliasedFreq*1.0)/samplingRate - fullAliases
    if(partAlias<0.5): # it is the [0,BW/2) region when partAlias = [0,0.5)
        aliasedFreq = samplingRate*partAlias
    else: # it has wrapped to the [-BW/2,0) region when partAlias = [0.5,1)
        aliasedFreq = samplingRate*partAlias - samplingRate
    return aliasedFreq

def convert_to_complex(vals):
    """
    DESCRIPTION:
        Takes in an array of tuples. The first component is real and
        the second is imaginary and returns a numpy array of complex values
    INPUTS:
        vals (list, numpy.ndarray) = Array of tuples Ex: [(0,0)...(2,45)]
                                     which represent the real and imaginary components
    RETURNS:
        Param1 (numpy.ndarray)     = Array of complex values Ex: [0+0j,...,2+45j]
    AFFECTS:
        None
    EXCEPTIONS:
        None
    DEPENDENCIES:
        numpy.array(...)
    """
    return numpy.array(map(lambda (re,im): re+im*1j, vals))

def track_frequency_shift(samplingRate, dataArray, algorithm='acfmethod', trackParam=None):
    """
    DESCRIPTION:
        Detects the frequency by which the average weight of the frequency spectra 
        is of from baseband and shifts the spectra to baseband.
    INPUTS:
        samplingRate (float)            = The sampling rate of the data in Hz
        dataArray (list, numpy.ndarray) = The voltage data array
        algorithm (string)              = Frequency tracking algorithm to use:
                                          'acfmethod' - use ACF based method
                                          'fftmethod' - use FFT based method
                                          Default = 'acfmethod'
        trackParam (int)                = Control parameter for frequency tracking calculation:
                                          Algorithm = 'acfmethod':  number of integer lags
                                          for the ACF to use (default = 5)
                                          Algorithm = 'fftmethod':  number of FFT bins to use
                                          (default = None - uses number of points in dataArray)
    RETURNS:
        Param1 (float)                  = The aliased frequency at which the data was found
                                          in Hz within [-BW/2,BW/2)
        Param2 (float)                  = The unaliased frequency at which the data 
                                          was recorded in Hz
    AFFECTS:
        None
    EXCEPTIONS:
        None
    DEPENDENCIES:
        None
    """
    
    # PJE 2009-01-30: some of the difficulty here is that the number of lags in the ACF
    # method are fixed.  This means that for very low Doppler shifts, we don't get far enough
    # out on the ACF phase slope and get an inaccurate / wrong final frequency measurement.
    #
    # What we should go to is a minimum Lag count, with the processor free to increase or
    # decrease the lags until we hit the real part going to some level such as 0.1 or something.
    #
    # This will work well for slow Dopplers.  For large Dopplers, the minimum lag will ensure
    # that we get a reasonable coarse frequency resolution (although we have to deal with phase
    # wraps - but I think the code to do this is already in place here).
    #
    # Implementing this should allow the tracking loop to use iterations to zero in on the 
    # final frequency.  We should monitor the number of lags chosen as a function of iteration
    # and ensure that it goes up as the Doppler shift goes down.
    #
    # This explains all the weird tracking shift estimates seen when things are close to zero
    # Doppler.
    
    
    
    if algorithm not in ('acfmethod', 'fftmethod'):
        raise ValueError, 'algorithm must be either acfmethod or fftmethod'
    if algorithm == 'acfmethod':
        # ACF tracking
        if trackParam == None:
            trackParam = 5
        aliasedFreq = __detect_shift_from_baseband_acfmethod__(samplingRate, dataArray,trackParam)
    elif algorithm == 'fftmethod':
        aliasedFreq = __detect_shift_from_baseband_fftmethod__(samplingRate, dataArray,trackParam)

    return aliasedFreq


def unwrap_phase(complexArr):
    """unwrap_phase takes an input complex numpy array, and then returns a double array of 
        unwrapped phases in radians.  That is, the returned array can have values above 2pi or below -2pi 
        depending on which way the phase progresses.  This method assumes the phase changes 
        slowly so that unwrapping is possible.  An error is raised if the phase change
        from one measurement to another is greater than 90 degrees.
        
    Inputs: complexArr - an array of complex values
    
    Returns: an array of floats, of length = len(complexArr).  The first value is 0.0 
        cand each following value is the value before plus the phase difference
        between the present and the preceeding.  An error is raised if two value differ by more 
        than pi/2
    """
    phaseArr = numpy.zeros(complexArr.shape[0], numpy.float64)
    
    phaseArr[0] = numpy.angle(phaseArr[0])
    
    for i in range(1, complexArr.shape[0]):
        diff = numpy.angle(complexArr[i]) - numpy.angle(complexArr[i-1])
        # force diff to be between -pi and pi
        while diff < -1.0*numpy.pi:
            diff += 2.0*numpy.pi
        while diff > 1.0*numpy.pi:
            diff -= 2.0*numpy.pi
        # verify not more than pi/2
        if abs(diff) > numpy.pi / 2.0:
            raise ValueError, 'phase change greater than 90 degrees between %f (pos %i) and %f (pos %i)' % \
                            (numpy.angle(complexArr[i-1]), i-1, numpy.angle(complexArr[i]), i)
                            
        phaseArr[i] = phaseArr[i-1] + diff
        
    return(phaseArr)


def test():
    """
    DESCRIPTION:
        Tests the Module SignalTracker.py, also provides an example of how to use it
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
        time, pylab, SignalGenerator.py
    """
    import time
    import pylab
    import SignalGenerator
    toneFrequency = 150.012e6
    toneStrength = -90.0
    toneSNR = 1
    samplesPerSecond = 80e3 # 100 kHz = sampling rate
    samplingTimePeriod = 0.01 # 1 milli second
    resistance = 50.0
    seed = -500
    totalWindows = int(15.0/samplingTimePeriod) # 30 seconds/samplingTimePeriod
    print totalWindows
    toneT0 = 150.024e6
    toneTN = 150.000e6
    inputData = [] # will have the input frequencies, format: (window#, inputfreq, numsamples)
    outputData = [] # will have the output (tracked) frequencies, format: (window# outputfreq, numsamples)
    totData  = []
    for i in range(0,totalWindows):
        toneFreq = (toneT0 - toneTN)/(1.0+math.exp(0.3*(i*samplingTimePeriod-(totalWindows*samplingTimePeriod-1)/2.0))) + toneTN
        totalPeriods = SignalGenerator.compute_total_periods(toneFreq,samplingTimePeriod)
        samplesPerPeriod = SignalGenerator.compute_samples_per_period(toneFreq, samplesPerSecond)
        obj = SignalGenerator.SignalGenerator(toneFreq,toneStrength,resistance,toneSNR,samplesPerPeriod,seed)
        totalSamples = int(obj.samplesPerPeriod*totalPeriods) # for given window
        inputData.append((i , toneFreq, totalSamples))
        
        s = 0
        sigList = []
        curr = time.time()
        while(s < totalSamples):
            sigVal = obj.toneSignal.next()
            sigList.append(sigVal)
            s = s + 1
        cArray = convert_to_complex(sigList)
        totData.extend(cArray)
        outFreq = track_frequency_shift(samplesPerSecond, cArray)
        unAliasedFreq = convert_to_unaliased_freq(toneFrequency, samplesPerSecond, outFreq)
        outputData.append((i, unAliasedFreq, len(cArray)))
        if(i % 100 == 0):
            print "Time for last %s windows till %s was %s." % (100, i, time.time()-curr)
            curr = time.time()

    inputData = numpy.array(inputData)
    outputData = numpy.array(outputData)
    pylab.figure(1)
    pylab.plot(inputData[:,0],inputData[:,1])
    pylab.title('Input Frequency Curve')
    pylab.xlabel('Window Number')
    pylab.ylabel('Frequency [Hz]')
    
    pylab.figure(2)
    pylab.plot(outputData[:,0],outputData[:,1])
    pylab.title('Output Frequency Curve')
    pylab.xlabel('Window Number')
    pylab.ylabel('Frequency [Hz]')
    
    pylab.figure(3)
    pylab.plot(inputData[:,0],inputData[:,1], label='Input')
    pylab.plot(outputData[:,0],outputData[:,1], label='Output')
    pylab.title('Input/Output Frequency Curve Overlay')
    pylab.xlabel('Window Number')
    pylab.ylabel('Frequency [Hz]')
    pylab.legend()
    
    pylab.figure(4)
    pylab.plot(inputData[:,0],inputData[:,1]-outputData[:,1])
    pylab.title('Difference between Input and Output')
    pylab.xlabel('Window Number')
    pylab.ylabel('Error (input-output)')
    
    pylab.figure(5)
    pylab.specgram(totData,Fs=samplesPerSecond)
    pylab.colorbar()
    pylab.title('Specgram of Input Tone')
    pylab.xlabel('Time [s]')
    pylab.ylabel('Frequency Hz')
        
    pylab.show()
      

if __name__ == "__main__":
    test()
