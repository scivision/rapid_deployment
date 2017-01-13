"""
MODULE:
    SignalFilter.py
DESCRIPTION:
    Module to low-pass and band-pass filter a given signal. 
    Module will also allow dropping and increasing the bandwidth of the signal.
EXCEPTIONS:
    ValueError  (python default exception)
DEPENDENCIES:
    math, numpy, scipy.signal
NOTE: 
    test() requires pylab and SignalGenerator.py

$Id: SignalFilter.py 2438 2010-01-15 21:26:05Z pje $
"""

import math
import numpy
import scipy.signal
import MillstoneScipy.filter

class SignalFilter:
    
    def __init__(self):
        self.filterObj = None

    def __lfilter_zi__(self, b,a):
        """
        DESCRIPTION:
            Compute the zi state from the filter parameters. see [Gust96].
        INPUTS:
            b (numpy.ndarray)      = The numerator coefficients of the filter
            a (numpy.ndarray)      = The denominator coefficients of the filter
        RETURNS:
            Param1 (numpy.ndarray) = The zi state for the filter
        AFFECTS:
            None
        EXCEPTIONS:
            None
        DEPENDENCIES:
            numpy
        ALGORITHM:
            Based on:
            [Gust96] Fredrik Gustafsson, Determining the initial states in forward-backward 
            filtering, IEEE Transactions on Signal Processing, pp. 988--992, April 1996, 
            Volume 44, Issue 4
        NOTE:
            * Code found at http://www.scipy.org/Cookbook/FiltFilt
            * NOT BEING USED...
        """
        n = max(len(a),len(b))
        
        zin = (numpy.eye(n-1) - numpy.hstack((-a[1:n,numpy.newaxis],
                                              numpy.vstack((numpy.eye(n-2),numpy.zeros(n-2)))))) 
        
        zid = b[1:n] - a[1:n]*b[0]
        
        zi_matrix=numpy.linalg.inv(zin)*(numpy.matrix(zid).transpose())
        zi_return=[]
    
        #convert the result into a regular array (not a matrix)
        for i in range(len(zi_matrix)):
            zi_return.append(float(zi_matrix[i][0]))
        
        return numpy.array(zi_return)
    
    
    def __filtfilt__(self, b,a,x):
        """
        DESCRIPTION:
            Implements a zero phase delay filter that processes the signal in the forward 
            and backward direction removing the phase delay. 
            The order of the filter is the double of the original filter order. 
            The function also computes the initial filter parameters in order to 
            provide a more stable response (via lfilter_zi).
        INPUTS:
            b (list, numpy.ndarray) = The numerator coefficients of the filter
            a (list, numpy.ndarray) = The denominator coefficients of the filter
            x (list, numpy.ndarray) = The data array to be filtered in time domain
        RETURNS:
            Param1 (numpy.ndarray)  = The filtered data array in time domain
        AFFECTS:
            None
        EXCEPTIONS:
            ValueError
        DEPENDENCIES:
            numpy, scipy.signal
        ALGORITHM:
            Based on:
            [Gust96] Fredrik Gustafsson, Determining the initial states in forward-backward 
            filtering, IEEE Transactions on Signal Processing, pp. 988--992, April 1996, 
            Volume 44, Issue 4
        NOTE:
            * Code found at http://www.scipy.org/Cookbook/FiltFilt
            * NOT BEING USED...
        """
        x = numpy.array(x)
        a = numpy.array(a)
        b = numpy.array(b)
        
        #For now only accepting 1d arrays
        ntaps = max(len(a),len(b))
        edge = ntaps*3
            
        if x.ndim != 1:
            raise ValueError, "Filiflit is only accepting 1 dimension arrays."
    
        #x must be bigger than edge
        if x.size < edge:
            raise ValueError, "Input vector needs to be bigger than 3 * max(len(a),len(b)."
            
        
        if len(b)<len(a):
            b=numpy.r_[b,numpy.zeros(len(a)-len(b))]
        elif len(b)>len(a):
            a=numpy.r_[a,numpy.zeros(len(b)-len(a))]
             
        
        zi = __lfilter_zi__(b,a)
        
        #Grow the signal to have edges for stabilizing 
        #the filter with inverted replicas of the signal
        s = numpy.r_[2*x[0]-x[edge:1:-1],x,2*x[-1]-x[-1:-edge:-1]]
        #in the case of one go we only need one of the extrems 
        # both are needed for filtfilt
        
        (y,zf) = scipy.signal.lfilter(b,a,s,-1,zi*s[0])
    
        (y,zf) = scipy.signal.lfilter(b,a,numpy.flipud(y),-1,zi*y[-1])
        
        return numpy.flipud(y[edge-1:-edge+1])
    
    
        
    
    def __insert_zero_samples__(self, upFactor, dataArray):
        """
        DESCRIPTION:
            Inserts zeroes between samples of a signal by a given factor
        INPUTS:
            upfactor (int)                  = factor by which to add zeroes
            dataArray (list, numpy.ndarray) = The array with the samples values
        RETURNS:
            Param1 (numpy.ndarray)          = The array with sampled values and added zeroes
        AFFECTS:
            None
        EXCEPTIONS:
            ValueError                      = When upFactor is not an int
        DEPENDENCIES:
            numpy
        NOTE:
            * If the downFactor is <=1 then the original dataArray is returned.
            * NOT BEING USED...scipy.signal.resample(...) being used for upsampling
        """
        if(type(1) == type(upFactor) or upFactor == int(upFactor)):
            upFactor = upFactor - 1
            if(upFactor <= 0):
                return numpy.array(dataArray)
            else:
                upArray = []
                for datum in dataArray:
                    j = upFactor
                    upArray.append(datum)
                    while (j>0):
                        upArray.append(0.0+0.0*1j)
                        j -= 1
                
                return numpy.array(upArray)
        else:
            raise ValueError, "The insert zero sample factor [%s] has to be an integer." % (upFactor) 
    
    def __decimate__(self, downFactor,dataArray):
        """
        DESCRIPTION:
            Drops samples by the given factor
        INPUTS:
            downfactor (int)                = factor by which to drop samples
            dataArray (list, numpy.ndarray) = The array with the samples values
        RETURNS:
            Param1 (numpy.ndarray)          = The array with dropped sample values
        AFFECTS:
            None
        DEPENDENCIES:
            numpy
        EXCEPTIONS:
            ValueError                      = When downFactor is not an int
        NOTE:
            * If the downFactor is <=1 then the original dataArray is returned.
            * NOT BEING USED...scipy.signal.resample(...) being used for downsampling
        """
        if(type(1) == type(downFactor) or downFactor == int(downFactor)):
            if(downFactor <= 1):
                return numpy.array(dataArray)
            else:
                return dataArray[::downFactor]
        else:
            raise ValueError, "The decimation factor [%s] has to be an integer." % (downFactor)
    
    
    def create_low_pass_filter(self, order, cutOffFreq,sampleRate, window='hamming'):
        """
        DESCRIPTION:
            Designs a low pass filter and returns the frequency response and the
            corresponding frequencies
        USAGE:
            (freqs,filter) = create_low_pass_filter(nTaps,somefreq, arrayofvalues, arrayoffreqs)
        INPUTS:
            order (int)            = The number of delays/taps used in the filter
            cutOffFreq (float)     = The frequency at which to cut the filter
                                     NOTE: Does not check if the freq is greater than half nyquist
            sampleRate (float)     = The frequency at which the signal was sampled
            window (string, tuple) = The type of window to use (refer to scipy.signal.window 
                                     for possible windows and format)
        RETURNS:
            Param1 (numpy.ndarray) = The coefficients of the fir filter
        AFFECTS:
            None
        EXCEPTIONS:
            None
        DEPENDENCIES:
            scipy.signal
        NOTE:
            * NOT BEING USED...
        """
        cutOff = cutOffFreq/sampleRate # normalize to be between 0 and 1, 1 corresponds to pi
        hcoef = scipy.signal.firwin(order,cutOff,window=window)
        return hcoef
    
    def apply_fir_filter(self, filterCoefs, dataArray):
        """
        DESCRIPTION:
            Applies the filter to the time domain array when the filter coefficients are passed
            The filter is a zero phase delay filter (Herbert FIR need phase delay?)
        INPUTS:
            filterCoefs (list, numpy.ndarray) = The coefficients of the filter as returned by 
                                                     create_low_pass_filter(...)
            dataArray (list, numpy.ndarray)   = The dataArray in time domain to filter
        RETURNS:
            Param1 (numpy.ndarray)            = The dataArray in time domain after filtering
            Param2 (numpy.ndarray)            = The final conditions of the filter
                                                Only returned if the intialConditions are not None
        AFFECTS:
            None
        EXCEPTIONS:
            None
        DEPENDENCIES:
            None
        Note:
            * NOT BEING USED...
        """
        #zf=__lfilter_zi__(numpy.array(filterCoefs),numpy.r_[1,numpy.zeros(len(filterCoefs)-1)])
        #(tmpArray,zf) = scipy.signal.lfilter(filterCoefs,1,dataArray,zi=zf)
        #(dataArray,zf) = scipy.signal.lfilter(filterCoefs,1,dataArray,zi=zf)
        dataArray = __filtfilt__(filterCoefs,numpy.array([1]),dataArray)
        return dataArray
    
    def up_convert(self, upFactor, dataArray):
        """
        DESCRIPTION:
            Increases the sampling rate of a given signal
        INPUTS:
            upFactor (int)                  = The factor by which to up convert the signal
            dataArray (list, numpy.ndarray) = The time domain signal 
        RETURNS:
            Param1 (numpy.ndarray)          = The up converted dataArray
        AFFECTS:
            None
        EXCEPTIONS:
            ValueError                      = When upFactor is not an int
        DEPENDENCIES:
            numpy, scipy.signal
        NOTE:
            * If the upFactor is <=1 then the original dataArray is returned.
            * Upsampling does not properly interpolate SNR since spectrum is treated
              as bandlimited when the noise is not bandlimited.
        """
        if(type(1) == type(upFactor) or upFactor == int(upFactor)):
            if(upFactor>1):
                upDataArray = scipy.signal.resample(dataArray,len(dataArray)*upFactor)
                return numpy.array(upDataArray)
            else:
                return numpy.array(dataArray)
        else:
            raise ValueError, "The up sample factor [%s] has to be an integer." % (upFactor) 
    
    def down_convert_deprecated(self, downFactor, dataArray):
        """
        DESCRIPTION:
            Decreases the sampling rate of a given signal
        INPUTS:
            downFactor (int)                = The factor by which to down convert the signal
            dataArray (list, numpy.ndarray) = The time domain signal 
            order (int)                     = The order of the anti-alias filter to be used
        RETURNS:
            Param1 (numpy.ndarray)          = The down converted dataArray
        AFFECTS:
            None
        EXCEPTIONS:
            ValueError                      = When downFactor is not an int
        DEPENDENCIES:
            numpy, scipy.signal
        NOTE:
            * If the downFactor is <=1 then the original dataArray is returned.
        """
        if(type(1) == type(downFactor) or downFactor == int(downFactor)):
            if(downFactor>1):
                downDataArray = scipy.signal.resample(dataArray,len(dataArray)/downFactor)
                return numpy.array(downDataArray)
            else:
                return numpy.array(dataArray)
        else:
            raise ValueError, "The down sample factor [%s] has to be an integer." % (downFactor) 
    
    def down_convert(self, downFactor, dataArray):
        """ A better version of down_convert which handles state if called repeatedly, as when the
        original sample stream is broken up into windows.
        
        Shamelessly stolen from MATLAB's resample function within the signal processing toolbox.
        
        Implements things the "slow" way - a faster way would be to use an equivalent of MATLAB's
        upfirdn function but this doesn't exist in official Python scipy yet.
        """
    
        if self.filterObj == None:
        
            # design parameter for Kaiser window LPF
            beta = 5
            # design parameter for LPF filter length
            N = 10
            
            # construct FIR low pass filter
            fc = 1.0/2/downFactor
            L = 2*N*downFactor
            self.filterCoef = scipy.signal.firwin(L-1, 2*fc, beta)
            # DEBUG - use an allpass filter
            #self.filterCoef[:] = 0.0
            #self.filterCoef[len(self.filterCoef)/2] = 1.0
            
            # initialize tracking filter object
            self.filterObj = MillstoneScipy.filter.lfilter(self.filterCoef)
        
        # apply anti-alias filter - note that this might return a zero length vector the first time it is called
        lpfData = self.filterObj.get_previous_lfilter(dataArray)

        # return decimated data (safely anti-alias-filtered) - this might be zero length
        return lpfData[::downFactor]
    
    
    def change_bandwidth(self, upFactor, downFactor, dataArray):
        """
        DESCRIPTION:
            Upsamples/Decimates the dataArray so that it is at a new bandwidth, but combines
            the interpolation and antialiasing filters (if both are necessary) so that 
            less noise is introduced . 
            Essentially calling up_convert(...) and down_convert(...) directly
        INPUTS:
            upFactor (int)                  = The factor by which to up convert the signal
            downFactor (int)                = The factor by which to down convert the signal
            dataArray (list, numpy.ndarray) = The data array
        RETURNS:
            Param1 (numpy.ndarray)          = The new bandwidth dataArray
        AFFECTS:
            None
        EXCEPTIONS:
            ValueError                      =  If upFactor and downFactor are not ints
        DEPENDENCIES:
            numpy, scipy.signal
        """
        if((type(1) == type(upFactor) or upFactor == int(upFactor))  # type(1)=<type 'int'>
           and (type(1) == type(downFactor) or downFactor == int(downFactor))):
            if((upFactor <= 1 and downFactor <=1) or upFactor == downFactor):
                return numpy.array(dataArray)
            elif(upFactor > 1 and downFactor <=1):
                return up_convert(upFactor, dataArray)
            elif(upFactor <=1 and downFactor >1):
                return down_convert(downFactor, dataArray)
            else:
                newDataArray = up_convert(upFactor, dataArray)
                newDataArray = down_convert(downFactor, newDataArray)
                return newDataArray
        else:
            raise ValueError, "The up sample factor [%s] and down sample factor [%s] have to be integers."\
                                                                     % (upFactor, downFactor) 
    
    def test(self):
        """
        DESCRIPTION:
            Tests the Module SignalFilter.py, also provides an example of how to use it
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
        toneFrequency = 10 #Hz
        toneStrength = 0 # dBm
        toneSNR = 30 # dB
        sampleRate = 200 # nyquist
        samplesPerPeriod = SignalGenerator.compute_samples_per_period(toneFrequency, sampleRate) 
        timePeriod = 1 # second
        totalPeriods = SignalGenerator.compute_total_periods(toneFrequency,timePeriod)
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
        print "Examine plots to check if the frequency spectrum of filters is correct."
    
        # Start Regular Signal
        timeArray = numpy.arange(0,signalGen.numSamples)
        timeArray = timeArray*((1.0*signalGen.numPeriods)/(signalGen.numSamples*signalGen.signalFreq))
        
        (freqArray,fftArray) = SignalGenerator.compute_ctft(cArray,
                                                            signalGen.signalSampleRate,
                                                            noShift=False)
        pylab.figure()
        pylab.plot(freqArray,numpy.log10(numpy.abs(fftArray)))
        pylab.title('Original Spectra')
        pylab.xlabel('Frequency (Hz)')
        pylab.ylabel('Magnitude of Frequency Response (log scale) (Frequency domain)')
    #    pylab.figure()
    #    pylab.plot(freqArray,numpy.angle(fftArray))
    #    pylab.title('Original Spectra')
    #    pylab.xlabel('Frequency (Hz)')
    #    pylab.ylabel('Phase of Frequency Response (Frequency domain)')
        
        pylab.figure()
        pylab.plot(timeArray,rArray)
        pylab.plot(timeArray,iArray)
        pylab.title('Real (Blue) & Imaginary (Green) Voltage Component of Original')
        pylab.xlabel('Time (s)')
        pylab.ylabel('Signal [V] (Time domain)')
        # End Regular Signal
        
        # Start Up Convert
        upfactor = 6.0
        upArray = up_convert(upfactor, cArray)
        (upfreqArray,upfftArray) = SignalGenerator.compute_ctft(upArray,
                                                                signalGen.signalSampleRate*upfactor,
                                                                noShift=False)
        pylab.figure()
        pylab.plot(upfreqArray,numpy.log10(numpy.abs(upfftArray)))
        pylab.title('Upsampled Spectra with Interpolation Filter')
        pylab.xlabel('Frequency (Hz)')
        pylab.ylabel('Magnitude of Frequency Response (log scale) (Frequency domain)')
    #    pylab.figure()
    #    pylab.plot(upfreqArray,numpy.angle(upfftArray))
    #    pylab.title('Upsampled Spectra with Interpolation Filter')
    #    pylab.xlabel('Frequency (Hz)')
    #    pylab.ylabel('Phase of Frequency Response (Frequency domain)')
    
                                                                    
        (timeArray1,cArray1) = SignalGenerator.compute_inv_ctft(upfftArray,
                                                                signalGen.signalSampleRate*upfactor,
                                                                notShifted=False)                                                      
        pylab.figure()
        pylab.plot(timeArray1,numpy.real(cArray1))
        pylab.plot(timeArray1,numpy.imag(cArray1))
        pylab.title('Real (Blue) & Imaginary (Green) Voltage Component of Original')
        pylab.xlabel('Time (s)')
        pylab.ylabel('Signal [V] (Time domain)')
        # End Up Convert
    
        # Start Down Convert
        downfactor = 5.0
        downArray = down_convert(downfactor, cArray)
        (downfreqArray,downfftArray) = SignalGenerator.compute_ctft(downArray,
                                                                signalGen.signalSampleRate*1.0/downfactor,
                                                                noShift=False)
        pylab.figure()
        pylab.plot(downfreqArray,numpy.log10(numpy.abs(downfftArray)))
        pylab.title('Downsampled Spectra with Anti-Aliasing Filter')
        pylab.xlabel('Frequency (Hz)')
        pylab.ylabel('Magnitude of Frequency Response (log scale) (Frequency domain)')
    #    pylab.figure()
    #    pylab.plot(downfreqArray,numpy.angle(downfftArray))
    #    pylab.title('Downsampled Spectra with Anti-Aliasing Filter')
    #    pylab.xlabel('Frequency (Hz)')
    #    pylab.ylabel('Phase of Frequency Response (Frequency domain)')
        
        (timeArray2,cArray2) = SignalGenerator.compute_inv_ctft(downfftArray,
                                                                signalGen.signalSampleRate*1.0/downfactor,
                                                                notShifted=False)                                                      
        pylab.figure()
        pylab.plot(timeArray2,numpy.real(cArray2))
        pylab.plot(timeArray2,numpy.imag(cArray2))
        pylab.title('Real (Blue) & Imaginary (Green) Voltage Component of Original')
        pylab.xlabel('Time (s)')
        pylab.ylabel('Signal [V] (Time domain)')
        # End Down Convert
        
        
        # Start Change
        upfactor = 7.0
        downfactor = 10.0
        newArray = change_bandwidth(upfactor,downfactor, cArray)
        (newfreqArray,newfftArray) = SignalGenerator.compute_ctft(newArray,
                                                                signalGen.signalSampleRate*(1.0*upfactor)/downfactor,
                                                                noShift=False)
        pylab.figure()
        pylab.plot(newfreqArray,numpy.log10(numpy.abs(newfftArray)))
        pylab.title('Resampled Spectra')
        pylab.xlabel('Frequency (Hz)')
        pylab.ylabel('Magnitude of Frequency Response (log scale) (Frequency domain)')
    #    pylab.figure()
    #    pylab.plot(newfreqArray,numpy.angle(newfftArray))
    #    pylab.title('Resampled Spectra')
    #    pylab.xlabel('Frequency (Hz)')
    #    pylab.ylabel('Phase of Frequency Response (Frequency domain)')
    
        (timeArray3,cArray3) = SignalGenerator.compute_inv_ctft(newfftArray,
                                                                signalGen.signalSampleRate*(1.0*upfactor)/downfactor,
                                                                notShifted=False)                                                      
        pylab.figure()
        pylab.plot(timeArray3,numpy.real(cArray3))
        pylab.plot(timeArray3,numpy.imag(cArray3))
        pylab.title('Real (Blue) & Imaginary (Green) Voltage Component of Original')
        pylab.xlabel('Time (s)')
        pylab.ylabel('Signal [V] (Time domain)')
        # End Change
        
        
        # other stuff
        seed = -500
        toneFrequency = 150.012e6 #Hz
        toneStrength = 0 # dBm
        toneSNR = 30 # dB
        sampleRate = 100e3 # nyquist
        timePeriod = 1 # second
        totalPeriods = SignalGenerator.compute_total_periods(toneFrequency,timePeriod)
        resistance = 50.0 # ohms
        signalGen1 = SignalGenerator.SignalGenerator(toneFrequency, toneStrength, resistance, 
                                                    toneSNR, 
                                                    SignalGenerator.compute_samples_per_period(toneFrequency, sampleRate), 
                                                    seed)
        signalGen2 = SignalGenerator.SignalGenerator(toneFrequency+10e3, toneStrength, resistance, 
                                                    toneSNR, 
                                                    SignalGenerator.compute_samples_per_period(5*toneFrequency, sampleRate), 
                                                    seed)
        signalGen3 = SignalGenerator.SignalGenerator(toneFrequency-70e3, toneStrength, resistance, 
                                                    toneSNR, 
                                                    SignalGenerator.compute_samples_per_period(3*toneFrequency, sampleRate), 
                                                    seed)
    
                                                  
        #del resistance, toneSNR, toneStrength, toneFrequency, seed
        
        totalSamples = int(signalGen1.samplesPerPeriod*totalPeriods) # samples to get full periods
        
        rArray = []
        iArray = []
        cArray = []  
        for s in range(0,totalSamples):
            sigVal1 = signalGen1.toneSignal.next()
            sigVal2 = signalGen2.toneSignal.next()
            sigVal3 = signalGen3.toneSignal.next()
            rArray.append(sigVal1[0]+sigVal2[0]+sigVal3[0])
            iArray.append(sigVal1[1]+sigVal2[1]+sigVal3[1])
            cArray.append(sigVal1[0]+sigVal2[0]+sigVal3[0]+(sigVal1[1]+sigVal2[1]+sigVal3[1])*1j)
            del sigVal1, sigVal2, sigVal3
    
        del totalSamples, totalPeriods
        import SignalMixer
        cArray = SignalMixer.mix_complex_exponential(-toneFrequency,sampleRate,cArray)
        (newfreqArray,newfftArray) = SignalGenerator.compute_ctft(cArray,
                                                                signalGen1.signalSampleRate,
                                                                noShift=False)
        pylab.figure()
        pylab.plot(newfreqArray,numpy.log10(numpy.abs(newfftArray)))
        pylab.title('Resampled Spectra')
        pylab.xlabel('Frequency (Hz)')
        pylab.ylabel('Magnitude of Frequency Response (log scale) (Frequency domain)')
        
        h = create_low_pass_filter(25, 5e3,sampleRate, window='hamming')
        cArray = apply_fir_filter(h,cArray)
        
        
        (newfreqArray,newfftArray) = SignalGenerator.compute_ctft(cArray,
                                                                signalGen1.signalSampleRate,
                                                                noShift=False)
        pylab.figure()
        pylab.plot(newfreqArray,numpy.log10(numpy.abs(newfftArray)))
        pylab.title('Resampled Spectra')
        pylab.xlabel('Frequency (Hz)')
        pylab.ylabel('Magnitude of Frequency Response (log scale) (Frequency domain)')
               
        pylab.show()

if __name__ == "__main__":
    
    SignalFilter().test()
