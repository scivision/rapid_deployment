"""
beacon_analysis.py is a script to calculate higher-level data from
processed beacon receiver voltages.

Possible parameters to calculate are S4, sigma-phi and absolute and
relative TEC. Each parameter has several different calculation
options; see the description for the __init__ function and each
function for more details. Additionally, the class can calculate
zenith TEC from line-of-sight TEC when given an ephemeris for the
satellite. Additionally, the class provides functions that can find
the power of a signal given its complex voltage and receiver
resistance.

Usage example:
    import HighLevelProc
    
    s4Array = HighLevelProc.s4(voltages, 100)
    sigmaPhiArray = HighLevelProc.sigma_phi(voltages, 100, ranges,
                                            frequencies)
    tecArray = HighLevelProc.slant_path_tec(vhfVoltages, vhfFreqs,
                                            uhfVoltages, uhfFreqs,
                                            lVoltages, lFreqs,
                                            ranges, 0.1, 0.1, 0.1) 
     
    
$Id: HighLevelProc.py 142 2007-08-10 19:54:09Z damiana $
"""

# TODO: - implement error in TEC calculation
#       - find test data and see if this actually works


import math

import numpy
import scipy
import scipy.interpolate

# speed of light in m/s
__c = 299792458.
# define NaN so that it can be returned
__nan = 1e300*1e300-1e300*1e300


def is_nan(tuple):
    """is_nan checks the values of each tuple to see any contain
    the IEEE-754 standard NaN value. If at least one does, True
    is returned; otherwise, False is returned. Algorithm by
    Bill Rideout.
    """

    for num in tuple: 
        if ( (num < 10.0) == False and (num > -10.0) == False ):
            return True
        
    # no NaN values found
    return False



def power(complexVoltage, resistance):
    """ power is a public funciton that calculates the signal power
    given a complex voltage and resistance.
    
    Inputs:
        complexVoltage - (V) a complex type containing the voltage
        resistance - (ohms) resistance of the receiver 

    Returns: the power of the signal (in watts)

    Algorithm: Calculate the square of the magnitude of the voltage by
             finding the sum of squares of the real and imaginary
             voltages. Find the power by dividing the voltage squared
             by the resistance.

    Exceptions: throws SyntaxError if a zero or negative resistance
                is given
    """
    # check for proper syntax

    if (resistance <= 0):
        raise SyntaxError,'Incorrect resistance of '+resistance+' ohms'

    # find voltage squared
    vsq = math.pow(abs(complexVoltage),2)

    # return power
    return vsq / resistance

def find_phases(voltageData, frequencies, ranges, phaseTolerance=math.pi/4):
    """ find_phases is a public function that finds the phases of
        a list of complex voltages, taking into account the change
        in range between the transmitted and receiver.

        Inputs:
            voltageData - a list of complex types, representing
                              the signal voltage
            frequencies - (Hz) a list of floats with the frequencies
                               at which the signal was received
            ranges - (m) a list of ranges between the satellite
                         and receiver
            phaseTolerance - (rad) the tolerance allowed to detect
                             a wrap-around of phase. Setting this
                             tolerance too low will result in impossible
                             sigma-phi values. 
            
        Returns - a list of floats containing the corrected phases
                  of the signal in radians, with a range of -pi to pi.

        Algorithm: The change in phase is calculated as follows:

              delta-phi = 2*pi * delta-S * f / c

              where phi is the phase, S is the range to the satellite,
              f is the received frequency and c is the speed of light.
              
              The phases are then redone to try to normalize wrap-arounds
              so that the phase range can exceed -pi to pi.

        Exceptions: none
    """
    # find normalized phases for all points
    phases = []
    for i in range(len(voltageData)):
        phases.append( math.atan2(voltageData[i].imag,voltageData[i].real) )
            
        # find normalization (if first element, only normalize to 0)
        if (i == 0):
            norm = 0.
        else:
            deltaS = ranges[i] - ranges[i-1]
            norm += 2 * math.pi * deltaS * frequencies[i] / __c 
            # shift normalization factor to a range of -pi to pi
            # (first rough, then smooth in order to save computation time)
            while (norm > 10*math.pi):
                norm -= 20*math.pi
            while (norm < -10*math.pi):
                norm += 20*math.pi
            while (norm > math.pi):
                norm -= 2*math.pi
            while (norm < -math.pi):
                norm += 2*math.pi

        # normalize and shift to a range of -pi to pi
        phases[i] -= norm
        while (phases[i] < -math.pi):
            phases[i] += 2*math.pi
        while (phases[i] > math.pi):
            phases[i] -= 2*math.pi

    # normalize out phase wrap-arounds        
    phases2 = [phases[0]]
    for i in range(len(phases)-1):
        if (not is_nan( (phases[i],phases[i+1]) )):
            diff = phases[i+1]-phases[i]
            if (abs(diff) > (2*math.pi - phaseTolerance)):
                diff = (diff + math.pi) % (2*math.pi) - math.pi
            phases2.append(phases2[-1] + diff)
        else:
            phases2.append(__nan)

    return phases2

def detrend_phases(phases, numPoints=25):
    """ detrend_phases fits a spline to a set of phases to cancel out any
        long-term drifts. It can be used directly after calling find_phases
        to make the phases easier to plot.

        Inputs:
            phases - a list/array of floats representing phases
            numPoints - the number of points that is used for each point
                        of the fitting. A smaller number is more accurate
                        but much more computationally intensive.

        Returns:
            a list of floats representing the corrected phases

        Exceptions: none
    """

    # find smoothing spline in order to cancel out long-term drift
    sparse = []
    for i in range(len(phases)):
        if (i%numPoints == 0):
            sparse.append(phases[i])
    tck = scipy.interpolate.splrep(numpy.arange(0,len(phases),numPoints),
                                   sparse)
    y = scipy.interpolate.splev(range(len(phases)),tck)
    for i in range(len(phases)):
        phases[i] = phases[i] - y[i]
        
    return phases

def __find_tec_wraps(tecList, tecAmbiguity, tecThreshold,
                     tecTolerance, maxBreak, minLookIndices,
                     maxLookIndices):
    """ __find_tec_wraps is a private function that takes in a list of
    TEC values and attempts to compensate for phase wrapping by
    detecting wrap-arounds and shifting the values by the appropriate
    number of wrap-arounds.

    Inputs:
        tecList: a list of floats representing a time series of TEC values
        tecAmbiguity: a float representing the specified TEC ambiguity. Must
                      be in the same units as tecList.
        tecThreshold: a float with the specified threshold at which
                      a discontinuity (which might be a wrap-around event)
                      is detected. Must be in same units as tecList.
        tecTolerance: a float with the specified tolerance so that the
                      function knows when to recognize a wrap-around.
                      Must be in the same units as tecList.
        maxBreak:    an integer with the maximum number of samples
                     containing NaN that are allowed before the TEC
                     is reset to zero.
        minLookIndices: an integer with the minumum allowed number of
                        indices to look back for a linear fit
        maxLookIndices: an integer representing the maximum number of 
                        indices to look back for a linear fit
    Returns: An array of floats representing the TEC values with
             wrap-arounds detected.

    Algorithm: Looks for candidate points that are discontinuous from
               previous points by the specified threshold, and does
               a linear fit to see if the new point will fit the ongoing
               trend (within tolerance) if the wrap-around is corrected

    Exceptions: none   
    """

    # copy over TEC values to a new array so that tecList does not
    # change
    tecArray = numpy.array(tecList)
    
    # find the first non-nan point
    i = 0
    while (i < len(tecArray)-1 and is_nan( (tecArray[i],) )):
        i += 1

    # the value to shift the TEC by (taking into account
    # wrap-arounds and normalization to zero)
    shift = -tecArray[i]
    
    while (i < len(tecArray)-1):
        # find next point to compare against
        next = i+1
        while (next < len(tecArray)-1 and is_nan( (tecArray[next],) )):
            next += 1
        # no more valid points to compare against: finish
        if (is_nan( (tecArray[next],) )):
            break
        
        # span is too big: give up, reset TEC to 0 and go to the next point
        if (next-i > maxBreak):
            shift = -tecArray[next]
            tecArray[next] = 0
            i = next
            continue

        # adjust the next point by the necessary number of wraps
        tecArray[next] += shift

        # point is within threshold: nothing special, move along
        if (abs(tecArray[next]-tecArray[i]) < tecThreshold):
            i = next
            continue

        # find minimum index for regression (ideally, this will
        # be i - maxLookIndices)
        min = i
        while (min > 0 and not is_nan( (tecArray[min-1],) ) and i-min < maxLookIndices):
            min -= 1

        # if not looking far back enough, give up and go to the first
        # non-nan point
        if (i-min < minLookIndices):
            i = next
            continue

        # do a linear regression
        (a,b) = scipy.polyfit(numpy.arange(min,i+1,1), # indices
                                          tecArray[min:i+1], # TEC values
                                          1) # linear regression
        
        # see if the next point follows the trend if shifted up or down
        # by one wrap
        if (abs(a*next+b + tecAmbiguity - tecArray[next]) < tecTolerance):
            shift -= tecAmbiguity
            tecArray[next] -= tecAmbiguity
            i = next
            continue

        elif (abs(a*next+b - tecAmbiguity - tecArray[next]) < tecTolerance):
            shift += tecAmbiguity
            tecArray[next] += tecAmbiguity
            i = next
            continue

        # no go: don't wrap
        else:
            i = next
            continue

    # return the new array
    return tecArray

def s4(voltageData, windowSize):
    """ s4 is a public function that finds the S4 index for a list of
    complex voltages, with a given number of points to average over
    for each index.

    Inputs:
        voltageData - a list of complex types, each containing the
                      signal voltage
        windowSize - the number of points to average over
                              to find each index

    Returns: a list of 2-tuples, each containing the S4 index and its
             associated standard error. The length of this list
             depends on windowSize. If windowSize does not divide
             evenly into the length of the voltage list, then the last
             few elements of the list will be still calculated but
             will have a larger associated error. Thus, the number of
             S4 points will equal the lenght of voltageData divided by
             windowSize and rounded to the next larger integer.

    Algorithm: Iterate through the list of voltage points. Find
               the average of signal power squared for a number of
               points, and the square of the average. Find S4
               according to the following formula:

               S4 = ((<P^2> - <P>^2) / (<P>^2))^(1/2)

               Find the standard error according to the
               following formula:

               S4_M = S4 / (n)^(1/2)

               where n is the number of data points used to find
               the S4 value.

    Exceptions: SyntaxError is raised if windowSize is 0 or negative
    """
    # check syntax
    if (windowSize <= 1):
        raise SyntaxError,'Incorrect window size of %s' % windowSize

    # find number of sequences for which to find S4
    s = int(len(voltageData) / windowSize)
    if ((len(voltageData) % windowSize) != 0):
        s += 1

    # initialize S4, error list
    s4list = []

    # loop through voltage series
    for i in range(s):
        sum = 0
        sumSquares = 0

        # find how many elements there are in the series (answer is either
        # windowSize or something smaller if at end of series)

        if (i == s-1):
            elem = len(voltageData) - i * windowSize
        else:
            elem = windowSize
        
        # get all voltages for window - discard those with a nan value
        windowVoltages = []
        for j in range(elem):
            if (not is_nan( (voltageData[i*windowSize+j].real, voltageData[i*windowSize+j].imag) )):
                windowVoltages.append(voltageData[i*windowSize+j])
        windowSamples = len(windowVoltages)

        # need at least two samples to find S4; if not enough,
        # give up for this window
        if (windowSamples < 2):
            s4list.append( (__nan, __nan) )
        else:
            # loop through each element in series
            for j in range(windowSamples):
                # find sum of power and sum of squares (resistance
                # is arbitrary when finding S4, as values are normalized)
                sum += power(windowVoltages[j],1)
                sumSquares += math.pow(power(windowVoltages[j],1),2)

            # find <P^2> and <P>^2
            averageSquarePower = sumSquares / windowSamples
            averagePowerSquared = math.pow(sum / windowSamples, 2)
 
            # find S4 and standard error
            s4 = math.sqrt(abs((averageSquarePower - averagePowerSquared) /
                               averagePowerSquared))
            standardError = s4 / math.sqrt(windowSamples)

            # append to list that will be returned
            s4list.append( (s4, standardError) )

    # return final list
    return s4list

def sigma_phi(voltageData, windowSize, ranges,
              frequencies, phaseTolerance= math.pi):
    """ sigma_phi is a public function that finds the standard
    deviation of the phase shift of the signal given lists of complex
    voltages, frequencies and ranges to the satellite. The number of
    points to average and correct over for each sigma-phi measurement
    (window size) is chosen by the user.
    
    Inputs:
        voltageData - a list of complex types, each containing
                      the signal voltage
        windowSize - the number of points to look at at a time
                     to get a sigma phi measurement
        ranges - a list of ranges (in meters) to the satellite
                 that corresponds to the voltage measurements. Must
                 have the same length as voltageData.
        frequencies - a list of frequencies (in Hz) at which the
                      satellite signal was found that corresponds
                      to the voltage measurements. Must have the
                      same length as voltageData.

    Returns: a list of 2-tuples, each containing the sigma-phi index
              in radians and its associated standard error. The length
              of this list depends on windowSize. If windowSize does
              not divide evenly into the length of the voltage list,
              then the last few elements of the list will be still
              calculated but will have a larger associated
              error. Thus, the number of sigma-phi points will equal
              the lenght of voltageData divided by windowSize and
              rounded to the next larger integer.

    Algorithm: The phases of all the data points are
               found. Additionally, they are normalized around
               zero phase and by taking the change in range and
               frequency into account. Each data point is
               shifted in phase by the following formula:

               delta-phi = 2*pi * delta-s * f / c

               where delta-phi is the difference in phase in radians,
               delta-s is the difference in satellite range from the
               previous point, f is the frequency and c is the speed
               of light. Additionally, a quadratic least-squares
               regression is performed on the points in the window,
               and the curve is subracted from the data.  The standard
               deviation of the phase difference list is then found
               according to the following formula:

               sigma-phi = (<phi^2> - <phi>^2)^(1/2)

               where <phi^2> is the average of the squares of the
               phase shifts and <phi>^2 is the square of the
               average of the phase shifts. The standard error is
               calculated according to the following formula:

               sigma-phi_M = sigma-phi / (n)^(1/2)

               where n is the number of data points used to find
               the sigma-phi value.

    Exceptions: SyntaxError is thrown if windowSize is not positive
                or if the lists have different lengths.
    """

    # check syntax
    if ((len(voltageData) != len(ranges)) or
        (len(voltageData) != len(frequencies))):
        raise SyntaxError,'Voltages, ranges and frequencies must have the same length'
    if (windowSize <= 1):
        raise SyntaxError,'Incorrect window size'

    phases = find_phases(voltageData, frequencies, ranges)
    
    # find number of sequences for which to find sigma-phi
    s = int(len(phases) / windowSize)
    if ((len(phases) % windowSize) != 0):
        s += 1

    # initialize sigma-phi, error list
    sigmaPhiList = []

    # loop through voltage series
    for i in range(s):
        sum = 0
        sumSquares = 0

        # find how many elements there are in the series (answer is either
        # windowSize or something smaller if at end of series)
        if (i == s-1):
            elem = len(phases) - i * windowSize
        else:
            elem = windowSize

        # copy over phase data in window to a new list so that
        # a quadratic regression can be made; keep note of
        # subscripts to help the regression
        windowPhases = []
        subscripts = []
        for j in range(elem):
            if (not is_nan( (phases[i*windowSize+j],) )):
                windowPhases.append(phases[i*windowSize+j])
                subscripts.append(j)
        windowSamples = len(windowPhases)
        
        # give up if there are less than two phase values
        # for the window
        if (windowSamples < 2):
            sigmaPhiList.append( (__nan, __nan) )

        else:
            # do quadratic regression to detrend phase
            coeffs = scipy.polyfit(subscripts, windowPhases, 2)   
            # detrend phase
            for j in range(windowSamples):
                # find point on quadratic regression
                y = coeffs[0]*math.pow(subscripts[j],2) + coeffs[1]*subscripts[j] + coeffs[2]
                windowPhases[j] = windowPhases[j] - y

            # loop through each element in series
            for j in range(windowSamples):
                # find sum of phase and sum of squares of phase
                sum += windowPhases[j]
                sumSquares += pow(windowPhases[j],2)

            # find <phi^2> and <phi>^2
            averageSquarePhase = sumSquares / windowSamples
            averagePhaseSquared = math.pow(sum / windowSamples,2)

            # find sigma-phi and standard error
            sigmaPhi = math.sqrt(abs(averageSquarePhase-averagePhaseSquared))
            standardError = sigmaPhi / math.sqrt(windowSamples)

            # append to list that will be returned
            sigmaPhiList.append( (sigmaPhi, standardError) )

    # return final list
    return sigmaPhiList

def slant_path_tec(vhfVoltageData, vhfFrequencies, uhfVoltageData,
                 uhfFrequencies, lVoltageData, lFrequencies,
                   ranges, vhfPhaseError, uhfPhaseError, lPhaseError,
                   phaseTolerance=math.pi,
                   vulWrap=(6e16, 0.8e16, 1000, 10, 1200),
                   vuWrap=(6e14, 2e14, 250, 3, 50),
                   vlWrap=(6e14, 1e14, 250, 3, 50),
                   ulWrap=(2e15, 3e14, 250, 3, 100)):

    """ slant_path_tec is a public function that finds the two- and
        three-frequency slant-path TEC based on voltage data for VHF, UHF
        and L-band frequencies taken concurrently. Additionally, it
        finds the associated errors.

    Inputs:
        vhfVoltageData, uhfVoltageData, lVoltageData:
            lists of 2-tuples, each containing the real and imaginary
            parts of the VHF (150 MHz), UHF (400 MHz) and L-band
            (1067 MHz) voltage signals, respectively. All three lists
            should be of the same length. If voltage data is unavailabe
            at a specific time for one frequency, then its
            value should be (nan, nan).
            
        vhfFrequencies, uhfFrequencies, lFrequencies:
            (Hz) lists of exact frequencies at which the signal was
            received for each frequency.
            
        ranges (m): list of floats of ranges to satellite

        vhfPhaseError, uhfPhaseError, lPhaseError:
           (rad) the phase errors, used to calculate the TEC errors
        
        phaseTolerance (rad): tolerance to finding a wrap in phases
        
        vulWrap, vuWrap, vlWrap, ulWrap: Each is a 4-tuple containing
            parameters that define how to look for TEC wrap-arounds.
            Each 5-tuple contains the following variables:
            tecThreshold (e- / m^2): a float with the specified
                threshold at which a discontinuity (which might be a
                wrap-around event) is detected
            tecTolerance (e- / m^2): a float with the specified
                tolerance so that the function knows when to recognize
                a wrap-around
            maxBreak: an integer with the maximum number of samples
                containing NaN that are allowed before the TEC is
                reset to zero
            minLookIndices: an integer with the minumum allowed number of
                indices to look back for a linear fit
            maxLookIndices: an integer representing the maximum number
                of indices to look back for a linear fit
                          
    Returns: an array of size (# samples) x 4 x 2, each row containing
             the three-frequency TEC and error, VHF & UHF TEC and
             error, VHF & L-band TEC and error, and UHF & L-band TEC
             and error. Units are in e- / m^2.

    Algorithm: The phases of the signals are found. Then, the
    differences in phase between VHF and UHF, VHF and L-band and UHF
    and L-band are found. The TEC is then computed from these
    differences using several different algorithms, depending on which
    sets of frequencies are available. If all three are available,
    then three-frequency TEC can be found and is computed from the following
    equation:

    TEC = 8.3165e16 * ( (P_13 * 7 - P_12 * 8) mod 1)
    
    where P_13 is the difference in phase between VHF and L-band in
    wavelengths and P_12 is the difference between VHF and UHF in
    wavelengths. Phase differences for two different frequencies are
    calculated as follows:

    P_ab = P_a - P_b * f_a / f_b

    where f_a / f_b is the ratio of frequencies used. Two-frequency
    TEC is calculated as follows:

    VHF and L-band: TEC = 1.134e15 * delta-P
    VHF and UHF: TEC = 1.256e15 * delta-P
    UHF and L-band: TEC = 3.906e15 * delta-P

    where delta-P is the phase difference between the frequencies in
    wavelengths.

    The error returned for all four TEC values is simply the phase error
    propagated through the TEC equations. 

    Wrap-arounds in phase between two successive samples are detected
    (as long as the phase shift is within the tolerance specified by
    phaseTolerance) and multiples of whole cycles are added to the
    TEC to attempt to reconstruct the absolute TEC as best as
    possible. If less than two frequencies are available, then the
    TEC values are marked as (nan, nan). If two much of a pause
    is present in any of the four TEC sequences, then it will reset
    to zero.

    Exceptions: throws SyntaxError if the voltage data lists are of
                differing lengths or if phaseTolerance is negative
    """
    # check for syntax
    if (len(vhfVoltageData) != len(uhfVoltageData) or
        len(vhfVoltageData) != len(lVoltageData)):
        raise SyntaxError, 'Voltage data lists are of mismatching length'
    if (phaseTolerance < 0.):
        raise SyntaxError, 'Incorrect phase tolerance of %s' % phaseTolerance

    samples = len(vhfVoltageData)

    vPhases = find_phases(vhfVoltageData, vhfFrequencies, ranges, phaseTolerance)
    uPhases = find_phases(uhfVoltageData, uhfFrequencies, ranges, phaseTolerance)
    lPhases = find_phases(lVoltageData, lFrequencies, ranges, phaseTolerance)

    # find differential phases and errors
    vuPhases = numpy.zeros(samples)   # VHF and UHF
    vlPhases = numpy.zeros(samples)   # VHF and L-band
    ulPhases = numpy.zeros(samples)   # UHF and L-band

    for i in range(samples):        
        vuPhases[i] = (1/(2*math.pi)*(vPhases[i]-uPhases[i]*vhfFrequencies[i]/uhfFrequencies[i]))
        vlPhases[i] = (1/(2*math.pi)*(vPhases[i]-lPhases[i]*vhfFrequencies[i]/lFrequencies[i])) 
        ulPhases[i] = (1/(2*math.pi)*(uPhases[i]-lPhases[i]*uhfFrequencies[i]/lFrequencies[i]))

    vuError = 1/(2*math.pi)*math.sqrt(math.pow(vhfPhaseError,2)+math.pow(uhfPhaseError,2))
    vlError = 1/(2*math.pi)*math.sqrt(math.pow(vhfPhaseError,2)+math.pow(lPhaseError,2))
    ulError = 1/(2*math.pi)*math.sqrt(math.pow(uhfPhaseError,2)+math.pow(lPhaseError,2))
    
    # arrays of TEC values and errors
    vulTECValues = numpy.zeros(samples)
    vuTECValues = numpy.zeros(samples)
    vlTECValues = numpy.zeros(samples)
    ulTECValues = numpy.zeros(samples)

    vulTECErrors = numpy.ones(samples) * 8.3165e16 * math.sqrt(math.pow(7*vlError,2)+math.pow(8*vuError,2))
    vuTECErrors = numpy.ones(samples) * 1.2995e15 * vuError
    vlTECErrors = numpy.ones(samples) * 1.1392e15 * vlError
    ulTECErrors = numpy.ones(samples) * 3.4652e15 * ulError

    # iterate through phase samples and attempt finding all four
    # different TEC values
    for i in range(samples):
        # first try three-frequency TEC
        if (not is_nan( (vuPhases[i],vlPhases[i]) )):
            vulTECValues[i] = 8.3165e16*((vlPhases[i]*7.-vuPhases[i]*8.) % 1.)
        else:
            vulTECValues[i] = __nan
            vulTECErrors[i] = __nan

        # VHF & UHF
        if (not is_nan( (vuPhases[i],))):
            vuTECValues[i] = 1.2995e15 * (-vuPhases[i] % 1.)
        else:
            vuTECValues[i] = __nan
            vuTECErrors[i] = __nan
        
        # VHF & L-band
        if (not is_nan( (vlPhases[i],))):
            vlTECValues[i] = 1.1392e15 * (-vlPhases[i] % 1.)
        else:
            vlTECValues[i] = __nan
            vlTECErrors[i] = __nan

        # UHF & L-band
        if (not is_nan( (ulPhases[i],))):
            ulTECValues[i] = 3.4652e15 * (-ulPhases[i] % 1.)
        else:
            ulTECValues[i] = __nan
            ulTECErrors[i] = __nan

    # normalize to zero and try to find wrap-arounds
    vulTECValues = __find_tec_wraps(vulTECValues, 8.3165e16, vulWrap[0],
                                    vulWrap[1], vulWrap[2], vulWrap[3],
                                    vulWrap[4])
    vuTECValues = __find_tec_wraps(vuTECValues, 1.2995e15, vuWrap[0],
                                   vuWrap[1], vuWrap[2], vuWrap[3], vuWrap[4])
    vlTECValues = __find_tec_wraps(vlTECValues, 1.1392e15, vlWrap[0],
                                   vlWrap[1], vlWrap[2], vlWrap[3], vlWrap[4])
    ulTECValues = __find_tec_wraps(ulTECValues, 3.4652e15, ulWrap[0],
                                   ulWrap[1], ulWrap[2], ulWrap[3], ulWrap[4])

    # put data from all eight arrays into an array of size
    # samples x 8
    tecArray = numpy.zeros((samples,4,2))
    for i in range(samples):
        tecArray[i] =  ( (vulTECValues[i],vulTECErrors[i]),
                         (vuTECValues[i],vuTECErrors[i]),
                         (vlTECValues[i],vlTECErrors[i]),
                         (ulTECValues[i],ulTECErrors[i])  )
    return tecArray
    
def tec_map(el):
    """ tec_map is a public function that returns the mapping function
    from slant-path to vertical TEC. The slant-path TEC divided by
    this mapping function gives the vertical TEC.

    Inputs:
        el: (deg) The elevation of the satellite with relation
            to the receiver

    Returns: the mapping function

    Algorithm: The following equation is used:

                   M = 1 / sqrt(1 - (fit * cos(el))^2)

               where fit is the fitting parameter.
    """
    fit = 0.95

    return 1. / math.sqrt(1. - pow(fit * math.cos(math.radians(el)), 2))
