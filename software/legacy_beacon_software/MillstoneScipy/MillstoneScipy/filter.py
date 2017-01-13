"""filter.py is a python module to filter data

$Id: filter.py 8778 2012-12-28 19:42:22Z brideout $
"""

import os, os.path, sys
import copy

import numpy
import scipy.signal

class lfilter:
    """ the lfilter class encapsulates the scipy.signal.lfilter method
    """
    
    def __init__(self, coef=None):
        """set up a default lfilter object.
        
        Input:  coef - numpy array specifying filter coefficients.  Default = None.
        
        Returns:  None
        
        Affects:  Sets self._coef to coefficients if specified.
        
        Exceptions:  Raises AttributeError if coef is not a numpy array.
        """
        if coef != None:
            if type(coef) != numpy.ndarray:
                raise AttributeError, 'coefficients must be a numpy array'
            # verify odd length
            if len(coef) % 2 != 1:
                raise ValueError, 'len of coef must be odd, not length %i' % (len(coef))
            self._delay = (len(coef)-1)/2
            
        else:
            self._delay = None
        self._coef = coef
        self._firstArray = None # first cached array
        
        
    def loadFromFilterFile(self, filename):
        """loadFromFilterFile loads self._coef from a filter file as found in
        prototypes/midasw_c/*.dat
        
        Each file has first line as number of filters, followed by real values,
        one per line
        
        Input: filename - path to filter file
        
        Returns: None
        
        Raises IOError if problem with filter file
        """
        f = open(filename)
        lines = f.readlines()
        f.close()
        
        numLines = None
        count = 0
        
        for i in range(len(lines)):
            line = lines[i]
            if len(line.strip()) == 0:
                continue
            if numLines == None:
                numLines = int(line.strip())
                self._coef = numpy.zeros((numLines), dtype=numpy.float64)
            else:
                self._coef[count] = float(line.strip())
                count += 1
                
        # verify
        if count != numLines:
            raise IOError, 'Expected %i values in file %s, got %i' % (numLines, filename, count)
        
        # verify odd length
        if len(self._coef) % 2 != 1:
            raise ValueError, 'len of coef must be odd, not length %i' % (len(self._coef))
        
        self._delay = (len(self._coef)-1)/2
        
        
    def get_previous_lfilter(self, array, reset=False):
        """get_previous_lfilter calls the scipy method lfilter.  It is a
        stateful method designed to allow repeated calls with windowed data.
        It always returns the lfilter result for the previous input, except
        for the first call, which returns a zero length array.
        It deals with edge effects by caching two
        previous arrays, so that array being filtered is always surrounded by
        previous and following arrays.  A zero array is used before the first
        input array.
        
        Inputs:
        
            array: a numpy complex64 array to apply the filter to
                    
        Outputs: filtered numpy array of the input array passed in by PREVIOUS call
                    
        Affects: resets self._firstArray and self._secArray
        
        Exception: if len(array) < len(filter)/2 + 1
        """
        if self._coef == None:
            raise ValueError, 'Coefficients not specified'
        
        if reset or (self._firstArray == None):
            self._firstArray = numpy.zeros((len(self._coef)), dtype=numpy.complex64)
            self._secArray = array
            return(numpy.zeros((0), dtype=numpy.complex64))
        
        concatArray = numpy.concatenate((self._firstArray, self._secArray, array[0:len(self._coef)]))
            
        #result = scipy.signal.lfilter(self._coef, 1.0, concatArray)
        result = numpy.convolve(concatArray, self._coef, mode='full')
        
        startIndex = self._delay + len(self._firstArray)
        endIndex = startIndex + len(self._secArray)
        
        if endIndex > len(result):
            raise ValueError, 'len of arrays too short to filter'
        
        finalResult = result[startIndex:endIndex]
        
        # reset state variables
        self._firstArray = copy.deepcopy(self._secArray)[-1*len(self._coef):]
        self._secArray =  copy.deepcopy(array)
        
        return(finalResult)
    
    
    def getCoefficients(self):
        """return coefficients from self._coef"""
        return(self._coef)
    
        
        
        
if __name__ == '__main__':
    # test
    filtObj = lfilter()
    filtObj.loadFromFilterFile('../unit_tests/lp128.0p33.txt')
    # test array
    test = numpy.ones((200))
    result = filtObj.get_previous_lfilter(test)
    print('first (should be empty): %s' % (str(result)))
    result = filtObj.get_previous_lfilter(test)
    print('second: %s' % (str(result)))
    print('done')
                
        
        