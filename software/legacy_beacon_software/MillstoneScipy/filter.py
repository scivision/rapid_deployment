"""filter.py is a python module to filter data

$Id: filter.py 2272 2009-11-13 18:02:45Z brideout $
"""

import os, os.path, sys

import numpy
import scipy.signal

class lfilter:
    """ the lfilter class encapsulates the scipy.signal.lfilter method
    """
    
    def __init__(self):
        """set up a default lfilter object
        """
        self._coef = None # this attribute will hold a numpy vector
        self._initCond = None # initial filter conditions
        self._numCoef = 0 # number of coefficients
        
        
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
                self._numCoef =  numLines
                self._coef = numpy.zeros((numLines), dtype=numpy.float64)
            else:
                self._coef[count] = float(line.strip())
                count += 1
                
        # verify
        if count != numLines:
            raise IOError, 'Expected %i values in file %s, got %i' % (numLines, filename, count)
        
        
    def lfilter(self, array, reset=False):
        """lfilter calls the scipy method lfilter to input numpy array, and returns resultant array.
        
        Inputs:
        
            array: a numpy array to apply the filter to
            
            reset: if False (the default), used initial conditions specified by self._initCond.  
                    If True, set initial conditions to all ones before applying filter.
                    
        Ouputs: resultand filters numpy array
                    
        Affects: resets self._initCond 
        """
        if reset or (self._initCond == None):
            self._initCond = numpy.zeros((len(self._coef)-1), dtype=numpy.float64)
            
        result, newInitCond = scipy.signal.lfilter(self._coef, 1.0, array,  zi=self._initCond)
        
        self._initCond = newInitCond
        
        return(result)
        
        
        
if __name__ == '__main__':
    # test
    filtObj = lfilter()
    filtObj.loadFromFilterFile('../../midasw_c/b112.30.dat')
    # test array
    test = numpy.ones((200))
    result = filtObj.lfilter(test)
    print('first: %s' % (str(result)))
    result = filtObj.lfilter(result)
    print('second: %s' % (str(result)))
    result = filtObj.lfilter(result, True)
    print('third with reset: %s' % (str(result)))
    print('done')
                
        
        