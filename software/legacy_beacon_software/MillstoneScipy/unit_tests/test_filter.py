"""test_filter.py tests filter.py using filter files lp128.0p33 and lp128.050

$Id: test_filter.py 7943 2012-10-17 17:32:40Z brideout $
"""

import os, os.path, sys

import numpy, scipy
import numpy.random
import matplotlib.pylab

import MillstoneScipy.filter

# create white noise
white = numpy.random.random_sample(1000) - 0.5

# first show fft of unfiltered white data
fft1 = scipy.fft(white)
print('fft of white noise')
p1 = matplotlib.pylab.plot(fft1)
matplotlib.pylab.title('fft of white noise')
matplotlib.pylab.savefig('white_noise.png')
matplotlib.pylab.clf()

# fftshifted x axis to make dc show in the middle
freqs = numpy.fft.fftfreq(1000, 0.00001)
xarr = numpy.fft.fftshift(freqs)

# now apply 50kHz filter
filter = MillstoneScipy.filter.lfilter()
filter.loadFromFilterFile('lp128.0p50.txt')
print('Coefficients are %s' % (filter.getCoefficients()))
filter.get_previous_lfilter(white)  # first call returns zero length array
s50 = filter.get_previous_lfilter(white)
f50 = scipy.fft(s50)
print('fft after 50kHz filter')
p1 = matplotlib.pylab.plot(xarr, numpy.fft.fftshift(f50))
matplotlib.pylab.title('fft after 50kHz filter')
matplotlib.pylab.savefig('50kHz.png')
matplotlib.pylab.clf()
filter = None

# now apply 33kHz filter
filter = MillstoneScipy.filter.lfilter()
filter.loadFromFilterFile('lp128.0p33.txt')
filter.get_previous_lfilter(white)  # first call returns zero length array
s33 = filter.get_previous_lfilter(white)
f33 = scipy.fft(s33)
print('fft after 33kHz filter')
p1 = matplotlib.pylab.plot(xarr, numpy.fft.fftshift(f33))
matplotlib.pylab.title('fft after 33kHz filter')
matplotlib.pylab.savefig('33kHz.png')
