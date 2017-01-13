"""test_filter_sin.py tests for windowing effects by low-pass
filtering a sin wave using several windows.

$Id$
"""

import os, os.path, sys

import numpy, scipy
import matplotlib.pylab

import MillstoneScipy.filter

# create and plot signal
time = numpy.arange(1500)
sin = numpy.sin(time*(numpy.pi/100.0))
p1 = matplotlib.pylab.plot(sin[:1000])
matplotlib.pylab.title('raw sin')
matplotlib.pylab.savefig('sin.png')
matplotlib.pylab.clf()

# now apply 50kHz filter
filter = MillstoneScipy.filter.lfilter()
filter.loadFromFilterFile('lp128.0p50.txt')
s50_1 = filter.get_previous_lfilter(sin[0:500]) # will be empty array
s50_2 = filter.get_previous_lfilter(sin[500:1000])
s50_3 = filter.get_previous_lfilter(sin[1000:1500])
s50 = numpy.concatenate((s50_1,s50_2,s50_3))
print('sin after 50kHz filter is len %i' % (len(s50)))
p1 = matplotlib.pylab.plot(s50)
matplotlib.pylab.title('sin after 50kHz filter')
matplotlib.pylab.savefig('50kHz.png')
matplotlib.pylab.clf()
