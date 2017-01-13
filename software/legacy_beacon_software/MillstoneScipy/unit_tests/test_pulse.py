"""test_pulse.py tests filter.py using filter file lp128.050 and lp128.033

This test is used to make sure filter delays are removed accurately.

$Id: test_filter.py 2331 2009-11-20 21:11:58Z brideout $
"""

import os, os.path, sys

import numpy, scipy
import numpy.random
import matplotlib.pylab

import MillstoneScipy.filter

# create pulse - pulse at beginning at small one near the end
pulse = numpy.ones(1000, numpy.float)
pulse[100:] = 0.0
pulse[980:990] = 1.0


matplotlib.pylab.plot(pulse)
matplotlib.pylab.title('original pulse')
matplotlib.pylab.savefig('org_pulse.png')
matplotlib.pylab.clf()


# now apply 50kHz filter
filter = MillstoneScipy.filter.lfilter()
filter.loadFromFilterFile('lp128.0p50.txt')
filter.get_previous_lfilter(pulse) # first call returns zero length array
s50 = filter.get_previous_lfilter(numpy.zeros(1000, numpy.float))

matplotlib.pylab.plot(s50)
matplotlib.pylab.title('pulse after 50kHz filter')
matplotlib.pylab.savefig('pulse_50kHz.png')
matplotlib.pylab.clf()

# now a second time
filter.get_previous_lfilter(s50)
s502 = filter.get_previous_lfilter(numpy.zeros(1000, numpy.float))

matplotlib.pylab.plot(s502)
matplotlib.pylab.title('pulse after 2nd 50kHz filter')
matplotlib.pylab.savefig('pulse_50kHz_2.png')
matplotlib.pylab.clf()

# now apply 33kHz filter
filter = MillstoneScipy.filter.lfilter()
filter.loadFromFilterFile('lp128.0p33.txt')
filter.get_previous_lfilter(pulse) # first call returns zero length array
s33 = filter.get_previous_lfilter(numpy.zeros(1000, numpy.float))

matplotlib.pylab.plot(s33)
matplotlib.pylab.title('pulse after 33kHz filter')
matplotlib.pylab.savefig('pulse_33kHz.png')
