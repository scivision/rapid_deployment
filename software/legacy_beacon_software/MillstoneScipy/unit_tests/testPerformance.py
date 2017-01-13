import time
import numpy
import MillstoneScipy.filter

lfilter = MillstoneScipy.filter.lfilter()
lfilter.loadFromFilterFile('/Library/Frameworks/Python.framework/Versions/2.7/config/lp128.0p10.txt')

array1 = numpy.arange(1000, dtype=numpy.complex64)
array2 = numpy.arange(10000, dtype=numpy.complex64)
array3 = numpy.arange(100000, dtype=numpy.complex64)
array4 = numpy.arange(1000000, dtype=numpy.complex64)

lfilter.get_previous_lfilter(array1)
t = time.time()
result = lfilter.get_previous_lfilter(array1)
total = time.time() - t
print('Size 1000 took %f seconds' % (total))

lfilter.get_previous_lfilter(array2)
t = time.time()
result = lfilter.get_previous_lfilter(array2)
total = time.time() - t
print('Size 10000 took %f seconds' % (total))

lfilter.get_previous_lfilter(array3)
t = time.time()
result = lfilter.get_previous_lfilter(array3)
total = time.time() - t
print('Size 100000 took %f seconds' % (total))

lfilter.get_previous_lfilter(array4)
t = time.time()
result = lfilter.get_previous_lfilter(array4)
total = time.time() - t
print('Size 1000000 took %f seconds' % (total))
