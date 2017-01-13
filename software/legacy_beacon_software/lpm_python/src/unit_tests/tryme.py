## Automatically adapted for numpy.oldnumeric Apr 08, 2010 by alter_code1.py

# Try out the numtest module

# $Id: tryme.py 3599 2010-04-08 15:28:15Z brideout $

import numpy.oldnumeric as Numeric
import PyLpm
import Lpm

datalength    = 10
fractionality = 1
integerLags   = 4
lpm = {}

strong4 = [
    "++++", "+-++", "++-+", "+--+",
    "+++-", "+-+-", "++--", "+---",
]
mhAltCode16   = [
    "++++++++++++++++", "+-++++---+--++-+",
    "++-+++--+++++-+-", "+--+++++-+--+---",
    "+++-++--+-+----+", "+-+-++++---+--++",
    "++--+++++-+--+--", "+---++-----+-++-",
    "++++-+--+---++--", "+-++-+++--+++++-",
    "++-+-++++---+--+", "+--+-+----+++-++",
    "+++--+++++-+--+-", "+-+--+---++-----",
    "++---+--++-+-+++", "+----+++-++--+-+",
    "+++++-+--+---++-", "+-+++--+++++-+--",
    "++-++--+-+----++", "+--++-+-++++---+",
    "+++-+--+---++---", "+-+-+-+-+-+-+-+-",
    "++--+-+----+++-+", "+---+--++-+-++++",
    "++++---+--++-+-+", "+-++--+-+----+++",
    "++-+--+---++----", "+--+---++-----+-",
    "+++---+--++-+-++", "+-+----+++-++--+",
    "++-----+-++-+++-", "+-----+-++-+++--",
]


codeseq = strong4

# Initialise data structure
for code in codeseq:
    lpm[code] = Lpm.LagProfileMatrix(0, fractionality, integerLags)
    # lpm[code] = Lpm.LagProfileMatrix(9, 2, 4)
    # lpm[code].lpm   = Numeric.ones(100)
    # lpm[code].weight = 1

print "Done initialising"

# Add some data
# import numpy.oldnumeric.random_array as RandomArray
# data = RandomArray.uniform(0., 1., (datalength))
data = (Numeric.pi + Numeric.e*1j) * Numeric.ones(datalength)

for code in codeseq:
    lpm[code].AcfMultiplyAccumulate(data)
    print('Result of AcfMultiplyAccumulate for code %s' % (str(code)))
    lpm[code].PrintLpm()

result = Lpm.ResultLagProfileMatrix(lpm[codeseq[0]])
print('Result of ResultLagProfileMatrix')
result.PrintLpm()

for code in codeseq:
    result.Decode(code, lpm[code])
    print('Result of Decode for code %s' % (str(code)))
    result.PrintLpm()
    
# run all test code in Lpm - commented out test that did not run as imported from CVS
print('about to run test_lpm')
Lpm.test_lpm()
# does not run - Lpm.test_lpm_walsh()
# produces non-deterministic result - Lpm.test_xlpm()
# does not run - Lpm.test_xcfmatrix()
# produces non-deterministic result - Lpm.test_acfmatrix()
print('about to run test_acfmatrix_deterministic')
Lpm.test_acfmatrix_deterministic()
# does not run - Lpm.test_bgsub()
# does not run - Lpm.test_PlotSpectra()





