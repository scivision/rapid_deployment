#!/usr/bin/env python
"""
    Code from edy555 github. Modified.

"""
import math
import struct

import fx2load

verbose = False

class ADF4351:
    DIV1 = 0
    DIV2 = 1
    DIV4 = 2
    DIV8 = 3
    DIV16 = 4
    DIV32 = 5
    DIV64 = 6

    REG0 = 0
    REG1 = 1
    REG2 = 2
    REG3 = 3
    REG4 = 4
    REG5 = 5

    REG0_DEF = ( ("unused", 1),
                 ("int", 16),
                 ("frac", 12),
                 ("control", 3, REG0))
    REG1_DEF = ( ("unused", 3),
                 ("phaseadjust", 1),
                 ("prescaler", 1),
                 ("phase", 12, 1),
                 ("mod", 12),
                 ("control", 3, REG1))
    REG2_DEF = ( ("unused", 1,0),
                 ("lownoisespurmode", 2, 0),
                 ("muxout", 3, 0),
                 ("refdoubler", 1, 0),
                 ("rdiv2", 1, 0),
                 ("rcounter", 10),
                 ("doublebuffer", 1, 0),
                 ("chargepump", 4, 0b0111),
                 ("ldf", 1, 0),
                 ("ldp", 1, 0),
                 ("pdpolar", 1, 1),
                 ("powerdown", 1, 0),
                 ("cpthreestate", 1, 0),
                 ("counterreset", 1, 0),
                 ("control", 3, REG2))
    REG3_DEF = ( ("unused", 8),
                 ("bandselectclockmode", 1, 1),
                 ("abp", 1, 0),
                 ("chargecancel", 1, 1),
                 ("unused", 2),
                 ("csr", 1, 0),
                 ("unused", 1),
                 ("clkdivmode", 2, 0),
                 ("clockdivider", 12, 0),
                 ("control", 3, REG3))
    REG4_DEF = ( ("unused", 8, 0),
                 ("feedbackselect", 1, 0),
                 ("rfdividerselect", 3, 0b010),
                 ("bandselectclockdivider", 8, 50),
                 ("vcopowerdown", 1, 0),
                 ("mtld", 1, 0),
                 ("auxoutputselect", 1, 0),
                 ("auxoutputenable", 1, 0),
                 ("auxoutputpower", 2, 0),
                 ("rfoutputenable", 1, 1),
                 ("outputpower", 2, 0b11),
                 ("control", 3, REG4))
    REG5_DEF = ( ("unused", 8, 0),
                 ("ldpinmode", 2, 1),
                 ("unused", 1, 0),
                 ("unused", 2, 0b11),
                 ("unused", 16),
                 ("control", 3, REG5))

    CONFIG_DEFAULT = { "int":160, "frac":0, "mod":2, "rcounter":1,
                       "pdpolar":1, "bandselectclockmode":1 }

    DEFAULT_PDFREQUENCY = 25e6

    def get_config_reg(self, conf, bitdef):
        #ary = bytearray()
        ary = []
        byte = 0
        i = 8
        for k in bitdef:
            s = 1
            d = 0
            if not isinstance(k, str):
                if len(k) > 2:
                    d = k[2]
                s = k[1]
                k = k[0]
            v = conf.get(k)
            if v is None:
                v = d
            while s > 0:
                if i >= s:
                    i -= s
                    s = 0
                    if v:
                        byte |= v << i
                else:
                    n = s - i
                    s -= i
                    i = 0
                    if v:
                        byte |= v >> n
                if i == 0:
                    ary.append(byte & 0xff)
                    byte = 0
                    i = 8
        return ary

    def __init__(self, conf = {}, h = None):
        c = {}
        c.update(ADF4351.CONFIG_DEFAULT)
        c.update(conf)
        self.conf = c
        self.pdfreq = ADF4351.DEFAULT_PDFREQUENCY

        fx2load.openfx2(0x0456,0xb40d)

    def close(self):
        pass
        #self.h.close()

    def setPDFrequency(self, pdfreq):
        self.pdfreq = pdfreq

    def write_reg(self, rd):
        print "write"
        buf=struct.pack ( 'B'*24, *rd)
        fx2load.f.ep_bulk( buf, 0x02, 24)

    def read_reg(self):
        print "read"
        buf='\x00'*24
        fx2load.f.ep_bulk( buf, 0x86, 24)
        rd = struct.unpack ( 'B'*24, buf )
        return rd


    def get_reg(self, conf, bitdef):
        reg = self.get_config_reg(conf, bitdef)
        #head = [0xc1, 32]
        #self.h.write(head + reg)
        if True:
            #print ''.join(["%02x"%x for x in reg])
            print ' '.join([hex(i) for i in reg])

        return reg


    def setup(self):
        conf = self.conf
        print "conf", conf

        print "reg5", ADF4351.REG5_DEF
        r5 = self.get_reg(conf, ADF4351.REG5_DEF)
        print "reg4", ADF4351.REG4_DEF
        r4 = self.get_reg(conf, ADF4351.REG4_DEF)
        print "reg3", ADF4351.REG3_DEF
        r3 = self.get_reg(conf, ADF4351.REG3_DEF)
        print "reg2", ADF4351.REG2_DEF
        r2 = self.get_reg(conf, ADF4351.REG2_DEF)
        print "reg1", ADF4351.REG1_DEF
        r1 = self.get_reg(conf, ADF4351.REG1_DEF)
        print "reg0", ADF4351.REG0_DEF
        r0 = self.get_reg(conf, ADF4351.REG0_DEF)


        rd = r5 + r4 + r3 + r2 + r1 + r0
        print len(rd),' '.join([hex(i) for i in rd])

        self.write_reg(rd)

        rd2 = self.read_reg()

        print len(rd2),' '.join([hex(i) for i in rd2])

    def setIntFrac(self, int, frac):
        self.conf.update({"int": int, "frac":frac })
        self.setup()

    def setFrequency(self, freq):
        f = self.pdfreq
        mod = self.conf["mod"]
        i = int(freq / f)
        diff = freq - i * f
        frac = int(diff * mod / f)
        print f, mod, i, diff, frac
        self.setIntFrac(i, frac)
        print "int:%d frac:%d"%(i, frac)

if __name__ == '__main__':
    import binascii
    from optparse import OptionParser
    def getDivFromFrequency(freq):
        n = int(math.log(4000e6 / freq, 2))
        if n > 6:
            n = 6
        if n < 0:
            n = 0
        return n

    parser = OptionParser(usage="%prog: [options]")
    parser.add_option("-v", dest="verbose", action="store_true", default=False,
                      help="print verbose")
    parser.add_option("-d", dest="divide", type="int",
                      help="rf divider value (note : maps power of two to bit value , 1 -> 0, 2 -> 1, 4 -> 2, etc.)")
    parser.add_option("-e", dest="prescaler", action="store_true", default=False,
                      help="select prescaler 8/9")
    parser.add_option("-f", dest="freq", type="float",default=1000.0E6, help="frequency (Hz)")
    parser.add_option("-l", dest="loop", action="store_true", default=False,
                      help="repeat executions")
    (options, args) = parser.parse_args()
    verbose = options.verbose
    freq = options.freq
    div = options.divide or getDivFromFrequency(freq)
    prescaler = options.prescaler

    try:
        if verbose:
            print "set frequency:%g div:%d"%(freq,div)
        pll = ADF4351(conf={"powerdown":0, "rfdividerselect": div, "prescaler": prescaler})
        #pll.setPDFrequency(1e6)
        #pll.setup()
        print "set freq"
        pll.setFrequency(freq)


        #pll.close()

    except IOError, ex:
        print ex
