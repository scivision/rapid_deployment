
import sys
from optparse import OptionParser

from fx2load import *

parser = OptionParser(usage="%prog: [options]")
parser.add_option("--vid", dest="vid", default='0x0456',
                help="select vendor id")
parser.add_option("--pid", dest="pid", default='0xb40d',
                help="select product id")
parser.add_option("--idx", dest="idx", default='0x0',
                help="select index")
parser.add_option("-f", dest="fname", default='build/pll.bix',
				help="select index")

(options, args) = parser.parse_args()

fname = options.fname
vid = int(options.vid,16)
pid = int(options.pid,16)
idx = int(options.idx,16)

print 'load ', fname, 'to device', hex(vid), hex(pid), hex(idx)

openfx2(vid=vid,pid=pid, idx = idx)
reset_bix(fname)
