"""
    ADF4351 PLL Load

    Uses a bulk transfer to copy PLL register data to the Cypress microcontroller.
    Register data is then written using serial encoding to the PLL chip. There
    is no feedback beyond reading back the data from the controller.

    Note: The custom firmware must be loaded prior to using this command!
    
"""
import struct

from fx2load import *

# open the Analog devices base unit.
openfx2(0x0456,0xb40d)


def do_bulk():
    # send 100 shorts to ep2

    buf=struct.pack ( 'H'*100, *[i for i in range(100)] )
    f.ep_bulk( buf, 0x02, 1000)

    # read them back out
    buf='\x00'*200
    f.ep_bulk( buf, 0x86, 1000)


    print struct.unpack ( 'H'*100, buf )

[do_bulk() for i in range(3)]
