#
# rps
# mit haystack obs
# 12/2/2014
#
# modified on 6/1/2016 to open/close usb device (to figure out why open fails)
#
# taken from: https://mail.python.org/pipermail/tutor/2008-March/060745.html
#
import sys
import os
import time
import usb


def PrintDevInfo(dev):
    """Print device information."""
    print "Device:", dev.filename
    print "  Device class:",dev.deviceClass
    print "  Device sub class:",dev.deviceSubClass
    print "  Device protocol:",dev.deviceProtocol
    print "  Max packet size:",dev.maxPacketSize
    print "  idVendor:",hex(dev.idVendor)
    print "  idProduct:",hex(dev.idProduct)
    print "  Device Version:",dev.deviceVersion
    print "  Device SerialNumber:",dev.iSerialNumber

    for config in dev.configurations:
        print "  Configuration:", config.value
        print "    Total length:", config.totalLength
        print "    selfPowered:", config.selfPowered
        print "    remoteWakeup:", config.remoteWakeup
        print "    maxPower:", config.maxPower
        for intf in config.interfaces:
            print "    Interface:",intf[0].interfaceNumber
            for alt in intf:
                print "    Alternate Setting:",alt.alternateSetting
                print "      Interface class:",alt.interfaceClass
                print "      Interface sub class:",alt.interfaceSubClass
                print "      Interface protocol:",alt.interfaceProtocol
                for ep in alt.endpoints:
                    print "      Endpoint:",hex(ep.address)
                    print "        Type:",ep.type
                    print "        Max packet size:",ep.maxPacketSize
                    print "        Interval:",ep.interval

				
if __name__ == '__main__':
  reset = 0x01
  endpoint = 8

  # find all of the USB busses
  busses = usb.busses()
  #print "busses=",busses # its an object not a list
  # Find one device
  ii = 0
  rdev = None
  for bus in busses:
    for dev in bus.devices:
        #PrintDevInfo(dev)
        if dev.idVendor == 0x0456 and dev.idProduct == 0xb40d:
            rdev = dev
            ii = ii + 1
            print "Analog Devices ADF435x found"
            #
            # expect two devices, how to distinguish?
            #
            PrintDevInfo(dev)
            print
            print "opening device:",ii
            try:
              handle = dev.open()
            except Exception as eobj:
              print "index=",ii, "Device open exception = ",eobj
            else:
              print "open successful" # there is no close   
        else:
            pass
            #print "dev.idVendor=",hex(dev.idVendor)
            #print "dev.idProduct=",hex(dev.idProduct)
            #print
  if rdev==None:
	print "Could not find a ADF435x\ndev.idVendor == 0x0456"
        print "and dev.idProduct == 0xb40d"       
	sys.exit()
  else:
	dev = rdev
  
  endpoint=8
  current_handle = dev.open()
  #
  # firmware write removed  
  #
  try:
    current_handle.releaseInterface()
  except Exception as eobj:
    pass  # probably requires sudo

  try:
    dev.close()
  except Exception as eobj:
    pass  # raises Exception: 
        # AttributeError: 'Device' object has no attribute 'close'


