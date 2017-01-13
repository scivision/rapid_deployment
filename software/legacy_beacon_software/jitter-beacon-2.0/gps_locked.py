#!/usr/bin/env python
##################################################
# Gnuradio Python Flow Graph
# Title: Top Block
# Generated: Mon Oct 24 13:02:28 2011
##################################################

from gnuradio import eng_notation
from gnuradio import gr
from gnuradio import uhd
from gnuradio import window
from gnuradio.eng_option import eng_option
from gnuradio.gr import firdes
from gnuradio.wxgui import forms
from gnuradio.wxgui import waterfallsink2
from grc_gnuradio import wxgui as grc_wxgui
from optparse import OptionParser
import wx

u = uhd.usrp_source(device_addr="addr=192.168.2.2",
		    io_type=uhd.io_type.COMPLEX_FLOAT32,
		    num_channels=1,
		    )
print u.get_mboard_sensor("gps_locked", 0)
print u.get_mboard_sensor("gps_gpgga", 0)
#print u.get_mboard_sensor("gps_gpgmc", 0)
#print u.get_mboard_sensor("gps_gpgsa", 0)
