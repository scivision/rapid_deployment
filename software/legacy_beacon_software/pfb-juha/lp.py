import numpy as np
import ctypes as C
import numpy
import time

# Stuff to set up access to C routine for median filter (fmedfil)
_lp = np.ctypeslib.load_library('liblp', './Ccode')
_lp.lp.restype = C.c_int
_lp.lp.argtypes = [C.POINTER(C.c_float), C.c_int, C.POINTER(C.c_float), C.c_int, C.POINTER(C.c_float)]
def lp(z, pfb_out, dec, window):
    return _lp.lp(z.ctypes.data_as(C.POINTER(C.c_float)),
                  len(z),
                  pfb_out.ctypes.data_as(C.POINTER(C.c_float)),
                  dec,
                  window.ctypes.data_as(C.POINTER(C.c_float)))
_lp = np.ctypeslib.load_library('liblp', './Ccode')
_lp.pfb.restype = C.c_int
_lp.pfb.argtypes = [C.POINTER(C.c_float), C.c_int, C.POINTER(C.c_float), C.c_int, C.POINTER(C.c_float), C.c_int, C.POINTER(C.c_int), C.c_int]
def pfb(z, n, pfb_out, dec, window, p, channels):
    if len(window) != dec*p:
        print("Error, wrong filter length")
        exit(0)
    return _lp.pfb(z.ctypes.data_as(C.POINTER(C.c_float)),
                   n, 
                   pfb_out.ctypes.data_as(C.POINTER(C.c_float)),
                   dec,
                   window.ctypes.data_as(C.POINTER(C.c_float)),
                   p,
                   channels.ctypes.data_as(C.POINTER(C.c_int)),
                   len(channels))
