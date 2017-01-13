#!/usr/bin/env python

"""
This is a python DistUtils setup file for the
MillstoneScipy Python API.

$Id: setupMillstoneScipy.py 8723 2012-12-14 20:34:44Z brideout $
"""

import os, sys

from distutils.core import setup, Extension

# verify scipy installed
try:
    import scipy
except:
    raise IOError, 'scipy must be installed before MillstoneScipy'

setup(name="MillstoneScipy",
        version="1.1",
        description="Millstone extension to scipy methods",
        author="Bill Rideout",
        author_email="wrideout@haystack.mit.edu",
        url="http://www.haystack.mit.edu/~brideout/",
        #package_dir = {'': 'src'},
        packages=['MillstoneScipy'])
