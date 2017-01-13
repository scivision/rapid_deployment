#
# This is a python DistUtils setup file for the
# PyLpm python module.
#
# This extension tests the use of the Lag Profile
# Matrix code. 
# 
# To build : 
# 
# python setupLpm.py build
#
# To install : 
# 
# python setupLpm.py install
#

# $Id: setupLpm.py 11044 2015-08-14 17:49:51Z brideout $

from distutils.core import setup,Extension
setup(name = "PyLpm", 
      version="1.1", 
      description="modules to support RemoteFileAccess via xmlrpc and RTFileAccess",
      author="Bill Rideout",
      author_email="wrideout@haystack.mit.edu",
      url="http://www.haystack.mit.edu/~brideout/",
      package_dir = {'': 'src'},
      py_modules = ["Lpm"],
      ext_modules = [Extension("PyLpm",["src/PyLpm.c"], libraries=["lpm"])])
