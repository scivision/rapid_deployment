#!/usr/bin/env python 

import datetime
import glob
import re
import time
import os

import beacon_conf as c

debug = True
def clean_files():
    dl = glob.glob("%s/????.??.??"%(c.datadir))
    dl.sort()

    ## go through all days
    for d in dl:
        date_string = re.search(".*/(..........)",d).group(1)
        dd = datetime.datetime.strptime(date_string,"%Y.%m.%d")
        tnow=time.time()
        tthen = time.mktime(dd.timetuple())
        if (tnow-tthen)/(60.0*60.0*24.0) > c.raw_file_max_age_days:
            if debug:
                print("deleting %s"%(d))
            os.system("rm -f %s/flight-*.bin"%(d))

if __name__ == "__main__":
    clean_files()

