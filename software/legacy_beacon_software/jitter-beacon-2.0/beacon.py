#!/usr/bin/env python
#
# Beacon satellite receiver 
#
# (c) 2014 Juha Vierinen
# 
#
import optparse, os, time, thread, datetime
import subprocess, math

import beacon_conf as c
import beacon_predict
import beacon_fetchtle

beacon_fetchtle.proxies = c.proxies

def deg2degminsec(x):
    deg = math.floor(x)
    minutes = math.floor(60.0*(x-deg))
    sec = math.floor(3600.0*(x-deg-minutes/60.0))
    coord = "%d:%d:%d"%(int(deg),int(minutes),int(sec))
    print(coord)
    return(coord)

def fetch_ephemeris(state):
    # check now and one hour from now
    times_to_check = [datetime.datetime.utcfromtimestamp(time.time()),
                      datetime.datetime.utcfromtimestamp(time.time()+3600.0)]
    state["last_tle_fetch"]=True
    for tnow in times_to_check:
        #   datetime.datetime(2014, 5, 9, 17, 32, 53, 818594)
        date_string = "%04d.%02d.%02d"%(tnow.year,tnow.month,tnow.day)

        # fetch new tle files if needed
        if not os.path.exists("%s/%s/beacon.tle"%(c.datadir,date_string)):
            print("Fetching files")
            os.system("mkdir -p %s/%s"%(c.datadir, date_string))
            success = beacon_fetchtle.download_tle_files("%s/%s"%(c.datadir, date_string))
            state["last_tle_fetch"]=success

        if not os.path.exists("%s/%s/passes.txt"%(c.datadir,date_string)):
            print("Predicting passes for %s"%(date_string))
            os.system("mkdir -p %s/%s"%(c.datadir, date_string))
            beacon_predict.write_pass_files(datestr="%04d/%02d/%02d 00:00:00"%(tnow.year,tnow.month,tnow.day),
                                            data_dir="%s/%s"%(c.datadir,date_string),
                                            latitude=deg2degminsec(c.station_latitude),
                                            longitude=deg2degminsec(c.station_longitude),
                                            altitude=c.station_elevation,
                                            sat_elev_cutoff=c.station_threshold_elevation)

def is_recorder_alive():
    cmd = "ps ax |grep \"./beacon -d\"|grep -v grep|awk '{print $1}'"
    p = os.popen(cmd)
    lines = p.readlines()
    p.close()
    if len(lines) < 1:
        return(False)
    else:
        return(True)

def do_housekeeping(state):
    fetch_ephemeris(state)
    rec_alive = is_recorder_alive()
    state["rec_alive"]=rec_alive
    p = os.popen("tail -1 %s/rout.log"%(c.datadir))
    last_lines = p.readlines()
    p.close()

    last_an = ""
    if os.path.exists("%s/last_an.txt"%(c.datadir)):
        p = os.popen("cat %s/last_an.txt"%(c.datadir))
        an_lines = p.readlines()
        if len(an_lines) > 0:
            last_an = an_lines[0].strip()
        p.close()

    status_string=""
    if len(last_lines) < 1:
        status_string="not alive"
    else:
        status_string = "%s tle %d recorder alive %d last analyzed %s %s"%(c.station, state["last_tle_fetch"],state["rec_alive"], last_an, last_lines[0].strip())
    print(status_string)
    status_file = file("%s/status.log"%(c.datadir),"w")
    status_file.write(status_string)
    status_file.close()

def run_beacon_calc():
    while True:
        os.system("./beacon_phasecurve.py 2> %s/calc.err > %s/calc.out"%(c.datadir, c.datadir))
        time.sleep(10)

# remove old binary files once a day
def cleanup():
    while True:
        os.system("./beacon_cleanfiles.py > %s/cleanup.out"%(c.datadir))
        time.sleep(60.0*60.0*24.0)

def run_beacon_receiver():
    os.system("killall beacon")
    time.sleep(2)
    cmd = "LC_ALL=\"C\"; ./beacon -d %s -e %d -g %d -G %s -a %s -r \"%s\" -o %s -s %1.20f -t %1.20f -u %d > %s/rout.log 2> %s/rerr.log"%(c.uhd_decimation,
                                                                                                                           c.uhd_refclock,
                                                                                                                           c.uhd_gain0,
                                                                                                                           c.uhd_gain1,
                                                                                                                           c.uhd_device_address,
                                                                                                                           c.uhd_subdev_string,
                                                                                                                           c.datadir,
                                                                                                                           c.uhd_rx_offset0,
                                                                                                                           c.uhd_rx_offset1,
                                                                                                                           0,
                                                                                                                           c.datadir,
                                                                                                                           c.datadir)
    print(cmd)
    os.system(cmd)

if __name__ == "__main__":
    parser = optparse.OptionParser()

    parser.add_option('-x', '--stop_all', dest='stop_all', action='store_true', help='Stop receiver')
    (op, args) = parser.parse_args()

    state = {}
    do_housekeeping(state)
    thread.start_new_thread(run_beacon_receiver, ())
    thread.start_new_thread(run_beacon_calc, ())
    thread.start_new_thread(cleanup, ())
    while True:
        do_housekeeping(state)
        time.sleep(5)


    
