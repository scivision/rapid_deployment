#!/usr/bin/env python
import os, datetime, re, time, calendar, sys
import optparse

import smtplib
from email.mime.text import MIMEText

#{"name":"SVA","cmd":"ssh tomo@158.39.49.141"},
stations = [{"name":"SVA","cmd":"ssh -o ConnectTimeout=5 tomo@158.39.49.141","scpcmd":"scp","scpdest":"tomo@158.39.49.141:/data/beacon"},
            {"name":"TAR","cmd":"ssh -o ConnectTimeout=5 -p 1003 tomo@193.40.1.103","scpcmd":"scp -P 1003","scpdest":"tomo@193.40.1.103:/data0/beacon"},
            {"name":"SOD","cmd":"ssh -o ConnectTimeout=5 tomo@193.167.135.86","scpcmd":"scp","scpdest":"tomo@193.167.135.86:/data/beacon"},
            {"name":"KEV","cmd":"ssh -o ConnectTimeout=5 tomo@192.168.90.30","scpcmd":"scp","scpdest":"tomo@192.168.90.30:/data/beacon"},
            {"name":"MEK","cmd":"ssh -o ConnectTimeout=5 tomo@193.167.65.123","scpcmd":"scp","scpdest":"tomo@193.167.65.123:/data/beacon"},
            {"name":"TRO","cmd":"ssh -o ConnectTimeout=5 tomo@129.242.30.160","scpcmd":"scp","scpdest":"tomo@129.242.30.160:/data/beacon"}]

warn_when_analysis_late = True
logfile = "/home/users/jvierine/tomo/beacon/monitor.log"
admin_email = ["jvierine@gmail.com","jussi.norberg@gmail.com"]

enable_email=True

def cmd_station(s,cmd):
    p = os.popen("%s \"%s\""%(s["cmd"],cmd))
    l = p.readlines()
    if len(l) > 0:
        l = l[0].strip()
    p.close()
    return(l)

def get_beacon_configuration(s):
    print(s["name"])
    os.system("%s \"cat /data*/beacon/beacon_conf.py\" > beacon_conf.py.%s"%(s["cmd"],s["name"]))

def send_beacon_configuration(s):
    print(s["name"])
    os.system("%s beacon_conf.py.%s %s/beacon_conf.py"%(s["scpcmd"],s["name"],s["scpdest"]))

def restart_beacon(s):
    cmd_station(s,"cd /data*/beacon ; ./stop_beacon.sh ; sleep 5 ; ./start_beacon.sh")

def get_beacon_pid(s):
    pids = cmd_station(s,"ps ax | grep 'python beacon.py'| grep -v grep")
    pid=-1
    if len(pids) > 0:
        pid = int(pids.split(" ")[0])
    return(pid)

def get_beacon_space(s):
    space = cmd_station(s,"df 2>/dev/null|grep /data")
    free = 0.0
    if len(space) > 0:
        space = space.strip()
        space = space.split()
        free = float(space[3])/float(space[1])
    return(free)

def email_error(error_msg,me="jvierine@plasma.mh.fmi.fi"):
    if enable_email:
        msg = MIMEText(error_msg)
        for rec in admin_email:
            # you == the recipient's email address
            msg['Subject'] = 'Beacon error %s'%(datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))
            msg['From'] = me
            msg['To'] = rec

            s = smtplib.SMTP('localhost')
            s.sendmail(me, [rec], msg.as_string())
            s.quit()

def get_beacon_status(s):
    stat = cmd_station(s,"tail -1 /data*/beacon/data/status.log")
    if len(stat) > 0:
        try:
            tlast = calendar.timegm(time.strptime(re.search("/data.*/beacon/data/(.*)/.*",stat).group(1),"%Y.%m.%d"))
            dt = time.time()-tlast
            if warn_when_analysis_late:
                if dt > 2.0*60.0*60.0*24.0:
                    email_error("Latest analyzed file more than 2 days from now in %s\n%s"%(s["name"],stat))
        except:
            email_error("Error reading last analyzed file @ %s."%(s["name"]))
            
    return(stat)

parser = optparse.OptionParser()
parser.add_option('-d', '--display', dest='display', action='store_true',
                  help='Display output')
parser.add_option('-g', '--get_config', dest='get_config', action='store_true',
                  help='Get configuration')
parser.add_option('-s', '--send_config', dest='send_config', action='store_true',
                  help='Send configuration')
parser.add_option('-r', '--restart_receivers', dest='restart_receivers', action='store_true',
                  help='Restart receivers')
parser.add_option('-R', '--restart_receiver', dest='restart_receiver', action='store',type='string',
                  help='Restart receiver',default=None)

(op, args) = parser.parse_args()
enable_email=not op.display

if op.get_config:
    for s in stations:
        get_beacon_configuration(s)
    exit(0)

if op.send_config:
    for s in stations:
        send_beacon_configuration(s)
    exit(0)

if op.restart_receivers:
    for s in stations:
        print("restarting %s"%(s["name"]))
        restart_beacon(s)
    exit(0)
if op.restart_receiver != None:
    for s in stations:
        if s["name"] == op.restart_receiver:
            print("restarting %s"%(s["name"]))
            restart_beacon(s)
    exit(0)
    

f=None
if op.display:
    f= sys.stdout
else:
    f = file(logfile,"w")        

for s in stations:
    pid = get_beacon_pid(s)
    if pid == -1:
        print("Error, no receiver running, or no connection to host.")
        email_error("Beacon receiver not running on %s, or no connection to host."%(s["name"]))
    else:
        free = get_beacon_space(s)
        stat = get_beacon_status(s)

        if free == None or free < 0.1:
            email_error("Less than 10%% disk space left on %s"%(s["name"]))

        f.write("%s pid %d free disk space %1.2f %% %s\n"%(s["name"],pid,free*100.0,stat))
f.close()

