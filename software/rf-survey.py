# Perform a survey over a range of frequncies and record data in DifitalRF format using thor3
# 2017 Gregory Allan

import time
import datetime
import os
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('dir')
parser.add_argument('-r', dest='sample_rate', default='10e6')
parser.add_argument('-c', dest='channels', default='cha,chb')
parser.add_argument('-d', dest='devices', default="\"A:RX1 A:RX2\"")
parser.add_argument('--f-start', dest='start_freq', default='50e6')
parser.add_argument('--f-end', dest='end_freq', default='860e6')
parser.add_argument('-i', dest='interval', default='5e6')
parser.add_argument('-g', dest='gain', default='0')

op = parser.parse_args()

thorcommand = "thor3.py"
#start_freq = 50e6
#end_freq = 860e6
#interval = 4e6
#sample_rate = 10e6
#channels = "cha,chb"
#devices = "\"A:RX1 A:RX2\""

center_freq = eval(op.start_freq)
interval = eval(op.interval)
while center_freq <= eval(op.end_freq):
    now = datetime.datetime.utcnow()
    now = now.replace(microsecond=0)
    starttime = (now + datetime.timedelta(seconds=10)).isoformat() + 'Z'
    endtime = (now + datetime.timedelta(seconds=20)).isoformat() + 'Z'
    bash_command = ' '.join([thorcommand, '-c', op.channels, '-d', op.devices, '-f', str(center_freq), '-g', op.gain, '-r', op.sample_rate, '-s', starttime, '-e', endtime, op.dir])
    print bash_command
    ret = os.system(bash_command)
    if ret:
        raise RuntimeError('thor3.py exited with non-zero status')
    center_freq += interval
    #time.sleep(20)


