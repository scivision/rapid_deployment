# Perform a survey over a range of frequncies and record data in DifitalRF format using thor3
# 2017 Gregory Allan

import time
import datetime
import os
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('dir', nargs=1)
parser.add_argument('-r', dest='sample_rate', default=10e6)
parser.add_argument('-c', dest='channels', default='cha,chb')
parser.add_argument('-d', dest='devices', default="\"A:RX1 A:RX2\"")
parser.add_argument('--f-start', dest='start_freq', default=50e6)
parser.add_argument('--f-end', dest='end_freq', default=860e6)
parser.add_argument('-i', dest='interval', default=4e6)

op = parser.parse_args()

thorcommand = "thor3.py"
#start_freq = 50e6
#end_freq = 860e6
#interval = 4e6
#sample_rate = 10e6
#channels = "cha,chb"
#devices = "\"A:RX1 A:RX2\""

center_freq = op.start_freq
while center_freq <= op.end_freq:
    now = datetime.utcnow()
    starttime = (now + timedelta(seconds=5)).isoformat()
    endtime = (now + timedelta(seconds=15)).isoformat()
    bash_command = [thorcommand, '-c', op.channels, '-d', op.devices, '-f', center_freq, '-g', 0, '-r', op.sample_rate, '-s', starttime, '-e', endtime, op.dir].join(' ')
    #os.system(bash_command)
    print bash_command
    center_freq += interval
    time.sleep(20)


