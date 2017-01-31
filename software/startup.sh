#!/bin/sh

python /home/midasop/midasmicro/scripts/activate_relays.py

sleep 10

sudo mount /dev/sdb1 /data/ssd0
sudo mount /dev/sdc1 /data/ssd1

sudo stty -F /dev/ttyUSB1 speed 115200

sudo systemctl restart gpsd
sudo systemctl restart ntp
