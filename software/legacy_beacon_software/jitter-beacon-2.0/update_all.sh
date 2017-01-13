#!/bin/bash

ssh tomo@158.39.49.141 "cd /data/beacon ; ./update_beacon.sh ; ./stop_beacon.sh ; ./start_beacon.sh &"
ssh -p 1003 tomo@193.40.1.103 "cd /data0/beacon ; ./update_beacon.sh ; ./stop_beacon.sh ; ./start_beacon.sh& "
ssh tomo@193.167.135.86 "cd /data/beacon ; ./update_beacon.sh ; ./stop_beacon.sh ; ./start_beacon.sh&"
ssh tomo@192.168.90.30 "cd /data/beacon ; ./update_beacon.sh ; ./stop_beacon.sh ; ./start_beacon.sh&"
ssh tomo@193.167.65.123 "cd /data/beacon ; ./update_beacon.sh ; ./stop_beacon.sh ; ./start_beacon.sh&"

