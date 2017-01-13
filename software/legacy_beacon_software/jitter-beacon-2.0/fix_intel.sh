#!/bin/bash
sudo rmmod e1000e
sudo modprobe e1000e InterruptThrottleRate=0,0,0,0
