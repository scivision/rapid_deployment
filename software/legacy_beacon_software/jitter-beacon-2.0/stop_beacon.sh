#!/bin/bash

kill `ps ax |grep beacon.py |grep -v grep |awk '{print $1}' |xargs`