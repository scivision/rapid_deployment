#!/bin/sh

while true; do
    ./cmd_all.sh "tail -1 /data*/beacon/data/status.log ; echo ; df -h |grep /data "
    echo "Done"
    sleep 5
done