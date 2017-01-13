#!/bin/bash
find . -name nanalyzed.log |grep $1 |sed -e 's/.*/rm \0/' |bash