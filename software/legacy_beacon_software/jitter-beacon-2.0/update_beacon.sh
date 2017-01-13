#!/bin/bash

rm beacon-2.0*.gz
wget http://www.haystack.mit.edu/~j/beacon-2.0.tar.gz 
tar xvfz beacon-2.0.tar.gz
cp ./beacon-2.0/* .
