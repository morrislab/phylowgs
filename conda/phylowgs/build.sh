#!/bin/bash

mkdir -p $PREFIX/bin

g++ -g -o phylowgs_mh -O3 mh.cpp  util.cpp `gsl-config --cflags --libs`

cp phylowgs_mh $PREFIX/bin/phylowgs_mh

$PYTHON setup.py install

