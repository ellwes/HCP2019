#!/bin/bash
export OMP_NUM_THREADS=$1
aprun -n 1 -d $1 ./stream.out
