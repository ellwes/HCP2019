#!/bin/bash
export OMP_NUM_THREADS=$1

make 
echo START
aprun -n 1 -d $1 ./bin/dft.out -i dat/input1.dat -v dat/input1_forward.dat -r 10 
echo END
