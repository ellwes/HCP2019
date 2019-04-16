#!/bin/bash $1

export OMP_NUM_THREADS=$1
rm -r res1
rm -r res2

make
mkdir res1 
aprun -n 1 -d $1 ./bin/nbody.out  -n 1000 -s 1 -t 10 -d 0.1 -e 0.01 -G 1 -o 1
mv 9.txt res1/9.txt

export OMP_NUM_THREADS=1

mkdir res2 
aprun -n 1 -d 1 ./bin/nbody.out  -n 1000 -s 1 -t 10 -d 0.1 -e 0.01 -G 1 -o 1 
mv 9.txt res2/9.txt

echo OUTPUT:

diff res1/9.txt res2/9.txt

rm *.txt
