#!/bin/bash
make
aprun -n 1 ./dat/generate.out $1 $1 $1 A.bin B.bin C.answer 1337

aprun -n $2 ./bin/matmul.out A.bin B.bin C_computed.answer

aprun -n 1 ./dat/verify.out C_computed.answer C.answer
