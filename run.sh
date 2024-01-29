#!/bin/bash

exec="./d2q9-bgk.llvm"
n=${1:-'128'}
input="input_${n}x${n}.params"
obst="obstacles_${n}x${n}.dat"
cmd="$exec $input $obst"
echo $cmd
sh -c "$cmd"
