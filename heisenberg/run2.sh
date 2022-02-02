#!/bin/bash

set=("$@")
str="Starting run number"
runNumber=0

for i in $(seq 64 1 127);
do
let "runNumber = $(( runNumber + 1 ))"
echo " "
echo "$str $runNumber"

./stTanhRamp $i

done
