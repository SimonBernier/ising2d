#!/bin/bash

set=("$@")
str="Starting run number"
runNumber=0

for i in $(seq 320 1 399);
do
let "runNumber = $(( runNumber + 1 ))"
echo " "
echo "$str $runNumber"

./stTanhRamp $i

done
