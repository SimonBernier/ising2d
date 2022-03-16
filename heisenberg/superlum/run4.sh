#!/bin/bash

set=("$@")
str="Starting run number"
runNumber=0

for i in $(seq 54 1 71);
do
let "runNumber = $(( runNumber + 1 ))"
echo " "
echo "$str $runNumber"

./stTanhRamp $i

done
