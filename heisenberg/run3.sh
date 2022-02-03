#!/bin/bash

set=("$@")
str="Starting run number"
runNumber=0

for i in $(seq 1728 1 1791);
do
let "runNumber = $(( runNumber + 1 ))"
echo " "
echo "$str $runNumber"

./stTanhRamp $i

done
