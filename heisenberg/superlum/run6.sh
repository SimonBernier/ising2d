#!/bin/bash

set=("$@")
str="Starting run number"
runNumber=0

for i in $(seq 208 1 251);
do
let "runNumber = $(( runNumber + 1 ))"
echo " "
echo "$str $runNumber"

./stTanhRamp $i

done
