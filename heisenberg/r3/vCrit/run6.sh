#!/bin/bash

set=("$@")
str="Starting run number"
runNumber=0

for i in $(seq 16 16 64);
do
let "runNumber = $(( runNumber + 1 ))"
echo " "
echo "$str $runNumber"

./vCrit $i 1.0

done
