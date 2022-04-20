#!/bin/bash

set=("$@")
str="Starting run number"
runNumber=0

for i in $(seq 0.75 0.05 1.0);
do
let "runNumber = $(( runNumber + 1 ))"
echo " "
echo "$str $runNumber"

./vCrit 64 $i

done
