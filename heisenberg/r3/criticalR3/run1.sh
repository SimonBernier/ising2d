#!/bin/bash

set=("$@")
str="Starting run number"
runNumber=0

for i in $(seq 64 16 192);
do
let "runNumber = $(( runNumber + 1 ))"
echo " "
echo "$str $runNumber"

./criticalR3 $i 0.1

done
