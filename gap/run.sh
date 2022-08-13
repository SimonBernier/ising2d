#!/bin/bash

set=("$@")
str="Starting run number"
runNumber=0

for i in $(seq 1.5 0.01 3);
do
let "runNumber = $(( runNumber + 1 ))"
echo " "
echo "$str $runNumber"

./gap_2dtfi 64 2 $i

done
