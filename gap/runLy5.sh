#!/bin/bash

set=("$@")
str="Starting run number"
runNumber=0

for i in $(seq 0.5 0.01 3.5);
do
let "runNumber = $(( runNumber + 1 ))"
echo " "
echo "$str $runNumber"

./gap_2dtfi 40 5 $i

done
