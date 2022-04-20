#!/bin/bash

set=("$@")
str="Starting run number"
runNumber=0
h=0

for i in $(seq 1.05 0.05 2.0);
do
let "runNumber = $(( runNumber + 1 ))"
echo " "
echo "$str $runNumber"

./ssc-1dtfi 100 $i

done
