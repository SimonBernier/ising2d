#!/bin/bash

set=("$@")
str="Starting run number"
runNumber=0
h=0

for i in $(seq 0.05 0.1 0.95);
do
let "runNumber = $(( runNumber + 1 ))"
echo " "
echo "$str $runNumber"

./ssc-1dtfi 100 $i 100 0.1 1

done
