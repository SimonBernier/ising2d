#!/bin/bash

set=("$@")
str="Starting run number"
runNumber=0
h=0

for i in $(seq 0.1 0.1 4.0);
do
let "runNumber = $(( runNumber + 1 ))"
echo " "
echo "$str $runNumber"

./ssc-1dtfi 3 32 $i 100 0.1 2

done
