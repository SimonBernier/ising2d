#!/bin/bash

set=("$@")
str="Starting run number"
runNumber=0

for i in $(seq 12 1 14);
do
let "runNumber = $(( runNumber + 1 ))"
echo " "
echo "$str $runNumber"

./spinspinSF $i

done
