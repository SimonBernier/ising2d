#!/bin/bash

set=("$@")
str="Starting run number"
runNumber=0

for i in $(seq 32 16 192);
do
let "runNumber = $(( runNumber + 1 ))"
echo " "
echo "$str $runNumber"

./heisCrit $i

done
