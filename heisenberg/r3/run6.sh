#!/bin/bash

set=("$@")
str="Starting run number"
runNumber=0

for i in $(seq 44 1 50);
do
let "runNumber = $(( runNumber + 1 ))"
echo " "
echo "$str $runNumber"

./superlumR3 $i

done
