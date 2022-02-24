#!/bin/bash

set=("$@")
str="Starting run number"
runNumber=0

for i in $(seq 160 1 239);
do
let "runNumber = $(( runNumber + 1 ))"
echo " "
echo "$str $runNumber"

./hyperbolicTanh $i

done
