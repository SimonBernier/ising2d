#!/bin/bash

set=("$@")
str="Starting run number"
runNumber=0

for i in $(seq 0 1 79);
do
let "runNumber = $(( runNumber + 1 ))"
echo " "
echo "$str $runNumber"

./hyperbolicTanh $i

done
