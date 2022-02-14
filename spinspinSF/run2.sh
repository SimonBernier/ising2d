#!/bin/bash

set=("$@")
str="Starting run number"
runNumber=0
h=0

for i in $(seq 3 1 5);
do
let "runNumber = $(( runNumber + 1 ))"
echo " "
echo "$str $runNumber"

./spinspinSF $i

done
