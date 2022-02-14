#!/bin/bash

set=("$@")
str="Starting run number"
runNumber=0
h=0

for i in $(seq 6 1 8);
do
let "runNumber = $(( runNumber + 1 ))"
echo " "
echo "$str $runNumber"

./spinspinSF $i

done
