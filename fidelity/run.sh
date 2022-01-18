#!/bin/bash

set=("$@")
str="Starting run number"
runNumber=0

for num in {0..57}
do
let "runNumber = $(( set * 58 + num))"  
echo " "
echo "$str $runNumber"

./fidelityPT $runNumber

done
