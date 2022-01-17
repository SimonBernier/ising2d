#!/bin/bash

set=("$@")
str="Starting run number"
runNumber=0

for num in {0..12}
do
let "runNumber = $(( set * 26 + 13 + num))"  
echo " "
echo "$str $runNumber"

./fidelityPT $runNumber

done
