#!/bin/bash

set=("$@")
str="Starting run number"
runNumber=0

for num in {0..4}
do
let "runNumber = $(( set * 13 + num))"  
echo "$str $runNumber"

./fidelityPT $runNumber

done
