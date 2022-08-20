#!/bin/bash

set=("$@")
str="Starting run number"
runNumber=0

for i in $(seq 0 1 26);
do
let "runNumber = $(( runNumber + 1 ))"
echo " "
echo "$str $runNumber"

./super_2dtfi input_2dtfi_Ly_3_runN_$i

done
