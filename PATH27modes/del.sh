#!/bin/bash
path=$(pwd)
for i in $( seq 0 19 )
do
rm -rf T$i/abinit
rm -rf T$i/*DDB
rm -rf T$i/*nc
done
