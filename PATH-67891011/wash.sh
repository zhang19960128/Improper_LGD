#!/bin/bash
for i in $(seq 5 19)
do
echo $i
rm T$i/*.nc
rm T$i/*.agr
done
