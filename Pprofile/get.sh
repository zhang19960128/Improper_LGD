#!/bin/bash
for i in $(seq 0 19 )
do
grep etotal2 scf${i}.abo
done
