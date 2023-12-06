#!/bin/bash
j=12;
path=$(pwd)
for i in $(seq -w -0.30 0.05 -0.05)
do
j=$(bc -l <<<"1+$j")
cp $path/../P$i/scf10 ./scf$j
cp $path/../P$i/scf10.abo ./scf${j}.abo
cp $path/../P$i/itenoe10 ./itenoe$j
cp $path/../P$i/itenoe10.abo ./itenoe${j}.abo
done
j=$(bc -l <<<"1+$j");
cp $path/../P-0.05/scf0 ./scf$j
cp $path/../P-0.05/scf0.abo ./scf${j}.abo
cp $path/../P-0.05/itenoe0 ./itenoe$j
cp $path/../P-0.05/itenoe0.abo ./itenoe${j}.abo

