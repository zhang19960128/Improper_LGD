#!/bin/bash
path=$(pwd)
j=0;
index=10;
for i in $( seq -w 0.02 0.02 0.24 )
do
j=$(bc -l <<<"1+$j");
cp $path/../P$i/scf${index} ./scf${j}
cp $path/../P$i/scf${index}.abo ./scf${j}.abo
cp $path/../P$i/itenoe${index} ./itenoe${j}
cp $path/../P$i/itenoe${index}.abo ./itenoe${j}.abo
done
