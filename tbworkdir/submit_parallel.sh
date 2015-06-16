#/bin/bash

for i in 7 8 11 12 13 16 17 18 22 24 
#for i in 13
do
  qsub -N HMM_$i -o yeah$i.log -vtagid=$i runscript.sh
done
