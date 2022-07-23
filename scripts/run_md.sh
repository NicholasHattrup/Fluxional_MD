#!/bin/bash

model=$1
xyz=$2
log=$3
samples=$4


# Run a trajectory run 
for (( n=1; n <= $samples; ++n )); do 
	python nequip/nequip/scripts/run_md.py $model $xyz $log --energy-units-to-eV 27.2114 --log-frequency 20 --save-frequency 20 --dt .5 --n-steps 100 --prefix $n
done

