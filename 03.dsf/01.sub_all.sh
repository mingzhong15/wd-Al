#!/bin/sh

for((i=1;i<31;i++))
do
	bsub -J zqy_i1.$i < $i.lsf
done

