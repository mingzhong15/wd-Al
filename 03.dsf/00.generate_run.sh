#!/bin/sh

for((i=1;i<41;i++))
do
	sed "s/python 1/python $i/g" run.lsf > $i.lsf
	sed "s/NUM = 1/NUM = $i/g" test.py > $i.py
done

