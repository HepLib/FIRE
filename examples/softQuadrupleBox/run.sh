#!/usr/bin/env bash

export tn=4

echo M/FIRE -c softQuadrupleBox -t1 $tn -t2 $tn
time ../../M/FIRE -c softQuadrupleBox -t1 $tn -t2 $tn
echo

