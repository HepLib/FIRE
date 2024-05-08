#!/usr/bin/env bash

# due to more variables, we can increase t3 with small lmt3

export tn=4

echo M/FIRE -c bananaUnequal -t1 $tn -t2 $tn -t3 $tn -lmt3 10
time ../../M/FIRE -c bananaUnequal -t1 $tn -t2 $tn -t3 $tn -lmt3 10
echo
