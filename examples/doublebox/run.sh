#!/usr/bin/env bash

export tn=4

echo M/FIRE -c doublebox -t1 $tn
time ../../M/FIRE -c doublebox -t1 $tn
echo

echo M/FIRE -c doublebox2 -t1 $tn
time ../../M/FIRE -c doublebox2 -t1 $tn
echo

echo M/FIRE -c doubleboxl -t1 $tn
time ../../M/FIRE -c doubleboxl -t1 $tn
echo

echo M/FIRE -c doubleboxr -t1 $tn
time ../../M/FIRE -c doubleboxr -t1 $tn
echo

echo M/FIRE -c doubleboxrp -t1 $tn
time ../../M/FIRE -c doubleboxrp -t1 $tn
echo

echo M/FIRE -c doubleboxrp_high -t1 $tn
time ../../M/FIRE -c doubleboxrp_high -t1 $tn
echo
