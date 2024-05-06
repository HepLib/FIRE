#!/usr/bin/env bash

export tn=8

echo M/FIRE -c doublebox -t $tn
time ../../M/FIRE -c doublebox -t $tn
echo

echo M/FIRE -c doublebox2 -t $tn
time ../../M/FIRE -c doublebox2 -t $tn
echo

echo M/FIRE -c doubleboxl -t $tn
time ../../M/FIRE -c doubleboxl -t $tn
echo

echo M/FIRE -c doubleboxr -t $tn
time ../../M/FIRE -c doubleboxr -t $tn
echo

echo M/FIRE -c doubleboxrp -t $tn
time ../../M/FIRE -c doubleboxrp -t $tn
echo

echo M/FIRE -c doubleboxrp_high -t $tn
time ../../M/FIRE -c doubleboxrp_high -t $tn
echo
