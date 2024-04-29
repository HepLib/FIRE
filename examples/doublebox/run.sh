#!/usr/bin/env bash

export tn=8

echo M/FIRE -c doublebox -t $tn
../../M/FIRE -c doublebox -t $tn
echo

echo M/FIRE -c doublebox2 -t $tn
../../M/FIRE -c doublebox2 -t $tn
echo

echo M/FIRE -c doubleboxl -t $tn
../../M/FIRE -c doubleboxl -t $tn
echo

echo M/FIRE -c doubleboxr -t $tn
../../M/FIRE -c doubleboxr -t $tn
echo

echo M/FIRE -c doubleboxrp -t $tn
../../M/FIRE -c doubleboxrp -t $tn
echo

echo M/FIRE -c doubleboxrp_high -t $tn
../../M/FIRE -c doubleboxrp_high -t $tn
echo
