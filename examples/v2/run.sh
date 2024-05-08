#!/usr/bin/env bash

export tn=8

echo M/FIRE -c v2 -t1 $tn
time ../../M/FIRE -c v2 -t1 $tn
echo

echo M/FIRE -c v2l -t1 $tn
time ../../M/FIRE -c v2l -t1 $tn
echo

echo M/FIRE -c v2l_high -t1 $tn
time ../../M/FIRE -c v2l_high -t1 $tn
echo
