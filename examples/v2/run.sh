#!/usr/bin/env bash

export tn=8

echo M/FIRE -c v2 -t $tn
time ../../M/FIRE -c v2 -t $tn
echo

echo M/FIRE -c v2l -t $tn
time ../../M/FIRE -c v2l -t $tn
echo

echo M/FIRE -c v2l_high -t $tn
time ../../M/FIRE -c v2l_high -t $tn
echo
