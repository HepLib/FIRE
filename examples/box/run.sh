#!/usr/bin/env bash

export tn=4

echo M/FIRE -c box -t1 $tn
time ../../M/FIRE -c box -t1 $tn
echo

echo M/FIRE -c box -t1 $tn
time ../../M/FIRE -c boxr -t1 $tn
echo
