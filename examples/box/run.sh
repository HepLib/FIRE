#!/usr/bin/env bash

export tn=8

echo M/FIRE -c box -t $tn
time ../../M/FIRE -c box -t $tn
echo

echo M/FIRE -c box -t $tn
time ../../M/FIRE -c boxr -t $tn
echo
