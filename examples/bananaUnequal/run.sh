#!/usr/bin/env bash

# due to more variables, we can increase t3 with small lmt3 or len

t3=8

echo M/FIRE -c bananaUnequal -t3 $t3 -len 100
time ../../M/FIRE -c bananaUnequal -t3 $t3 -len 100
echo
