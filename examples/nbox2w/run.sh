#!/usr/bin/env bash

# due to more variables, we can increase -tp with small -lmt/-len

tp=16 # thread pool size

echo M/FIRE -c nbox2w -tp $tp -len 100
time ../../M/FIRE -c nbox2w -tp $tp -len 100
echo
