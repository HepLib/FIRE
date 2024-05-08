#!/usr/bin/env bash

#!/usr/bin/env bash

# due to more variables, we can increase t3 with small lmt3

export tn=4
export t3=16

echo M/FIRE -c nbox2w -t1 $tn -t2 $tn -t3 $t3 -lmt3 10
time ../../M/FIRE -c nbox2w -t1 $tn -t2 $tn -t3 $t3 -lmt3 10
echo
