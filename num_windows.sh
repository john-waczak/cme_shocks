#!/bin/bash




for i in `seq 30 300`;
do
    echo $i
    python correctSimulation.py -n $i
    test=$(python chiSquared.py)
    result="$i $test"
    echo $result
    echo $result >> N_windows_test.txt
done
