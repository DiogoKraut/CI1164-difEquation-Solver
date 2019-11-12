#!/bin/bash
rm scripts/time.dat scripts/L2.dat
for i in {32,50}; do
    sudo likwid-perfctr -C 2 -g L2CACHE -m ./pdeSolver -nx $i -ny $i > out2.dat
    cat out2.dat | grep "RDTSC Runtime" | awk '{print $6}' >> scripts/time.dat
    cat out2.dat | grep "miss ratio"    | awk '{print $6}' >> scripts/L2.dat
done
