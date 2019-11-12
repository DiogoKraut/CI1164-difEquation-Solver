#!/bin/bash
rm scripts/L3.dat
for i in {32,50}; do
    sudo likwid-perfctr -C 2 -g L3 -m ./pdeSolver -nx $i -ny $i > out2.dat
    cat out2.dat | grep "L3 bandwidth" | awk '{print $6}' >> scripts/L3.dat
done
