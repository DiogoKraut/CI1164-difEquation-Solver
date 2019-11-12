#!/bin/bash
rm scripts/flops.dat
for i in {32,200}; do
    sudo likwid-perfctr -C 2 -g FLOPS_DP -m ./pdeSolver -nx $i -ny $i > out2.dat
    cat out2.dat | grep "DP MFLOPS" | awk '{print $6}' >> scripts/flops.dat
    cat out2.dat | grep "AVX DP"    | awk '{print $6}' >> scripts/avx.dat

done
