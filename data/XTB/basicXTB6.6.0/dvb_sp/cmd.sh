#!/bin/sh

[ -f xtbrestart ] && rm xtbrestart
OMP_NUM_THREADS=1 \
    MKL_NUM_THREADS=1 \
    xtb --input dvb_sp.inp --sp dvb_sp.xyz > dvb_sp.out
