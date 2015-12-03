#!/bin/bash

# 1x1x1 with @ 1 threads = 1 cores
mpirun -np 1 \
  liggghts \
  -in in_omp.boxInBox \
  -var PROCX 1 \
  -var PROCY 1 \
  -var PROCZ 1 \
  -var NTHREADS 1

# 1x1x1 with @ 2 threads = 2 cores
mpirun -np 1 \
  liggghts \
  -in in_omp.boxInBox \
  -var PROCX 1 \
  -var PROCY 1 \
  -var PROCZ 1 \
  -var NTHREADS 2

# 1x1x1 with @ 4 threads = 4 cores
mpirun -np 1 \
  liggghts \
  -in in_omp.boxInBox \
  -var PROCX 1 \
  -var PROCY 1 \
  -var PROCZ 1 \
  -var NTHREADS 4

# 1x1x1 with @ 8 threads = 8 cores
mpirun -np 1 \
  liggghts \
  -in in_omp.boxInBox \
  -var PROCX 1 \
  -var PROCY 1 \
  -var PROCZ 1 \
  -var NTHREADS 8

# 1x1x1 with @ 16 threads = 16 cores
mpirun -np 1 \
  liggghts \
  -in in_omp.boxInBox \
  -var PROCX 1 \
  -var PROCY 1 \
  -var PROCZ 1 \
  -var NTHREADS 16

# 1x1x1 with @ 32 threads = 32 cores
mpirun -np 1 \
  liggghts \
  -in in_omp.boxInBox \
  -var PROCX 1 \
  -var PROCY 1 \
  -var PROCZ 1 \
  -var NTHREADS 32
