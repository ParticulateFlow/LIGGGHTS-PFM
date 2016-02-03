#!/bin/bash

# 2x1x1 @ 2 threads = 4 cores
mpirun -np 2 \
  liggghts \
  -in in_omp.boxInBox \
  -var PROCX 2 \
  -var PROCY 1 \
  -var PROCZ 1 \
  -var NTHREADS 2 

# 2x2x1 @ 2 threads = 8 cores
mpirun -np 4 \
  liggghts \
  -in in_omp.boxInBox \
  -var PROCX 2 \
  -var PROCY 2 \
  -var PROCZ 1 \
  -var NTHREADS 2 

# 2x1x1 @ 4 threads = 8 cores
mpirun -np 2 \
  liggghts \
  -in in_omp.boxInBox \
  -var PROCX 2 \
  -var PROCY 1 \
  -var PROCZ 1 \
  -var NTHREADS 4 

# 2x2x1 @ 4 threads = 16 cores
mpirun -np 4 \
  liggghts \
  -in in_omp.boxInBox \
  -var PROCX 2 \
  -var PROCY 2 \
  -var PROCZ 1 \
  -var NTHREADS 4 

# 2x2x1 @ 8 threads = 32 cores
mpirun -np 4 \
  liggghts \
  -in in_omp.boxInBox \
  -var PROCX 2 \
  -var PROCY 2 \
  -var PROCZ 1 \
  -var NTHREADS 8
