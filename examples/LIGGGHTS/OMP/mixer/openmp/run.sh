#!/bin/bash

##############################################################
# initial run: creates restart file
##############################################################
mpirun -np 1 liggghts \
  -in in_omp.mixer_init \
  -var INSERT_TIMESTEPS 50000 \
  -var TOTAL_TIMESTEPS 100000 \
  -var XPROC 1 \
  -var YPROC 1 \
  -var ZPROC 1 \
  -var IJK xyz \
  -var NTHREADS 8

##############################################################
# benchmark run: simulate 50,000 timesteps
##############################################################

# 1x1x1 @ 1 threads = 1 cores
mpirun -np 1 \
  liggghts \
  -in in_omp.mixer_run \
  -var IJK xyz \
  -var INSERT_TIMESTEPS 50000 \
  -var LB 0 \
  -var LB_DIM x \
  -var NTHREADS 1 \
  -var TOTAL_TIMESTEPS 100000 \
  -var XPROC 1 \
  -var YPROC 1 \
  -var ZPROC 1

# 1x1x1 @ 2 threads = 2 cores
mpirun -np 1 \
  liggghts \
  -in in_omp.mixer_run \
  -var IJK xyz \
  -var INSERT_TIMESTEPS 50000 \
  -var LB 0 \
  -var LB_DIM x \
  -var NTHREADS 2 \
  -var TOTAL_TIMESTEPS 100000 \
  -var XPROC 1 \
  -var YPROC 1 \
  -var ZPROC 1

# 1x1x1 @ 4 threads = 4 cores
mpirun -np 1 \
  liggghts \
  -in in_omp.mixer_run \
  -var IJK xyz \
  -var INSERT_TIMESTEPS 50000 \
  -var LB 0 \
  -var LB_DIM x \
  -var NTHREADS 4 \
  -var TOTAL_TIMESTEPS 100000 \
  -var XPROC 1 \
  -var YPROC 1 \
  -var ZPROC 1

# 1x1x1 @ 8 threads = 8 cores
mpirun -np 1 \
  liggghts \
  -in in_omp.mixer_run \
  -var IJK xyz \
  -var INSERT_TIMESTEPS 50000 \
  -var LB 0 \
  -var LB_DIM x \
  -var NTHREADS 8 \
  -var TOTAL_TIMESTEPS 100000 \
  -var XPROC 1 \
  -var YPROC 1 \
  -var ZPROC 1

# 1x1x1 @ 16 threads = 16 cores
mpirun -np 1 \
  liggghts \
  -in in_omp.mixer_run \
  -var IJK xyz \
  -var INSERT_TIMESTEPS 50000 \
  -var LB 0 \
  -var LB_DIM x \
  -var NTHREADS 16 \
  -var TOTAL_TIMESTEPS 100000 \
  -var XPROC 1 \
  -var YPROC 1 \
  -var ZPROC 1

# 1x1x1 @ 32 threads = 32 cores
mpirun -np 1 \
  liggghts \
  -in in_omp.mixer_run \
  -var IJK xyz \
  -var INSERT_TIMESTEPS 50000 \
  -var LB 0 \
  -var LB_DIM x \
  -var NTHREADS 32 \
  -var TOTAL_TIMESTEPS 100000 \
  -var XPROC 1 \
  -var YPROC 1 \
  -var ZPROC 1
