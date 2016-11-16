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

# 2x1x1 with MPI-LB along x @ 8 threads = 16 cores
mpirun -report-bindings \
  -hostfile myhostfile \
  -np 2 \
  -rf myrankfile \
  liggghts \
  -in in_omp.mixer_run \
  -var IJK xyz \
  -var INSERT_TIMESTEPS 50000 \
  -var LB 1 \
  -var LB_DIM x \
  -var NTHREADS 8 \
  -var TOTAL_TIMESTEPS 100000 \
  -var XPROC 2 \
  -var YPROC 1 \
  -var ZPROC 1

# 4x1x1 with MPI-LB along x @ 8 threads = 32 cores
mpirun -report-bindings \
  -hostfile myhostfile \
  -np 4 \
  -rf myrankfile \
  liggghts \
  -in in_omp.mixer_run \
  -var IJK xyz \
  -var INSERT_TIMESTEPS 50000 \
  -var LB 1 \
  -var LB_DIM x \
  -var NTHREADS 8 \
  -var TOTAL_TIMESTEPS 100000 \
  -var XPROC 4 \
  -var YPROC 1 \
  -var ZPROC 1

# 2x2x1 with MPI-LB along x @ 8 threads = 32 cores
mpirun -report-bindings \
  -hostfile myhostfile \
  -np 4 \
  -rf myrankfile \
  liggghts \
  -in in_omp.mixer_run \
  -var IJK xyz \
  -var INSERT_TIMESTEPS 50000 \
  -var LB 1 \
  -var LB_DIM x \
  -var NTHREADS 8 \
  -var TOTAL_TIMESTEPS 100000 \
  -var XPROC 2 \
  -var YPROC 2 \
  -var ZPROC 1

# 8x1x1 with MPI-LB along x @ 8 threads = 64 cores
mpirun -report-bindings \
  -hostfile myhostfile \
  -np 8 \
  -rf myrankfile \
  liggghts \
  -in in_omp.mixer_run \
  -var IJK xyz \
  -var INSERT_TIMESTEPS 50000 \
  -var LB 1 \
  -var LB_DIM x \
  -var NTHREADS 8 \
  -var TOTAL_TIMESTEPS 100000 \
  -var XPROC 8 \
  -var YPROC 1 \
  -var ZPROC 1

# 4x2x1 with MPI-LB along x @ 8 threads = 64 cores
mpirun -report-bindings \
  -hostfile myhostfile \
  -np 8 \
  -rf myrankfile \
  liggghts \
  -in in_omp.mixer_run \
  -var IJK xyz \
  -var INSERT_TIMESTEPS 50000 \
  -var LB 1 \
  -var LB_DIM x \
  -var NTHREADS 8 \
  -var TOTAL_TIMESTEPS 100000 \
  -var XPROC 4 \
  -var YPROC 2 \
  -var ZPROC 1

# 4x1x1 with MPI-LB along x @ 4 threads = 16 cores
mpirun -report-bindings \
  -hostfile myhostfile \
  -np 4 \
  -rf myrankfile \
  liggghts \
  -in in_omp.mixer_run \
  -var IJK xyz \
  -var INSERT_TIMESTEPS 50000 \
  -var LB 1 \
  -var LB_DIM x \
  -var NTHREADS 4 \
  -var TOTAL_TIMESTEPS 100000 \
  -var XPROC 4 \
  -var YPROC 1 \
  -var ZPROC 1

# 1x4x1 with MPI-LB along x @ 8 threads = 32 cores
mpirun -report-bindings \
  -hostfile myhostfile \
  -np 4 \
  -rf myrankfile \
  liggghts \
  -in in_omp.mixer_run \
  -var IJK xyz \
  -var INSERT_TIMESTEPS 50000 \
  -var LB 1 \
  -var LB_DIM y \
  -var NTHREADS 4 \
  -var TOTAL_TIMESTEPS 100000 \
  -var XPROC 1 \
  -var YPROC 4 \
  -var ZPROC 1

# 8x1x1 with MPI-LB along x @ 4 threads = 32 cores
mpirun -report-bindings \
  -hostfile myhostfile \
  -np 4 \
  -rf myrankfile \
  liggghts \
  -in in_omp.mixer_run \
  -var IJK xyz \
  -var INSERT_TIMESTEPS 50000 \
  -var LB 1 \
  -var LB_DIM x \
  -var NTHREADS 4 \
  -var TOTAL_TIMESTEPS 100000 \
  -var XPROC 8 \
  -var YPROC 1 \
  -var ZPROC 1

# 2x1x1 with MPI-LB along x @ 2 threads = 4 cores
mpirun -np 2 \
  liggghts \
  -in in_omp.mixer_run \
  -var IJK xyz \
  -var INSERT_TIMESTEPS 50000 \
  -var LB 1 \
  -var LB_DIM x \
  -var NTHREADS 2 \
  -var TOTAL_TIMESTEPS 100000 \
  -var XPROC 2 \
  -var YPROC 1 \
  -var ZPROC 1

# 2x1x1 with MPI-LB along x @ 2 threads = 8 cores
mpirun -np 2 \
  liggghts \
  -in in_omp.mixer_run \
  -var IJK xyz \
  -var INSERT_TIMESTEPS 50000 \
  -var LB 1 \
  -var LB_DIM x \
  -var NTHREADS 4 \
  -var TOTAL_TIMESTEPS 100000 \
  -var XPROC 2 \
  -var YPROC 1 \
  -var ZPROC 1

# 4x1x1 with MPI-LB along x @ 2 threads = 8 cores
mpirun -np 4 \
  liggghts \
  -in in_omp.mixer_run \
  -var IJK xyz \
  -var INSERT_TIMESTEPS 50000 \
  -var LB 1 \
  -var LB_DIM x \
  -var NTHREADS 2 \
  -var TOTAL_TIMESTEPS 100000 \
  -var XPROC 4 \
  -var YPROC 1 \
  -var ZPROC 1

# 16x1x1 with MPI-LB along x @ 8 threads = 128 cores
mpirun -report-bindings \
  -hostfile myhostfile128 \
  -np 16 \
  -rf myrankfile \
  liggghts \
  -in in_omp.mixer_run \
  -var IJK xyz \
  -var INSERT_TIMESTEPS 50000 \
  -var LB 1 \
  -var LB_DIM xy \
  -var NTHREADS 8 \
  -var TOTAL_TIMESTEPS 100000 \
  -var XPROC 16 \
  -var YPROC 1 \
  -var ZPROC 1

# 8x2x1 with MPI-LB along x @ 8 threads = 128 cores
mpirun -report-bindings \
  -hostfile myhostfile128 \
  -np 16 \
  -rf myrankfile \
  liggghts \
  -in in_omp.mixer_run \
  -var IJK xyz \
  -var INSERT_TIMESTEPS 50000 \
  -var LB 1 \
  -var LB_DIM xy \
  -var NTHREADS 8 \
  -var TOTAL_TIMESTEPS 100000 \
  -var XPROC 8 \
  -var YPROC 2 \
  -var ZPROC 1

# 4x4x1 with MPI-LB along x @ 8 threads = 128 cores
mpirun -report-bindings \
  -hostfile myhostfile128 \
  -np 16 \
  -rf myrankfile \
  liggghts \
  -in in_omp.mixer_run \
  -var IJK xyz \
  -var INSERT_TIMESTEPS 50000 \
  -var LB 1 \
  -var LB_DIM xy \
  -var NTHREADS 8 \
  -var TOTAL_TIMESTEPS 100000 \
  -var XPROC 4 \
  -var YPROC 4 \
  -var ZPROC 1
