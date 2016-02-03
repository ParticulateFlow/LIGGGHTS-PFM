#!/bin/bash

##############################################################
# initial run: creates restart file
##############################################################
mpirun -np 8 liggghts \
  -in in.mixer_init \
  -var INSERT_TIMESTEPS 50000 \
  -var TOTAL_TIMESTEPS 100000 \
  -var XPROC 8 \
  -var YPROC 1 \
  -var ZPROC 1 \
  -var IJK xyz

##############################################################
# benchmark run: simulate 50,000 timesteps
##############################################################

# 1x1x1 with MPI-LB along xyz @ 1 threads = 1 cores
mpirun -np 1 \
  liggghts \
  -in in.mixer_run \
  -var IJK xyz \
  -var INSERT_TIMESTEPS 50000 \
  -var LB 0 \
  -var LB_DIM xy \
  -var TOTAL_TIMESTEPS 100000 \
  -var XPROC 1 \
  -var YPROC 1 \
  -var ZPROC 1

# 2x1x1 with MPI-LB along x @ 1 threads = 2 cores
mpirun -np 2 \
  liggghts \
  -in in.mixer_run \
  -var IJK xyz \
  -var INSERT_TIMESTEPS 50000 \
  -var LB 1 \
  -var LB_DIM x \
  -var TOTAL_TIMESTEPS 100000 \
  -var XPROC 2 \
  -var YPROC 1 \
  -var ZPROC 1

# 4x1x1 with MPI-LB along x @ 1 threads = 4 cores
mpirun -np 4 \
  liggghts \
  -in in.mixer_run \
  -var IJK xyz \
  -var INSERT_TIMESTEPS 50000 \
  -var LB 1 \
  -var LB_DIM x \
  -var TOTAL_TIMESTEPS 100000 \
  -var XPROC 4 \
  -var YPROC 1 \
  -var ZPROC 1

# 8x1x1 with MPI-LB along x @ 1 threads = 8 cores
mpirun -np 8 \
  liggghts \
  -in in.mixer_run \
  -var IJK xyz \
  -var INSERT_TIMESTEPS 50000 \
  -var LB 1 \
  -var LB_DIM x \
  -var TOTAL_TIMESTEPS 100000 \
  -var XPROC 8 \
  -var YPROC 1 \
  -var ZPROC 1

# 8x2x1 with MPI-LB along xy @ 1 threads = 16 cores
mpirun -np 16 \
  liggghts \
  -in in.mixer_run \
  -var IJK xyz \
  -var INSERT_TIMESTEPS 50000 \
  -var LB 1 \
  -var LB_DIM xy \
  -var TOTAL_TIMESTEPS 100000 \
  -var XPROC 8 \
  -var YPROC 2 \
  -var ZPROC 1

# 8x4x1 with MPI-LB along xy @ 1 threads = 32 cores
mpirun -np 32 \
  liggghts \
  -in in.mixer_run \
  -var IJK xyz \
  -var INSERT_TIMESTEPS 50000 \
  -var LB 1 \
  -var LB_DIM xy \
  -var TOTAL_TIMESTEPS 100000 \
  -var XPROC 8 \
  -var YPROC 4 \
  -var ZPROC 1

# 4x16x1 with MPI-LB along xy @ 1 threads = 64 cores
mpirun -report-bindings \
  -hostfile myhostfile64 \
  -np 64 \
  liggghts \
  -in in.mixer_run \
  -var IJK xyz \
  -var INSERT_TIMESTEPS 50000 \
  -var LB 1 \
  -var LB_DIM xy \
  -var TOTAL_TIMESTEPS 100000 \
  -var XPROC 4 \
  -var YPROC 16 \
  -var ZPROC 1

# 8x2x4 with MPI-LB along xz @ 1 threads = 64 cores
mpirun -report-bindings \
  -hostfile myhostfile64 \
  -np 64 \
  liggghts \
  -in in.mixer_run \
  -var IJK xyz \
  -var INSERT_TIMESTEPS 50000 \
  -var LB 1 \
  -var LB_DIM xz \
  -var TOTAL_TIMESTEPS 100000 \
  -var XPROC 8 \
  -var YPROC 2 \
  -var ZPROC 4

# 8x8x1 with MPI-LB along xz @ 1 threads = 64 cores
mpirun -report-bindings \
  -hostfile myhostfile64 \
  -np 64 \
  liggghts \
  -in in.mixer_run \
  -var IJK xyz \
  -var INSERT_TIMESTEPS 50000 \
  -var LB 1 \
  -var LB_DIM xy \
  -var TOTAL_TIMESTEPS 100000 \
  -var XPROC 8 \
  -var YPROC 8 \
  -var ZPROC 1

# 16x4x1 with MPI-LB along xy @ 1 threads = 64 cores
mpirun -report-bindings \
  -hostfile myhostfile64 \
  -np 64 \
  liggghts \
  -in in.mixer_run \
  -var IJK xyz \
  -var INSERT_TIMESTEPS 50000 \
  -var LB 1 \
  -var LB_DIM xy \
  -var TOTAL_TIMESTEPS 100000 \
  -var XPROC 16 \
  -var YPROC 4 \
  -var ZPROC 1

# 8x4x2 with MPI-LB along xyz @ 1 threads = 64 cores
mpirun -report-bindings \
  -hostfile myhostfile64 \
  -np 64 \
  liggghts \
  -in in.mixer_run \
  -var IJK xyz \
  -var INSERT_TIMESTEPS 50000 \
  -var LB 1 \
  -var LB_DIM xyz \
  -var TOTAL_TIMESTEPS 100000 \
  -var XPROC 8 \
  -var YPROC 4 \
  -var ZPROC 2

# 16x8x1 with MPI-LB along xy @ 1 threads = 128 cores
mpirun -report-bindings \
  -hostfile myhostfile128 \
  -np 128 \
  liggghts \
  -in in.mixer_run \
  -var IJK xyz \
  -var INSERT_TIMESTEPS 50000 \
  -var LB 1 \
  -var LB_DIM xy \
  -var TOTAL_TIMESTEPS 100000 \
  -var XPROC 16 \
  -var YPROC 8 \
  -var ZPROC 1
