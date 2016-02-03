#!/bin/bash

##############################################################
# initial run: creates restart file with 1.5 million particles
# uses 8 MPI processes and MPI load-balancing
##############################################################
mpirun -np 8 liggghts \
  -in in.funnel_init \
  -var PROCX 1 \
  -var PROCY 1 \
  -var PROCZ 8 \
  -var NPARTICLES 1500000 \
  -var LB 1

##############################################################
# benchmark run: simulate 100,000 timesteps
##############################################################

# 1x1x1 with MPI-LB along xyz @ 1 threads = 1 cores
mpirun -np 1 \
  liggghts \
  -in in_omp_weights.funnel_run \
  -var LB 0 \
  -var LB_DIM xyz \
  -var NPARTICLES 1500000 \
  -var PROCX 1 \
  -var PROCY 1 \
  -var PROCZ 1 \
  -var NTHREADS 1 \
  -var TIMESTEPS 100000

# 1x1x1 with MPI-LB along xyz @ 2 threads = 2 cores
mpirun -np 1 \
  liggghts \
  -in in_omp_weights.funnel_run \
  -var LB 0 \
  -var LB_DIM xyz \
  -var NPARTICLES 1500000 \
  -var PROCX 1 \
  -var PROCY 1 \
  -var PROCZ 1 \
  -var NTHREADS 2 \
  -var TIMESTEPS 100000

# 1x1x1 with MPI-LB along xyz @ 4 threads = 4 cores
mpirun -np 1 \
  liggghts \
  -in in_omp_weights.funnel_run \
  -var LB 0 \
  -var LB_DIM xyz \
  -var NPARTICLES 1500000 \
  -var PROCX 1 \
  -var PROCY 1 \
  -var PROCZ 1 \
  -var NTHREADS 4 \
  -var TIMESTEPS 100000

# 1x1x1 with MPI-LB along xyz @ 8 threads = 8 cores
mpirun -np 1 \
  liggghts \
  -in in_omp_weights.funnel_run \
  -var LB 0 \
  -var LB_DIM xyz \
  -var NPARTICLES 1500000 \
  -var PROCX 1 \
  -var PROCY 1 \
  -var PROCZ 1 \
  -var NTHREADS 8 \
  -var TIMESTEPS 100000

# 1x1x1 with MPI-LB along xyz @ 16 threads = 16 cores
mpirun -np 1 \
  liggghts \
  -in in_omp_weights.funnel_run \
  -var LB 0 \
  -var LB_DIM xyz \
  -var NPARTICLES 1500000 \
  -var PROCX 1 \
  -var PROCY 1 \
  -var PROCZ 1 \
  -var NTHREADS 16 \
  -var TIMESTEPS 100000

# 1x1x1 with MPI-LB along xyz @ 32 threads = 32 cores
mpirun -np 1 \
  liggghts \
  -in in_omp_weights.funnel_run \
  -var LB 0 \
  -var LB_DIM xyz \
  -var NPARTICLES 1500000 \
  -var PROCX 1 \
  -var PROCY 1 \
  -var PROCZ 1 \
  -var NTHREADS 32 \
  -var TIMESTEPS 100000

