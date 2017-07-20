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

# 2x1x1 with MPI-LB along x @ 8 threads = 16 cores
mpirun -report-bindings \
  -hostfile myhostfile \
  -np 2 \
  -rf myrankfile \
  liggghts \
  -in in_omp_weights.funnel_run \
  -var LB 1 \
  -var LB_DIM x \
  -var NPARTICLES 1500000 \
  -var PROCX 2 \
  -var PROCY 1 \
  -var PROCZ 1 \
  -var NTHREADS 8 \
  -var TIMESTEPS 100000

# 1x1x2 with MPI-LB along z @ 8 threads = 16 cores
mpirun -report-bindings \
  -hostfile myhostfile \
  -np 2 \
  -rf myrankfile \
  liggghts \
  -in in_omp_weights.funnel_run \
  -var LB 1 \
  -var LB_DIM z \
  -var NPARTICLES 1500000 \
  -var PROCX 1 \
  -var PROCY 1 \
  -var PROCZ 2 \
  -var NTHREADS 8 \
  -var TIMESTEPS 100000

# 2x2x1 with MPI-LB along xy @ 4 threads = 16 cores
mpirun -report-bindings \
  -hostfile myhostfile \
  -np 4 \
  -rf myrankfile \
  liggghts \
  -in in_omp_weights.funnel_run \
  -var LB 1 \
  -var LB_DIM xy \
  -var NPARTICLES 1500000 \
  -var PROCX 2 \
  -var PROCY 2 \
  -var PROCZ 1 \
  -var NTHREADS 4 \
  -var TIMESTEPS 100000

# 2x2x1 with MPI-LB along xy @ 8 threads = 32 cores
mpirun -report-bindings \
  -hostfile myhostfile \
  -np 4 \
  -rf myrankfile \
  liggghts \
  -in in_omp_weights.funnel_run \
  -var LB 1 \
  -var LB_DIM xy \
  -var NPARTICLES 1500000 \
  -var PROCX 2 \
  -var PROCY 2 \
  -var PROCZ 1 \
  -var NTHREADS 8 \
  -var TIMESTEPS 100000

# 2x1x2 with MPI-LB along z @ 8 threads = 32 cores
mpirun -report-bindings \
  -hostfile myhostfile \
  -np 4 \
  -rf myrankfile \
  liggghts \
  -in in_omp_weights.funnel_run \
  -var LB 1 \
  -var LB_DIM z \
  -var NPARTICLES 1500000 \
  -var PROCX 2 \
  -var PROCY 1 \
  -var PROCZ 2 \
  -var NTHREADS 8 \
  -var TIMESTEPS 100000

# 1x1x4 with MPI-LB along z @ 8 threads = 32 cores
mpirun -report-bindings \
  -hostfile myhostfile \
  -np 4 \
  -rf myrankfile \
  liggghts \
  -in in_omp_weights.funnel_run \
  -var LB 1 \
  -var LB_DIM z \
  -var NPARTICLES 1500000 \
  -var PROCX 1 \
  -var PROCY 1 \
  -var PROCZ 4 \
  -var NTHREADS 8 \
  -var TIMESTEPS 100000

# 2x2x2 with MPI-LB along z @ 8 threads = 64 cores
mpirun -report-bindings \
  -hostfile myhostfile128 \
  -np 8 \
  -rf myrankfile \
  liggghts \
  -in in_omp_weights.funnel_run \
  -var LB 1 \
  -var LB_DIM z \
  -var NPARTICLES 1500000 \
  -var PROCX 2 \
  -var PROCY 2 \
  -var PROCZ 2 \
  -var NTHREADS 8 \
  -var TIMESTEPS 100000

# 2x2x4 with MPI-LB along z @ 8 threads = 128 cores
mpirun -report-bindings \
  -hostfile myhostfile128 \
  -np 16 \
  -rf myrankfile \
  liggghts \
  -in in_omp_weights.funnel_run \
  -var LB 1 \
  -var LB_DIM z \
  -var NPARTICLES 1500000 \
  -var PROCX 2 \
  -var PROCY 2 \
  -var PROCZ 4 \
  -var NTHREADS 8 \
  -var TIMESTEPS 100000

# 4x4x1 with MPI-LB along xy @ 8 threads = 128 cores
mpirun -report-bindings \
  -hostfile myhostfile128 \
  -np 16 \
  -rf myrankfile \
  liggghts \
  -in in_omp_weights.funnel_run \
  -var LB 1 \
  -var LB_DIM xy \
  -var NPARTICLES 1500000 \
  -var PROCX 4 \
  -var PROCY 4 \
  -var PROCZ 1 \
  -var NTHREADS 8 \
  -var TIMESTEPS 100000
