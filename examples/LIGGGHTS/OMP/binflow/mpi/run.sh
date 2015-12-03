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
  -in in.funnel_run \
  -var LB 0 \
  -var LB_DIM xyz \
  -var NPARTICLES 1500000 \
  -var PROCX 1 \
  -var PROCY 1 \
  -var PROCZ 1 \
  -var TIMESTEPS 100000

# 1x1x2 with MPI-LB along xyz @ 1 threads = 2 cores
mpirun -np 2 \
  liggghts \
  -in in.funnel_run \
  -var LB 1 \
  -var LB_DIM xyz \
  -var NPARTICLES 1500000 \
  -var PROCX 1 \
  -var PROCY 1 \
  -var PROCZ 2 \
  -var TIMESTEPS 100000

# 2x1x1 with MPI-LB along x @ 1 threads = 2 cores
mpirun -np 2 \
  liggghts \
  -in in.funnel_run \
  -var LB 1 \
  -var LB_DIM x \
  -var NPARTICLES 1500000 \
  -var PROCX 2 \
  -var PROCY 1 \
  -var PROCZ 1 \
  -var TIMESTEPS 100000

# 1x2x1 with MPI-LB along y @ 1 threads = 2 cores
mpirun -np 2 \
  liggghts \
  -in in.funnel_run \
  -var LB 1 \
  -var LB_DIM y \
  -var NPARTICLES 1500000 \
  -var PROCX 1 \
  -var PROCY 2 \
  -var PROCZ 1 \
  -var TIMESTEPS 100000

# 1x1x4 with MPI-LB along z @ 1 threads = 4 cores
mpirun -np 4 \
  liggghts \
  -in in.funnel_run \
  -var LB 1 \
  -var LB_DIM z \
  -var NPARTICLES 1500000 \
  -var PROCX 1 \
  -var PROCY 1 \
  -var PROCZ 4 \
  -var TIMESTEPS 100000

# 2x1x2 with MPI-LB along xz @ 1 threads = 4 cores
mpirun -np 4 \
  liggghts \
  -in in.funnel_run \
  -var LB 1 \
  -var LB_DIM xz \
  -var NPARTICLES 1500000 \
  -var PROCX 2 \
  -var PROCY 1 \
  -var PROCZ 2 \
  -var TIMESTEPS 100000

# 4x1x1 with MPI-LB along x @ 1 threads = 4 cores
mpirun -np 4 \
  liggghts \
  -in in.funnel_run \
  -var LB 1 \
  -var LB_DIM x \
  -var NPARTICLES 1500000 \
  -var PROCX 4 \
  -var PROCY 1 \
  -var PROCZ 1 \
  -var TIMESTEPS 100000

# 4x1x1 with MPI-LB along x @ 1 threads = 4 cores
mpirun -np 4 \
  liggghts \
  -in in.funnel_run \
  -var LB 1 \
  -var LB_DIM x \
  -var NPARTICLES 1500000 \
  -var PROCX 4 \
  -var PROCY 1 \
  -var PROCZ 1 \
  -var TIMESTEPS 100000

# 1x1x8 with MPI-LB along z @ 1 threads = 8 cores
mpirun -np 8 \
  liggghts \
  -in in.funnel_run \
  -var LB 1 \
  -var LB_DIM z \
  -var NPARTICLES 1500000 \
  -var PROCX 1 \
  -var PROCY 1 \
  -var PROCZ 8 \
  -var TIMESTEPS 100000

# 2x2x2 with MPI-LB along z @ 1 threads = 8 cores
mpirun -np 8 \
  liggghts \
  -in in.funnel_run \
  -var LB 1 \
  -var LB_DIM z \
  -var NPARTICLES 1500000 \
  -var PROCX 2 \
  -var PROCY 2 \
  -var PROCZ 2 \
  -var TIMESTEPS 100000

# 2x1x4 with MPI-LB along xz @ 1 threads = 8 cores
mpirun -np 8 \
  liggghts \
  -in in.funnel_run \
  -var LB 1 \
  -var LB_DIM xz \
  -var NPARTICLES 1500000 \
  -var PROCX 2 \
  -var PROCY 1 \
  -var PROCZ 4 \
  -var TIMESTEPS 100000

# 4x1x2 with MPI-LB along xz @ 1 threads = 8 cores
mpirun -np 8 \
  liggghts \
  -in in.funnel_run \
  -var LB 1 \
  -var LB_DIM xz \
  -var NPARTICLES 1500000 \
  -var PROCX 4 \
  -var PROCY 1 \
  -var PROCZ 2 \
  -var TIMESTEPS 100000

# 4x2x1 with MPI-LB along xy @ 1 threads = 8 cores
mpirun -np 8 \
  liggghts \
  -in in.funnel_run \
  -var LB 1 \
  -var LB_DIM xy \
  -var NPARTICLES 1500000 \
  -var PROCX 4 \
  -var PROCY 2 \
  -var PROCZ 1 \
  -var TIMESTEPS 100000

# 2x2x4 with MPI-LB along z @ 1 threads = 16 cores
mpirun -np 16 \
  liggghts \
  -in in.funnel_run \
  -var LB 1 \
  -var LB_DIM z \
  -var NPARTICLES 1500000 \
  -var PROCX 2 \
  -var PROCY 2 \
  -var PROCZ 4 \
  -var TIMESTEPS 100000
  
# 4x2x2 with MPI-LB along xz @ 1 threads = 16 cores
mpirun -np 16 \
  liggghts \
  -in in.funnel_run \
  -var LB 1 \
  -var LB_DIM xz \
  -var NPARTICLES 1500000 \
  -var PROCX 4 \
  -var PROCY 2 \
  -var PROCZ 2 \
  -var TIMESTEPS 100000

# 2x1x8 with MPI-LB along xz @ 1 threads = 16 cores
mpirun -np 16 \
  liggghts \
  -in in.funnel_run \
  -var LB 1 \
  -var LB_DIM xz \
  -var NPARTICLES 1500000 \
  -var PROCX 2 \
  -var PROCY 1 \
  -var PROCZ 8 \
  -var TIMESTEPS 100000

# 8x1x2 with MPI-LB along xz @ 1 threads = 16 cores
mpirun -np 16 \
  liggghts \
  -in in.funnel_run \
  -var LB 1 \
  -var LB_DIM xz \
  -var NPARTICLES 1500000 \
  -var PROCX 8 \
  -var PROCY 1 \
  -var PROCZ 2 \
  -var TIMESTEPS 100000

# 4x4x2 with MPI-LB along xyz @ 1 threads = 32 cores
mpirun -np 32 \
  liggghts \
  -in in.funnel_run \
  -var LB 1 \
  -var LB_DIM xyz \
  -var NPARTICLES 1500000 \
  -var PROCX 4 \
  -var PROCY 4 \
  -var PROCZ 2 \
  -var TIMESTEPS 100000

# 2x2x8 with MPI-LB along z @ 1 threads = 32 cores
mpirun -np 32 \
  liggghts \
  -in in.funnel_run \
  -var LB 1 \
  -var LB_DIM z \
  -var NPARTICLES 1500000 \
  -var PROCX 2 \
  -var PROCY 2 \
  -var PROCZ 8 \
  -var TIMESTEPS 100000

# 4x2x4 with MPI-LB along xz @ 1 threads = 32 cores
mpirun -np 32 \
  liggghts \
  -in in.funnel_run \
  -var LB 1 \
  -var LB_DIM xz \
  -var NPARTICLES 1500000 \
  -var PROCX 4 \
  -var PROCY 2 \
  -var PROCZ 4 \
  -var TIMESTEPS 100000

# 4x2x8 with MPI-LB along z @ 1 threads = 64 cores
mpirun -report-bindings \
  -hostfile myhostfile64 \
  -np 64 \
  liggghts \
  -in in.funnel_run \
  -var LB 1 \
  -var LB_DIM z \
  -var NPARTICLES 1500000 \
  -var PROCX 4 \
  -var PROCY 2 \
  -var PROCZ 8 \
  -var TIMESTEPS 100000

# 4x2x8 with MPI-LB along xyz @ 1 threads = 64 cores
mpirun -report-bindings \
  -hostfile myhostfile64 \
  -np 64 \
  liggghts \
  -in in.funnel_run \
  -var LB 1 \
  -var LB_DIM xyz \
  -var NPARTICLES 1500000 \
  -var PROCX 4 \
  -var PROCY 2 \
  -var PROCZ 8 \
  -var TIMESTEPS 100000

# 4x2x8 with MPI-LB along xz @ 1 threads = 64 cores
mpirun -report-bindings \
  -hostfile myhostfile64 \
  -np 64 \
  liggghts \
  -in in.funnel_run \
  -var LB 1 \
  -var LB_DIM xz \
  -var NPARTICLES 1500000 \
  -var PROCX 4 \
  -var PROCY 2 \
  -var PROCZ 8 \
  -var TIMESTEPS 100000

# 4x4x4 with MPI-LB along xyz @ 1 threads = 64 cores
mpirun -report-bindings \
  -hostfile myhostfile64 \
  -np 64 \
  liggghts \
  -in in.funnel_run \
  -var LB 1 \
  -var LB_DIM xyz \
  -var NPARTICLES 1500000 \
  -var PROCX 4 \
  -var PROCY 4 \
  -var PROCZ 4 \
  -var TIMESTEPS 100000

# 4x4x8 with MPI-LB along xyz @ 1 threads = 128 cores
mpirun -report-bindings \
  -hostfile myhostfile128 \
  -np 128 \
  liggghts \
  -in in.funnel_run \
  -var LB 1 \
  -var LB_DIM xyz \
  -var NPARTICLES 1500000 \
  -var PROCX 4 \
  -var PROCY 4 \
  -var PROCZ 8 \
  -var TIMESTEPS 100000