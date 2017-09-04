#!/bin/bash

# 1x1x1 with MPI-LB along z @ 1 threads = 1 cores
mpirun -np 1 \
  liggghts \
  -in in.boxInBox \
  -var PROCX 1 \
  -var PROCY 1 \
  -var PROCZ 1 \
  -var LB 0 \
  -var LB_DIM z

# 2x1x1 with MPI-LB along x @ 1 threads = 2 cores
mpirun -np 2 \
  liggghts \
  -in in.boxInBox \
  -var PROCX 2 \
  -var PROCY 1 \
  -var PROCZ 1 \
  -var LB 1 \
  -var LB_DIM x

# 1x1x2 with MPI-LB along z @ 1 threads = 2 cores
mpirun -np 2 \
  liggghts \
  -in in.boxInBox \
  -var PROCX 1 \
  -var PROCY 1 \
  -var PROCZ 2 \
  -var LB 1 \
  -var LB_DIM z

# 2x2x1 1 threads = 4 cores
mpirun -np 4 \
  liggghts \
  -in in.boxInBox \
  -var PROCX 2 \
  -var PROCY 2 \
  -var PROCZ 1 \
  -var LB 0 \
  -var LB_DIM z

# 2x1x2 with MPI-LB along xz @ 1 threads = 4 cores
mpirun -np 4 \
  liggghts \
  -in in.boxInBox \
  -var PROCX 2 \
  -var PROCY 1 \
  -var PROCZ 2 \
  -var LB 1 \
  -var LB_DIM xz

# 2x2x2 with MPI-LB along z @ 1 threads = 8 cores
mpirun -np 4 \
  liggghts \
  -in in.boxInBox \
  -var PROCX 2 \
  -var PROCY 2 \
  -var PROCZ 2 \
  -var LB 1 \
  -var LB_DIM z

# 4x2x1 with MPI-LB along xy @ 1 threads = 8 cores
mpirun -np 8 \
  liggghts \
  -in in.boxInBox \
  -var PROCX 4 \
  -var PROCY 2 \
  -var PROCZ 1 \
  -var LB 1 \
  -var LB_DIM xy

# 2x2x4 with MPI-LB along z @ 1 threads = 16 cores
mpirun -np 16 \
  liggghts \
  -in in.boxInBox \
  -var PROCX 2 \
  -var PROCY 2 \
  -var PROCZ 4 \
  -var LB 1 \
  -var LB_DIM z

# 4x4x1 with MPI-LB along xy @ 1 threads = 16 cores
mpirun -np 16 \
  liggghts \
  -in in.boxInBox \
  -var PROCX 4 \
  -var PROCY 4 \
  -var PROCZ 1 \
  -var LB 1 \
  -var LB_DIM xy

# 4x4x2 with MPI-LB along z @ 1 threads = 32 cores
mpirun -np 32 \
  liggghts \
  -in in.boxInBox \
  -var PROCX 4 \
  -var PROCY 4 \
  -var PROCZ 2 \
  -var LB 1 \
  -var LB_DIM z

# 4x4x2 with MPI-LB along xyz @ 1 threads = 32 cores
mpirun -np 32 \
  liggghts \
  -in in.boxInBox \
  -var PROCX 4 \
  -var PROCY 4 \
  -var PROCZ 2 \
  -var LB 1 \
  -var LB_DIM xyz

# 4x4x2 @ 1 threads = 32 cores
mpirun -np 32 \
  liggghts \
  -in in.boxInBox \
  -var PROCX 4 \
  -var PROCY 4 \
  -var PROCZ 2 \
  -var LB 0 \
  -var LB_DIM xyz

# 2x2x8 with MPI-LB along xyz @ 1 threads = 32 cores
mpirun -np 32 \
  liggghts \
  -in in.boxInBox \
  -var PROCX 2 \
  -var PROCY 2 \
  -var PROCZ 8 \
  -var LB 1 \
  -var LB_DIM xyz

# 2x2x8 @ 1 threads = 32 cores
mpirun -np 32 \
  liggghts \
  -in in.boxInBox \
  -var PROCX 2 \
  -var PROCY 2 \
  -var PROCZ 8 \
  -var LB 0 \
  -var LB_DIM xyz
