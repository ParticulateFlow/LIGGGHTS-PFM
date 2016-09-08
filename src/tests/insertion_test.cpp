#include "gtest/gtest.h"
#include <mpi.h>
#include "atom.h"
#include "input.h"
#include "lammps.h"

using namespace LAMMPS_NS;

TEST(insertion, insertPack_1) {
  const char * argv[3] = {"liggghts", "-in", "scripts/in.insertPack"};
  LAMMPS lammps(3, const_cast<char**>(argv), MPI_COMM_WORLD);
  lammps.input->file();

  lammps.input->one("fix ins all insert/pack seed 100001 distributiontemplate pdd1 vel constant 0. 0. 0. insert_every once overlapcheck yes all_in yes particles_in_region 1 region reg");
  lammps.input->one("run 1");

  EXPECT_EQ(1, lammps.atom->natoms);
}

TEST(insertion, insertPack_150k) {
  const char * argv[3] = {"liggghts", "-in", "scripts/in.insertPack"};
  LAMMPS lammps(3, const_cast<char**>(argv), MPI_COMM_WORLD);
  lammps.input->file();

  lammps.input->one("fix ins all insert/pack seed 100001 distributiontemplate pdd1 vel constant 0. 0. 0. insert_every once overlapcheck yes all_in yes particles_in_region 150000 region reg");
  lammps.input->one("run 1");

  EXPECT_EQ(150000, lammps.atom->natoms);
}
