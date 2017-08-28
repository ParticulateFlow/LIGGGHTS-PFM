#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include <mpi.h>
#include "atom.h"
#include "input.h"
#include "lammps.h"
#include <cmath>
#include <vector>
#include <functional>
#include <algorithm>
#include <iostream>
#include <random>
#include "random_park.h"

using namespace LAMMPS_NS;
using namespace testing;

TEST(particle_distribution_discrete, uniform) {
  const int nrolls=10000;  // number of experiments
  const int nstars=95;     // maximum number of stars to distribute
  const int nintervals=10; // number of intervals

  std::default_random_engine generator;
  std::uniform_real_distribution<double> distribution(0.0,1.0);

  std::vector<int> p(nintervals, 0);
  std::vector<int> starCount(nintervals, 0);

  for (int i=0; i<nrolls; ++i) {
    double number = distribution(generator);
    ++p[int(nintervals*number)];
  }

  std::cout << "uniform_real_distribution (0.0,1.0):" << std::endl;
  std::cout << std::fixed; std::cout.precision(1);

  for (int i=0; i<nintervals; ++i) {
    std::cout << float(i)/nintervals << "-" << float(i+1)/nintervals << ": ";
    std::cout << std::string(p[i]*nstars/nrolls,'*') << std::endl;
  }
}


TEST(particle_distribution_discrete, uniform_park) {
  const int nrolls=10000;  // number of experiments
  const int nstars=95;     // maximum number of stars to distribute
  const int nintervals=10; // number of intervals

  const char * argv[3] = {"liggghts", "-in", "scripts/in.insertPackCustom"};
  LAMMPS lammps(3, const_cast<char**>(argv), MPI_COMM_WORLD);
  RanPark r(&lammps, 1);

  std::vector<int> p(nintervals, 0);
  std::vector<int> starCount(nintervals, 0);

  for (int i=0; i<nrolls; ++i) {
    double number = r.uniform();
    ++p[int(nintervals*number)];
  }

  std::cout << "RanPark (0.0,1.0):" << std::endl;
  std::cout << std::fixed; std::cout.precision(1);

  for (int i=0; i<nintervals; ++i) {
    std::cout << float(i)/nintervals << "-" << float(i+1)/nintervals << ": ";
    std::cout << std::string(p[i]*nstars/nrolls,'*') << std::endl;
    starCount[i] = static_cast<int>(p[i]*nstars/nrolls);
  }

  ASSERT_THAT(starCount, Each(starCount[0]));
}


/*TEST(particle_distribution_discrete, tube) {
  int nprocs = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  std::vector<const char*> args;
  args.push_back("liggghts");
  args.push_back("-in");
  args.push_back("scripts/in.tube");
  args.push_back("-var");
  args.push_back("DCYLDP");
  args.push_back("10");

  switch(nprocs) {
    case 1:
      args.push_back("-var"); args.push_back("XPROCS"); args.push_back("1");
      args.push_back("-var"); args.push_back("YPROCS"); args.push_back("1");
      args.push_back("-var"); args.push_back("ZPROCS"); args.push_back("1");
      break;

    case 32:
      args.push_back("-var"); args.push_back("XPROCS"); args.push_back("4");
      args.push_back("-var"); args.push_back("YPROCS"); args.push_back("4");
      args.push_back("-var"); args.push_back("ZPROCS"); args.push_back("2");
      break;

    default:
      FAIL() << "unsupported nprocs";
  }

  LAMMPS lammps(args.size(), const_cast<char**>(&args[0]), MPI_COMM_WORLD);
  lammps.input->file();

  Atom * atom = lammps.atom;
  std::vector<int> count(6, 0);
  std::vector<double> radii;
  radii.push_back(0.001);
  radii.push_back(0.0008);
  radii.push_back(0.000625);
  radii.push_back(0.0005);
  radii.push_back(0.000315);
  radii.push_back(0.0001575);

  for(int i = 0; i < atom->nlocal; ++i) {
    for(size_t j = 0; j < radii.size(); ++j) {
      if(std::abs(atom->radius[i] - radii[j]) < 0.00001) {
        ++count[j];
      }
    }
  }

  int total = 0;
  for(size_t j = 0; j < count.size(); ++j) total += count[j];

  EXPECT_EQ(atom->nlocal, total);

  for(size_t j = 0; j < count.size(); ++j) {
    EXPECT_TRUE(count[j] > 0) << "radii[" << j << "] = " << radii[j] << " does not appear";
  }
}*/


TEST(particle_distribution_discrete, two_types) {
  const char * argv[3] = {"liggghts", "-in", "scripts/in.insertPackCustom"};
  LAMMPS lammps(3, const_cast<char**>(argv), MPI_COMM_WORLD);
  lammps.input->file();

  lammps.input->one("fix pts1 all particletemplate/sphere 1 atom_type 1 density constant 2500 radius constant 0.0025");
  lammps.input->one("fix pts2 all particletemplate/sphere 1 atom_type 1 density constant 2500 radius constant 0.0050");
  lammps.input->one("fix pdd1 all particledistribution/discrete 1.  2 pts1 0.5 pts2 0.5");

  lammps.input->one("fix ins all insert/pack seed 100001 distributiontemplate pdd1 vel constant 0. 0. 0. "
                    "insert_every once overlapcheck yes all_in yes particles_in_region 10000 region reg");
  lammps.input->one("run 1");

  EXPECT_EQ(10000, lammps.atom->natoms);

  Atom * atom = lammps.atom;
  int countA = 0;
  int countB = 0;

  for(int i = 0; i < atom->nlocal; ++i) {
    if(std::abs(atom->radius[i] - 0.0025) < 0.0001) {
      countA++;
    } else if(std::abs(atom->radius[i] - 0.0050) < 0.0001) {
      countB++;
    }
  }

  double VA = (4*M_PI/3.0) * 0.0025*0.0025*0.0025;
  double VB = (4*M_PI/3.0) * 0.005*0.005*0.005;
  double mass_A = countA * 2500 * VA;
  double mass_B = countB * 2500 * VB;

  MPI_Allreduce(MPI_IN_PLACE, &mass_A, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &mass_B, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  double mass_total = mass_A + mass_B;
  double mass_ratioA = mass_A / mass_total;
  double mass_ratioB = mass_B / mass_total;

  EXPECT_TRUE(countA > 0);
  EXPECT_TRUE(countB > 0);
  EXPECT_EQ(atom->nlocal, countA+countB);
  EXPECT_NEAR(0.5, mass_ratioA, 0.01);
  EXPECT_NEAR(0.5, mass_ratioB, 0.01);
}



TEST(particle_distribution_discrete, three_types) {
  const char * argv[3] = {"liggghts", "-in", "scripts/in.insertPackCustom"};
  LAMMPS lammps(3, const_cast<char**>(argv), MPI_COMM_WORLD);
  lammps.input->file();

  lammps.input->one("fix pts1 all particletemplate/sphere 1 atom_type 1 density constant 2500 radius constant 0.0025");
  lammps.input->one("fix pts2 all particletemplate/sphere 1 atom_type 1 density constant 2500 radius constant 0.0050");
  lammps.input->one("fix pts3 all particletemplate/sphere 1 atom_type 1 density constant 2500 radius constant 0.0075");
  lammps.input->one("fix pdd1 all particledistribution/discrete 1.  3 pts1 0.33 pts2 0.33 pts3 0.34");

  lammps.input->one("fix ins all insert/pack seed 100001 distributiontemplate pdd1 vel constant 0. 0. 0. "
                    "insert_every once overlapcheck yes all_in yes particles_in_region 10000 region reg");
  lammps.input->one("run 1");

  EXPECT_EQ(10000, lammps.atom->natoms);

  Atom * atom = lammps.atom;
  int countA = 0;
  int countB = 0;
  int countC = 0;

  for(int i = 0; i < atom->nlocal; ++i) {
    if(std::abs(atom->radius[i] - 0.0025) < 0.0001) {
      countA++;
    } else if(std::abs(atom->radius[i] - 0.0050) < 0.0001) {
      countB++;
    }  else if(std::abs(atom->radius[i] - 0.0075) < 0.0001) {
      countC++;
    }
  }

  double VA = (4*M_PI/3.0) * 0.0025*0.0025*0.0025;
  double VB = (4*M_PI/3.0) * 0.005*0.005*0.005;
  double VC = (4*M_PI/3.0) * 0.0075*0.0075*0.0075;
  double mass_A = countA * 2500 * VA;
  double mass_B = countB * 2500 * VB;
  double mass_C = countC * 2500 * VC;

  MPI_Allreduce(MPI_IN_PLACE, &mass_A, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &mass_B, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &mass_C, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  double mass_total = mass_A + mass_B + mass_C;
  double mass_ratioA = mass_A / mass_total;
  double mass_ratioB = mass_B / mass_total;
  double mass_ratioC = mass_C / mass_total;

  EXPECT_TRUE(countA > 0);
  EXPECT_TRUE(countB > 0);
  EXPECT_TRUE(countC > 0);
  EXPECT_EQ(atom->nlocal, countA+countB+countC);
  EXPECT_NEAR(0.33, mass_ratioA, 0.01);
  EXPECT_NEAR(0.33, mass_ratioB, 0.01);
  EXPECT_NEAR(0.34, mass_ratioC, 0.01);
}

