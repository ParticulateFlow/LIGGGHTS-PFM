/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   Christoph Kloss, christoph.kloss@cfdem.com
   Copyright 2009-2012 JKU Linz
   Copyright 2012-     DCS Computing GmbH, Linz

   LIGGGHTS is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(mesh/surface/stress/contact,FixMeshSurfaceStressContact)

#else

#ifndef LMP_FIX_MESH_SURFACE_STRESS_CONTACT_H
#define LMP_FIX_MESH_SURFACE_STRESS_CONTACT_H

#include "fix_mesh_surface_stress.h"

namespace LAMMPS_NS
{
  class FixMeshSurfaceStressContact : public FixMeshSurfaceStress
  {
      public:

        FixMeshSurfaceStressContact(LAMMPS *lmp, int narg, char **arg);
        virtual ~FixMeshSurfaceStressContact();

        virtual void post_create();
        virtual void init();
        virtual int setmask();

        void pre_force(int vflag);
        void final_integrate();

        void add_particle_contribution(int ip, double *frc,
                            double *delta, int iTri, double *v_wall);

      protected:

        // inline access

        inline double& contactAreaStep(int i)
        { return (*contact_area_step_)(i); }

        inline double& contactArea(int i)
        { return (*contact_area_)(i); }

      private:

        void init_area_correction();
        void calc_contact_time_average();

        class FixPropertyAtom *fix_wallcontacttime_;
        ScalarContainer<double> *contact_area_;
        ScalarContainer<double> *contact_area_step_;

        //NP time for averaging
        double T_;

        //NP time step for starting averageing
        int step_ave_start_;

        bool area_correction_;
        double const* const* deltan_ratio_;
  };

} /* namespace LAMMPS_NS */

#endif
#endif
