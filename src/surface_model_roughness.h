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

/* ----------------------------------------------------------------------
   Contributing authors:
   Christoph Kloss (JKU Linz, DCS Computing GmbH, Linz)
   Richard Berger (JKU Linz)
------------------------------------------------------------------------- */
#ifdef SURFACE_MODEL
SURFACE_MODEL(SURFACE_ROUGHNESS,roughness,1)
#else
#ifndef SURFACE_ROUGHNESS_H_
#define SURFACE_ROUGHNESS_H_
#include "contact_models.h"
#include "math.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "fix_roughness.h"

namespace LIGGGHTS {
namespace ContactModels
{
  template<typename Style>
  class SurfaceModel<SURFACE_ROUGHNESS, Style> : protected Pointers
  {
  public:
    static const int MASK = CM_COLLISION;

    SurfaceModel(LAMMPS * lmp, IContactHistorySetup * hsetup) : Pointers(lmp)
    {
      history_offset = hsetup->add_history_value("deltaGamma", "0");
      hsetup->add_history_value("psi", "0");
      hsetup->add_history_value("contactPX", "0");
      hsetup->add_history_value("contactPY", "0");
      hsetup->add_history_value("contactPZ", "0");
      hsetup->add_history_value("contactN1X", "0");
      hsetup->add_history_value("contactN1Y", "0");
      hsetup->add_history_value("contactN1Z", "0");
      /*NL*/ if(comm->me == 0) fprintf(screen, "SURFACE/ROUGHNESS loaded\n");
    }

    inline void registerSettings(Settings&) {}

    inline void connectToProperties(PropertyRegistry &)
    {
      //Check if incompatible rolling friction model is used
      //if (Style::ROLLING == ROLLING_EPSD) {
      //  error->all(FLERR,"PairGranHookeHistoryRoughness cannot handle an epsd-type rolling friction model.");
      //}

      //Get handle to roughness fix
      int i = modify->find_fix("roughness");
      if (i < 0) {
        error->all(FLERR,"Illegal roughness command, need a fix called 'roughness'");
      }

      fix_roughness_ = static_cast<FixRoughness*>(modify->fix[i]);
    }

    inline void collision(CollisionData & cdata, ForceData&, ForceData&)
    {
      static const double piHalf = M_PI / 2.0;

      double enx = cdata.en[0];
      double eny = cdata.en[1];
      double enz = cdata.en[2];
      double deltan = cdata.deltan;

      // relative translational velocity
      const double vr1 = cdata.v_i[0] - cdata.v_j[0];
      const double vr2 = cdata.v_i[1] - cdata.v_j[1];
      const double vr3 = cdata.v_i[2] - cdata.v_j[2];

      // normal component
      double vn = vr1 * enx + vr2 * eny + vr3 * enz;
      double vn1 = vn * enx;
      double vn2 = vn * eny;
      double vn3 = vn * enz;

      // tangential component
      double vt1;
      double vt2;
      double vt3;

      double cosDeltaGamma = 1.0;
      if(cdata.touch) *cdata.touch |= TOUCH_SURFACE_MODEL;
      double * const c_history = &cdata.contact_history[history_offset];

      //****************Handling of Roughness Effect ****************
      if (fix_roughness_->haveNonZeroDeltaGamma())
      {
        if(cdata.is_wall) {
          // Check deltaGamma for collision calculation
          double * x = atom->x[cdata.i];
          double dx = cdata.delta[0];
          double dy = cdata.delta[1];
          double dz = cdata.delta[2];
          double nX = enx;
          double nY = eny;
          double nZ = enz;
          const double nPlaneX = nX;
          const double nPlaneY = nY;
          const double nPlaneZ = nZ;
          double deltaGamma  = c_history[0];
          double psi         = c_history[1];
          double * contactP  = &c_history[2];
          double * contactN1 = &c_history[5];
          double nLength;
          bool isFirstContact = false;

          // tangential component
          vt1 = vr1 - vn1;
          vt2 = vr2 - vn2;
          vt3 = vr3 - vn3;

          // Calculate the impact angle
          double cosAlpha = -vn / sqrt(vr1*vr1 + vr2*vr2 + vr3*vr3);
          if (cosAlpha>1.0) cosAlpha = 1.0; //ensure cosAlpha is smaller or equal unity
          double gamma = piHalf - acos(cosAlpha);
          if (gamma > 0.0) gamma = 0.0; //ensure gamma is larger than zero

          //Set the deltaGamma Angle
          if (deltaGamma == 0.0) //if there was no contact before, we need to randomize the collision angle and save to history
          {
            //return if there is no phyiscal overlap with the (flat) wall
            if(deltan < 0.0) return;

            isFirstContact = true;
            //Generate an appropriate deltaGamma from a Gaussian Distribution
            deltaGamma = fix_roughness_->generateDeltaGamma(gamma);
            psi        = fix_roughness_->generatePsi();

            //Generate the vectors in a plane normal to the contact vector
            //use tangential velocity in the normal plane as the first normal vector
            //and save this vector to contact History
            vt1 += 2e-32; vt2 += 2e-32; vt3 += 1e-32; //add small number to avoid division by 0
            double normVt = sqrt(vt1*vt1 + vt2*vt2 + vt3*vt3);
            contactN1[0]  = -vt1 / normVt;
            contactN1[1]  = -vt2 / normVt;
            contactN1[2]  = -vt3 / normVt;
            nLength = sqrt(contactN1[0]*contactN1[0]
                          +contactN1[1]*contactN1[1]
                          +contactN1[2]*contactN1[2]);  //calculate the length

            // FIXME: Note to original implementor. Fix your code! Floating-point comparison?!?
            if(nLength!=1.0)
              error->all(FLERR,"nLength not equal to unity");
          }

          double n1[3];
          n1[0] = contactN1[0];
          n1[1] = contactN1[1];
          n1[2] = contactN1[2];

          double n2[3];
          n2[0] = nY*n1[2] - nZ*n1[1]; //use cross product to determine second normal vector
          n2[1] =-nX*n1[2] + nZ*n1[0];
          n2[2] = nX*n1[1] - nY*n1[0];
          nLength = sqrt(n2[0]*n2[0] + n2[1]*n2[1] + n2[2]*n2[2]);  //calculate the length

          // FIXME: Note to original implementor. Fix your code! Floating-point comparison?!?
          if(nLength!=1.0)
          {
            fprintf(screen, "deltaGamma, psi: %e %e, contactP: %e %e %e n1: %e %e %e \n",
                          deltaGamma, psi,
                          contactP[0], contactP[1], contactP[2],
                          contactN1[0],  contactN1[1], contactN1[2]);
            error->all(FLERR,"nLength not equal to unity");
          }
          //Generate the vectors in the normal plane for the shift
          const double cosPsi = cos(psi);
          const double sinPsi = sin(psi);
          const double cosDeltaGamma  = cos(deltaGamma);
          const double sinDeltaGamma  = sin(deltaGamma);

          double deltaR[3];
          deltaR[0] = sinDeltaGamma*(sinPsi*n1[0] + cosPsi*n2[0]);
          deltaR[1] = sinDeltaGamma*(sinPsi*n1[1] + cosPsi*n2[1]);
          deltaR[2] = sinDeltaGamma*(sinPsi*n1[2] + cosPsi*n2[2]);

          //Rotate the contact normal to take roughness into account
          //in order to obtain the normal vector for the collision
          //in a (randomly generated) collision plane
          nX =  cosDeltaGamma * nX + deltaR[0];
          nY =  cosDeltaGamma * nY + deltaR[1];
          nZ =  cosDeltaGamma * nZ + deltaR[2];

          if(isFirstContact)
          {
            //Save contact information on the corrected contact point
            dx =  nX * cdata.r;
            dy =  nY * cdata.r;
            dz =  nZ * cdata.r;
            //save a point in the (randomly generated) collision plane
            contactP[0]   = x[0]-dx;
            contactP[1]   = x[1]-dy;
            contactP[2]   = x[2]-dz;
          }

          //calculate the normal distance between the particle
          //and the (randomly generated) collision plane
          dx  =  x[0]-contactP[0];
          dy  =  x[1]-contactP[1];
          dz  =  x[2]-contactP[2];
          //Project into normal direction of n to get
          //the distance to the (randomly generated) collision plane
          const double dxMag = dx*nX + dy*nY + dz*nZ;
          dx    = dxMag*nX;
          dy    = dxMag*nY;
          dz    = dxMag*nZ;

          //re-calculate the distance and the overlap
          cdata.rsq = dx*dx + dy*dy + dz*dz;
          cdata.r = sqrt(cdata.rsq);
          deltan = cdata.radi - cdata.r;

          // store new normal vector and delta vector
          enx = nX;
          eny = nY;
          enz = nZ;
          cdata.en[0] = enx;
          cdata.en[1] = eny;
          cdata.en[2] = enz;
          cdata.delta[0] = dx;
          cdata.delta[1] = dy;
          cdata.delta[2] = dz;

          //Re-Calculate the normal relative velocity
          // normal component
          vn = vr1*enx + vr2*eny + vr3*enz;
          vn1 = vn * enx;
          vn2 = vn * eny;
          vn3 = vn * enz;

              //DEBUG: Save history to file
      #if 0
           FILE * history = fopen("fixWallGranHookeHistoryRoughness.txt","a+");
           fprintf(history,"ip: %d , dx dy dz %g %g %g, r / deltan %g %g c_history %e %e %e extraHistory %g %g %g %g %g %g %g %g \n",
                                         ip,
                                         dx, dy, dz,
                                         r, deltan,
                                         c_history[0],c_history[1],c_history[2],
                                         c_history[3], c_history[4],
                                         c_history[5],c_history[6],c_history[7],
                                         c_history[8],c_history[9],c_history[10] );
           fclose(history);
      #endif

      #if 0
           FILE * history = fopen("fixWallGranHookeHistoryRoughness_2.txt","a+");
           fprintf(history,"i: %d, nX nY nZ %f %f %f, vn %f, vn1 vn2 vn3 %f %f %f, n1/n2/nLength %f %f %f  %f %f %f %f cosAlpha/gamma %f %f, deltaR %f %f %f, deltaGamma / psi %f %f contactP %g %g %g\n",
                                         ip,
                                         nX, nY, nZ,
                                         vn , vn1,vn2,vn3,
                                         n1[0],n1[1],n1[2],
                                         n2[0],n2[1],n2[2], nLength,
                                         cosAlpha, gamma,
                                         deltaR[0],deltaR[1],deltaR[2],
                                         deltaGamma[0], psi[0],
                                         contactP[0],contactP[1],contactP[2] );
           fclose(history);
      #endif
          if(deltan < 0.0)
          {
            //Check if normal velocity is towards wall; if the case,
            //then ensure that particle moves away from the wall
            if (vn > 0.0)
            {
                cdata.v_i[0] += 2.0*vn * nPlaneX;
                cdata.v_i[1] += 2.0*vn * nPlaneY;
                cdata.v_i[2] += 2.0*vn * nPlaneZ;
            }
            return;
          }
        } else {
          // Check deltaGamma for collision calculation
          double nX = -enx;
          double nY = -eny;
          double nZ = -enz;
          double & deltaGamma = c_history[0];
          double & psi        = c_history[1];

          //Calculate the impact angle
          double cosAlpha = -1.0/vn * (vn1*nX + vn2*nY + vn3*nZ);
          if (cosAlpha>1.0) cosAlpha = 1.0; //ensure cosAlpha is smaller or equal unity
          const double gamma =piHalf-acos(cosAlpha);

          //Set the deltaGamma Angle
          if (deltaGamma == 0.0) //if there was no contact before, we need to randomize the collision angle and save to history
          {
            //Generate an appropriate deltaGamma from a Gaussian Distribution
            deltaGamma = fix_roughness_->generateDeltaGamma(gamma);
            psi        = fix_roughness_->generatePsi();
          }

          double n1[3];

          //Generate the vectors in a plane normal to the contact vector
          if (fabs(nZ)>0.1)
          {
            n1[0] = 1.0;
            n1[1] = 0.0;
            n1[2] = -nX/nZ;
          }
          else if (fabs(nY)>0.1)
          {
            n1[0] = 1.0;
            n1[2] = 0.0;
            n1[1] = -nX/nY;
          }
          else
          {
            n1[1] = 1.0;
            n1[2] = 0.0;
            n1[0] = -nY/nX;
          }

          const double nLength = sqrt(n1[0]*n1[0]+n1[1]*n1[1]+n1[2]*n1[2]);  //calculate the length

          //normalize
          n1[0] /= nLength;
          n1[1] /= nLength;
          n1[2] /= nLength;

          //use cross product to determine second normal vector
          double n2[3];
          n2[0] =  nY*n1[2] - nZ*n1[1];
          n2[1] = -nX*n1[2] + nZ*n1[0];
          n2[2] =  nX*n1[1] - nY*n1[0];

          //Generate the vectors in the normal plane for the shift
          const double cosPsi = cos(psi);
          const double sinPsi = sin(psi);
          cosDeltaGamma = cos(deltaGamma);
          const double sinDeltaGamma = sin(deltaGamma);

          double deltaR[3];
          deltaR[0] = sinDeltaGamma*(cosPsi*n1[0] + sinPsi*n2[0]);
          deltaR[1] = sinDeltaGamma*(cosPsi*n1[1] + sinPsi*n2[1]);
          deltaR[2] = sinDeltaGamma*(cosPsi*n1[2] + sinPsi*n2[2]);

          //Rotate the contact normal to take roughness into account
          nX =  cosDeltaGamma * nX + deltaR[0];
          nY =  cosDeltaGamma * nY + deltaR[1];
          nZ =  cosDeltaGamma * nZ + deltaR[2];
          enx = -nX;
          eny = -nY;
          enz = -nZ;

          // store new normal vector and delta vector
          cdata.en[0] = enx;
          cdata.en[1] = eny;
          cdata.en[2] = enz;
          cdata.delta[0] = enx * cdata.r;
          cdata.delta[1] = eny * cdata.r;
          cdata.delta[2] = enz * cdata.r;

          //Re-Calculate the normal relative velocity
          // normal component
          vn = vr1*enx + vr2*eny + vr3*enz;
          vn1 = vn * enx;
          vn2 = vn * eny;
          vn3 = vn * enz;

        //DEBUG: Save history to file
//      FILE * history = fopen("pairGranHookeHistoryRoughness.txt","a+");
////    fprintf(history,"i/j: %d %d,dnum_pairgran: %d , touch[jj] %d, nX nY nZ %f %f %f, vn %f, vn1 vn2 vn3 %f %f %f gamma %f, shear %f %f %f ,  deltaGamma / psi %f %f\n",
////                                   i, j,
////                                   dnum_pairgran,
////                                   touch[jj],
////                                   nX, nY, nZ,
////                                   vn , vn1,vn2,vn3,
////                                   gamma,
////                                   shear[0], shear[1], shear[2], deltaGamma[0], psi[0] );
//      fprintf(history,"i/j: %d %d, nX nY nZ %f %f %f, n1 %f %f %f n2 %f %f %f,  deltaGamma / psi %f %f\n",
//                                   i, j,
//                                   nX, nY, nZ,
//                                   n1[0],n1[1],n1[2],
//                                   n2[0],n2[1],n2[2],
//                                   deltaGamma[0], psi[0] );
//      fclose(history);
        }
      }

      // tangential component
      vt1 = vr1 - vn1;
      vt2 = vr2 - vn2;
      vt3 = vr3 - vn3;

      // relative rotational velocity
      deltan = cdata.radsum - cdata.r;
      const double dx = cdata.delta[0];
      const double dy = cdata.delta[1];
      const double dz = cdata.delta[2];
      const double rinv = cdata.rinv;
      double wr1, wr2, wr3;

      if(cdata.is_wall) {
        // in case of wall contact, r is the contact radius
        const double cr = cdata.radi - 0.5*cdata.deltan;
        wr1 = cr * cdata.omega_i[0] * rinv;
        wr2 = cr * cdata.omega_i[1] * rinv;
        wr3 = cr * cdata.omega_i[2] * rinv;
        cdata.cri = cr;
      } else {
        const double cri = cdata.radi - 0.5 * deltan;
        const double crj = cdata.radj - 0.5 * deltan;
        wr1 = (cri * cdata.omega_i[0] + crj * cdata.omega_j[0]) * rinv;
        wr2 = (cri * cdata.omega_i[1] + crj * cdata.omega_j[1]) * rinv;
        wr3 = (cri * cdata.omega_i[2] + crj * cdata.omega_j[2]) * rinv;
        cdata.cri = cri;
        cdata.crj = crj;
      }

      // relative velocities
      const double vtr1 = vt1 - (dz * wr2 - dy * wr3);
      const double vtr2 = vt2 - (dx * wr3 - dz * wr1);
      const double vtr3 = vt3 - (dy * wr1 - dx * wr2);

      cdata.vn = vn;
      cdata.deltan = deltan;
      cdata.wr1 = wr1;
      cdata.wr2 = wr2;
      cdata.wr3 = wr3;
      cdata.vtr1 = vtr1;
      cdata.vtr2 = vtr2;
      cdata.vtr3 = vtr3;

    }

    inline void noCollision(ContactData&, ForceData&, ForceData&){}
    void beginPass(CollisionData&, ForceData&, ForceData&){}
    void endPass(CollisionData&, ForceData&, ForceData&){}

  protected:
    int history_offset;

    FixRoughness * fix_roughness_;
  };
}
}
#endif // SURFACE_ROUGHNESS
#endif
