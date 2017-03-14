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
   Stefan Amberger (JKU Linz)
------------------------------------------------------------------------- */

#include "modified_andrew.h"
#include "mpi_liggghts.h"
#include "vector_liggghts.h"
#include "error.h"
#include "comm.h"
/*NL*/ // #include "update.h"

#include <algorithm>
#include <cmath>
#include <vector>

using MODIFIED_ANDREW_AUX::Circle;
using MODIFIED_ANDREW_AUX::Point;
using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ModifiedAndrew::ModifiedAndrew(LAMMPS *lmp) : Pointers(lmp) {
  npoints_per_circle_ = 6;
  directions_ = new double*[npoints_per_circle_];
  for (int i = 0; i < npoints_per_circle_; ++i)
  {
    directions_[i] = new double[2];
  }
  double angle = (2.0*M_PI)/((double)npoints_per_circle_);
  for (int i = 0; i < npoints_per_circle_; ++i)
  {
    directions_[i][0] = sin(((double)i)*angle);
    directions_[i][1] = cos(((double)i)*angle);
  }

  /*NL*/ // // DEBUG
  /*NL*/ // for (int j = 0; j < comm->nprocs; ++j)
  /*NL*/ // {
  /*NL*/ //   if (comm->me == j){
  /*NL*/ //     printf("directions proc %d:\n",comm->me);
  /*NL*/ //     for (int i = 0; i < npoints_per_circle_; ++i)
  /*NL*/ //       printf("(%f %f)\n", directions_[i][0],directions_[i][1]);
  /*NL*/ //     }
  /*NL*/ //   MPI_Barrier(world);
  /*NL*/ // }
  /*NL*/ // // END DEBUG

}
ModifiedAndrew::~ModifiedAndrew(){
  for (int i = 0; i < npoints_per_circle_; ++i)
    delete []directions_[i];
  delete []directions_;
}

/* ---------------------------------------------------------------------- */

double ModifiedAndrew::area()
{
    double A;
    vector<Point> hull_c = convex_hull(contacts_);

    /*NL*/ //if (screen) {
    /*NL*/ //fprintf(screen,"hull_c size %d\n",hull_c.size());
    /*NL*/ //for(int i = 0; i < hull_c.size(); i++)
    /*NL*/ //  fprintf(screen,"element %d: %f %f %f\n",i,hull_c[i].x,hull_c[i].y,hull_c[i].r);
    /*NL*/ //}

    // multi-proc case
    if(1 < comm->nprocs) //(1)//
    {
      /*NL*/ // // DEBUG
      /*NL*/ // for (int i = 0; i < comm->nprocs; ++i)
      /*NL*/ // {
      /*NL*/ //   if (comm->me == i)
      /*NL*/ //   {
      /*NL*/ //     printf("contacts_%d=[\n", i);
      /*NL*/ //     for (int j = 0; j < contacts_.size(); ++j)
      /*NL*/ //     {
      /*NL*/ //       printf("(%f, %f),\n", contacts_[j].x, contacts_[j].y);
      /*NL*/ //     }
      /*NL*/ //     printf("]\n");
      /*NL*/ //     fflush(stdout);
      /*NL*/ //   }
      /*NL*/ //   MPI_Barrier(world);
      /*NL*/ // }
      /*NL*/ //   MPI_Barrier(world);
      /*NL*/ // // END DEBUG

      /*NL*/ //  // DEBUG
      /*NL*/ //  for (int i = 0; i < comm->nprocs; ++i)
      /*NL*/ //  {
      /*NL*/ //    if (comm->me == i)
      /*NL*/ //    {
      /*NL*/ //      printf("hull_c%d=[\n", i);
      /*NL*/ //      for (int j = 0; j < hull_c.size(); ++j)
      /*NL*/ //      {
      /*NL*/ //        printf("(%f, %f),\n", hull_c[j].x, hull_c[j].y);
      /*NL*/ //      }
      /*NL*/ //      printf("]\n");
      /*NL*/ //      fflush(stdout);
      /*NL*/ //    }
      /*NL*/ //    MPI_Barrier(world);
      /*NL*/ //  }
      /*NL*/ //    MPI_Barrier(world);
      /*NL*/ //  // END DEBUG


        // perform allgather
        double *data = 0,*data0 = 0;
        int ndata = 0, ndata0 = 0;

        //NP erase first element which is duplicate
        if(hull_c.size() > 2)
            hull_c.erase(hull_c.begin());

        // copy from STL vector to C array to safely use MPI
        ndata = construct_data(hull_c,data);
        ndata0 = MPI_Gather0_Vector(data,ndata,data0,world);

        /*NL*///if (screen) printVecN(screen,"data",data,ndata);
        /*NL*///if (screen) printVecN(screen,"data0",data0,ndata0);

        // proc0 calculates convex hull
        if(0 == comm->me) {

            vector<Point> hull_c_allreduced = construct_hull_c_all(data0,ndata0);
            vector<Point> hull_c_global = convex_hull(hull_c_allreduced);

            /*NL*/ //  // DEBUG
            /*NL*/ //      printf("hull_c_global=[\n");
            /*NL*/ //      for (int j = 0; j < hull_c_global.size(); ++j)
            /*NL*/ //        printf("(%f, %f),\n", hull_c_global[j].x, hull_c_global[j].y);
            /*NL*/ //      printf("]\n");
            /*NL*/ //      fflush(stdout);
            /*NL*/ //  // END DEBUG


            /*NL*/ //if (screen) {
            /*NL*/ //fprintf(screen,"hull_c_allreduced size %d\n",hull_c_allreduced.size());
            /*NL*/ //for(int i = 0; i < hull_c_allreduced.size(); i++)
            /*NL*/ //  fprintf(screen,"element %d: %f %f %f\n",i,hull_c_allreduced[i].x,hull_c_allreduced[i].y,hull_c_allreduced[i].r);
            /*NL*/ //}

            if (hull_c_global.size() < 3)
                A = 0.;
            else
                A = area(hull_c_global);
        }
        /*NL*/ // // DEBUG
        /*NL*/ // MPI_Barrier(world);
        /*NL*/ // // END DEBUG

        if(data) delete []data;
        if(data0) delete []data0;

        /*NL*/ // // DEBUG
        /*NL*/ // for (int i = 0; i < comm->nprocs; ++i)
        /*NL*/ // {
        /*NL*/ //   if (comm->me == i)
        /*NL*/ //     printf("A%d=%f\n", i, A);
        /*NL*/ //   MPI_Barrier(world);
        /*NL*/ // }
        /*NL*/ //   MPI_Barrier(world);
        /*NL*/ // // END DEBUG


        // broadcast result
        MPI_Bcast(&A,1,MPI_DOUBLE,0,world);
    }
    // single proc case
    else
    {
        if (hull_c.size() < 3)
           A = 0.;
        else
            A = area(hull_c);
    }

    // clear contact data
    contacts_.clear();

    /*NL*/ // DEBUG
    /*NL*/ // if (update->ntimestep == 5000)
    /*NL*/ //   error->all(FLERR, "reached timestep 5000");
    /*NL*/ // END DEBUG

    return A;
}

/* ---------------------------------------------------------------------- */

void ModifiedAndrew::add_contact(Circle c){
  // generate npoints_per_circle_ points for this circle
  Point p;
  for (int i = 0; i < npoints_per_circle_; ++i)
  {
    p.x = c.x+directions_[i][0]*c.r;
    p.y = c.y+directions_[i][1]*c.r;
    contacts_.push_back(p);
  }
}

/* ---------------------------------------------------------------------- */

int ModifiedAndrew::construct_data(vector<Point> hull_c, double *&data)
{
    int size = hull_c.size();
    int datasize = 2*size;
    data = new double[datasize];

    for(int i = 0; i < size; i++)
    {
        data[2*i+0] = hull_c[i].x;
        data[2*i+1] = hull_c[i].y;
    }

    return datasize;
}

/* ---------------------------------------------------------------------- */

vector<Point> ModifiedAndrew::construct_hull_c_all(double *data0, int ndata0)
{
    vector<Point> result;

    Point c;

    for(int i = 0; i < ndata0/2; i++)
    {
        c.x = data0[i*2+0];
        c.y = data0[i*2+1];
        result.push_back(c);
    }

    /*NL*/ //if (screen) fprintf(screen,"ndata0 %d result.size() %d\n",ndata0,result.size());

    return result;
}

/* ---------------------------------------------------------------------- */

// area of triangle spanned by three points p,m,q
// based on cross product
double ModifiedAndrew::area(Point &p, Point &m, Point &q){
return 0.5*(-(m.y*p.x) + m.x*p.y +
      m.y*q.x - p.y*q.x - m.x*q.y + p.x*q.y);
}

/* ---------------------------------------------------------------------- */

// area of convex hull
double ModifiedAndrew::area(vector<Point> H){
  double a = 0.0;
  Point m = mean_point(H);

  for (size_t i = 1; i < H.size(); i++){
    a += area(H[i-1],m,H[i]);
  }

  return a;
}

/* ---------------------------------------------------------------------- */

// 2D cross product of OA and OB vectors, i.e. z-component of their 3D cross product.
// Returns a positive value, if OAB makes a counter-clockwise turn,
// negative for clockwise turn, and zero if the points are collinear.
double ModifiedAndrew::cross(Point O, Point A, Point B)
{
  return (A.x - O.x)*(B.y - O.y) - (A.y - O.y)*(B.x - O.x);
}


/* ---------------------------------------------------------------------- */

// find the mean of a vector of points (excluding the first)
Point ModifiedAndrew::mean_point(vector<Point> P){
  Point mean;

  double x_accum = 0.0, y_accum = 0.0;

  for (size_t i = 1; i < P.size(); i++){
    x_accum += P[i].x;
    y_accum += P[i].y;
  }

  mean.x = x_accum/(static_cast<double>(P.size())-1.);
  mean.y = y_accum/(static_cast<double>(P.size())-1.);

  return mean;
}

/* ---------------------------------------------------------------------- */

// Returns a list of points on the convex hull in counter-clockwise order.
// Note: the last circle in the returned list is the same as the first one.
vector<Point> ModifiedAndrew::convex_hull(vector<Point> P)
{
    int n = P.size(), k = 0;
    vector<Point> H(2*n);

    // Sort points lexicographically
    sort(P.begin(), P.end());

    // Build lower hull
    for (int i = 0; i < n; i++) {
        while (k >= 2 && cross(H[k-2], H[k-1], P[i]) <= 0) k--;
        H[k++] = P[i];
    }

    // Build upper hull
    for (int i = n-2, t = k+1; i >= 0; i--) {
        while (k >= t && cross(H[k-2], H[k-1], P[i]) <= 0) k--;
        H[k++] = P[i];
    }

    H.resize(k);
    return H;
}
