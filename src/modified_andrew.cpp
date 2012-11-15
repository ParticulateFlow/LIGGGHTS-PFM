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

#include <algorithm>
#include <cmath>
#include <vector>

using MODIFIED_ANDREW_AUX::Circle;
using MODIFIED_ANDREW_AUX::Line;
using MODIFIED_ANDREW_AUX::Point;
using namespace LAMMPS_NS;

ModifiedAndrew::ModifiedAndrew(){}
ModifiedAndrew::~ModifiedAndrew(){}

double ModifiedAndrew::area(vector<Circle> C){

	//NP TODO: Size of 4 is good?
	if (C.size() < 4){
		return 0.0;
	}

  vector<Circle> hull_c = convex_hull(C);
  vector<Point> hull_p = hull_points(hull_c);

  return area(hull_p);
}

// area of triangle spanned by three points p,m,q
// based on cross product
double ModifiedAndrew::area(Point &p, Point &m, Point &q){
return 0.5*(-(m.y*p.x) + m.x*p.y +
      m.y*q.x - p.y*q.x - m.x*q.y + p.x*q.y);
}

// area of convex hull
double ModifiedAndrew::area(vector<Point> H){
  double a = 0.0;
  Point m = mean_point(H);

  for (int i = 1; i < H.size(); i++){
    a += area(H[i-1],m,H[i]);
  }

  return a;
}

// 2D cross product of OA and OB vectors, i.e. z-component of their 3D cross product.
// Returns a positive value, if OAB makes a counter-clockwise turn,
// negative for clockwise turn, and zero if the points are collinear.
double ModifiedAndrew::turn(const Circle &O, const Circle &A, const Circle &B)
{
  // calculate lines
  Line l1 = find_right_tangent(O,A);
  Line l2 = find_right_tangent(O,B);

  return (-l1.a) * (-l2.b) - (-l1.b) * (-l2.a);
}

// find the right tangent of two circles
Line ModifiedAndrew::find_right_tangent(const Circle &p1, const Circle &p2){
  Line l;
  double d;
  double x,y,r,X,Y,R;

  x = p2.x - p1.x;
  y = p2.y - p1.y;
  r = p2.r - p1.r;
  d = sqrt(x*x + y*y);

  X = x/d;
  Y = y/d;
  R = r/d;

  l.a = R*X - Y*sqrt(1-R*R);
  l.b = R*Y + X*sqrt(1-R*R);
  l.c = p1.r - (l.a*p1.x + l.b*p1.y);

  return l;
}

// find intersection point of two lines that are not parallel
Point ModifiedAndrew::intersectLines(Line &l, Line &j){
  double x,y;

  x = -((-(j.b*l.c) + l.b*j.c) / (j.a*l.b - l.a*j.b));
  y = -((j.a*l.c - l.a*j.c) / (j.a*l.b - l.a*j.b));

  Point p;
  p.x = x;
  p.y = y;

  return p;
}

// find the mean of a vector of points (excluding the first)
Point ModifiedAndrew::mean_point(vector<Point> P){
  Point mean;

  double x_accum, y_accum;

  for (int i = 1; i < P.size(); i++){
    x_accum += P[i].x;
    y_accum += P[i].y;
  }

  mean.x = x_accum/P.size();
  mean.y = y_accum/P.size();

  return mean;
}

// Returns a list of points on the convex hull in counter-clockwise order.
// Note: the last circle in the returned list is the same as the first one.
vector<Circle> ModifiedAndrew::convex_hull(vector<Circle> P)
{
    int n = P.size(), k = 0;
    vector<Circle> H(2*n);

    // Sort points lexicographically
    sort(P.begin(), P.end());

    // Build lower hull
    for (int i = 0; i < n; i++) {
        while (k >= 2 && turn(H[k-2], H[k-1], P[i]) <= 0) k--;
        H[k++] = P[i];
    }

    // Build upper hull
    for (int i = n-2, t = k+1; i >= 0; i--) {
        while (k >= t && turn(H[k-2], H[k-1], P[i]) <= 0) k--;
        H[k++] = P[i];
    }

    H.resize(k);
    return H;
}


// find exterior points of convex hull, i.e. the points
// that are outside the circles, such that their convex hull
// includes all circles, and is the minimal polygon
// with as many edges as hull-circles to do so.
vector<Point> ModifiedAndrew::hull_points(vector<Circle> C){
  double x,y;
  vector<Point> points(C.size());
  Line l1, l2;


  // find intersection around first sphere
  l1 = find_right_tangent(C[C.size()-2],C[0]);
  l2 = find_right_tangent(C[0],C[1]);
  points[0] = intersectLines(l1,l2);

  int N = C.size();
  for (int i = 2; i < N; i++){

    // find tangents
    l1 = find_right_tangent(C[i-2], C[i-1]);
    l2 = find_right_tangent(C[i-1], C[i]);

    // find intersection between tangents
    points[i-1] = intersectLines(l1,l2);
  }

  // copy first point, which equals last point
  N = C.size()-1;
  points[N].x = points[0].x;
  points[N].y = points[0].y;

  return points;
}
