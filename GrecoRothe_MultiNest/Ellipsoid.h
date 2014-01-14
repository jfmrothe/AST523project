/***************************************************

File: Ellipsoid.h 

Description:
Header file for the Ellipsoid class. Ellipsoid objects
are used to represent MultiNest's clustering of the 
active points. They keep track of their own position
and dimensions in MultiNest's native hypercube parameter
space and own a list of the Point objects they enclose.
Important methods are the Mahalanobis distance needed 
for optimization of the clustering, and a sampling method
that returns a random coordinate uniformly distributed
in the ellipsoid.

Programmers: Johnny Greco & Johannes Rothe
Contacts: greco@princeton.edu, jrothe@princeton.edu

****************************************************/
#ifndef ELLIPSOID_H
#define ELLIPSOID_H

#include "Point.h"

float boxmuller();
float quadr();
void unisphere(float *, int);

class Ellipsoid {
public:
  Ellipsoid(int D, gsl_vector * center, gsl_matrix * C, double f, vector<Point *>&);
  ~Ellipsoid();
  int getD();
  gsl_vector * getCenter();
  gsl_matrix * getCovMat();
  gsl_matrix * getCinv();
  gsl_matrix * getA();
  gsl_matrix * getAinv();
  void SampleEllipsoid();
  bool IsMember(gsl_vector *);
  bool intersect(Ellipsoid& other);
  double getEnlFac();
  void setEnlFac(double f);
  double getVol();
  double mdist(Point *);
  void printout();
  gsl_vector * get_newcoor() {return newcoor_;}
private:
  int D_;
  gsl_vector * center_;
  gsl_matrix * Cinv_;
  gsl_matrix * A_;
  gsl_matrix * Ainv_;
  gsl_matrix * covMat_;
  gsl_vector * newcoor_;
  double f_;
  double vol_;
  vector<Point *> ell_pts_;

};
#endif
