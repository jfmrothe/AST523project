#ifndef ELLIPSOID_H
#define ELLIPSOID_H

#include <list>
#include <fstream>
#include <string>
#include <float.h>
#include <math.h>
#include <vector>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
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
  void RescaleToCatch();
  bool IsMember(gsl_vector *);
  bool intersect(Ellipsoid& other);
  void fetchPoints(Ellipsoid& other);
  double getEnlFac();
  void setEnlFac(double f);
  double getVol();
  double mdist(Point *);
  void printout();
  gsl_vector * get_newcoor() {return newcoor_;}
  vector<Point *> ell_pts_;
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

};
#endif
