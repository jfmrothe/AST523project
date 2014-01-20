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
  void SampleEllipsoid(gsl_vector * target);
  void RescaleToCatch();
  bool IsMember(gsl_vector *);
  bool intersect(Ellipsoid& other);
  void fetchPoints(Ellipsoid& other);
  double getEnlFac();
  void setEnlFac(double f);
  double getVol();
  double mdist(Point *);
  void printout();
  vector<Point *> ell_pts_;
private:
  int D_;
  gsl_vector * center_;
  gsl_matrix * Cinv_;
  gsl_matrix * A_;
  gsl_matrix * Ainv_;
  gsl_matrix * covMat_;
  double f_;
  double vol_;

};
#endif
