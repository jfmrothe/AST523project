# ifndef _AST523_GSL_LIB_
# define _AST523_GSL_LIB_
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
# endif

#ifndef ELLIPSOID_H_
#define ELLIPSOID_H_

class Ellipsoid {
public:
  Ellipsoid(int D, gsl_vector * center, gsl_matrix * C, double f);
  ~Ellipsoid();
  int getD();
  gsl_vector * getCenter();
  gsl_matrix * getCovMat();
  gsl_matrix * getCinv();
  gsl_matrix * getA();
  gsl_matrix * getAinv();
  double getEnlFac();
  void setEnlFac(double f);
  double getVol();
  double mdist(gsl_vector * coor);
  void printout();
  bool intersect(Ellipsoid& other);
private:
  int D_;
  gsl_vector * center_;
  gsl_matrix * covMat_;
  double f_;
  double vol_;
  gsl_matrix * Cinv_;
  gsl_matrix * A_;
  gsl_matrix * Ainv_;
};

int TestMdist();
int TestIntersect();
#endif
