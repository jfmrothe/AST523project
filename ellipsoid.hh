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
  int GetD();
  gsl_vector * GetCenter();
  gsl_matrix * GetCovMat();
  double GetEnlFac();
  void SetEnlFac(double f);
private:
  int D_;
  gsl_vector * center_;
  gsl_matrix * covMat_;
  double f_;
};

#endif
