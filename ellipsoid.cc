#include <ellipsoid.hh>
#include <math.h>
//#include <assert.h>


Ellipsoid::Ellipsoid(int D, gsl_vector * center, gsl_matrix * C, double f) {
  int i;
  D_ = D;
  f_ = f;
  center_ = gsl_vector_alloc(D_);
  gsl_vector_memcpy(center_,center);
  covMat_ = gsl_matrix_alloc(D_,D_);
  gsl_matrix_memcpy(covMat_,C);

  gsl_matrix * tmpMat = gsl_matrix_alloc(D,D);
  float pi = 2*acos(0);
  double det;
  int signum = +1;
  gsl_permutation * p = gsl_permutation_calloc(D);
  // use gsl to compute determinant
  gsl_matrix_memcpy(tmpMat,covMat_);
  gsl_linalg_LU_decomp(tmpMat, p, &signum);
  det =  gsl_linalg_LU_det (tmpMat, signum);

  vol_ = 4.0/3.0*pi*sqrt(pow(f_,D)*det);

  gsl_matrix_free(tmpMat); 
  gsl_permutation_free(p); 
}

Ellipsoid::~Ellipsoid() {
  gsl_vector_free(center_);
  gsl_matrix_free(covMat_);
}

int Ellipsoid::GetD() {
  return D_;
}

gsl_vector * Ellipsoid::GetCenter() {
  return center_;
}

gsl_matrix * Ellipsoid::GetCovMat() {
  return covMat_;
}

double Ellipsoid::GetEnlFac() {
  return f_;
}

void Ellipsoid::SetEnlFac(double f) {
  f_ = f;
}

double Ellipsoid::GetVol() {
  return vol_;
}
