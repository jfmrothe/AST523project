#include <ellipsoid.hh>
#include <assert.h>

Ellipsoid::Ellipsoid(int D, gsl_vector * center, gsl_matrix * C, double f) {
  int i;
  D_ = D;
  f_ = f;
  center_ = gsl_vector_alloc(D_);
  gsl_vector_memcpy(center_,center);
  covMat_ = gsl_matrix_alloc(D_,D_);
  gsl_matrix_memcpy(covMat_,C);
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
