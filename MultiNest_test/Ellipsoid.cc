#include "Ellipsoid.h"

float boxmuller()
{
  //returns univariate, zero-mean gaussian sample
  float u1 = UNIFORM;
  float u2 = UNIFORM;
  float pi = 2*acos(0);
  while( u1 == 0.0 ){
    u1 = UNIFORM;
  }
  return sqrt(-2.0*log(u1))*cos(2*pi*u2);
}

float quadr()
{
  //returns sample from quadratic distribution between 0 and 1
  return pow(UNIFORM,1.0/3.0);
}

void unisphere(float * coor, int D)
{
  //returns pseudorandom number uniformly distributed in D-sphere (r=1)
  int i;
  float sample[D];
  float r=0;
  for(i=0;i<D;i++){
    sample[i] = boxmuller();
    r+=pow(sample[i],2);
  }
  r=sqrt(r);
  r/=quadr();
  for(i=0;i<D;i++){
    coor[i] = sample[i]/r;
  }
}

Ellipsoid::Ellipsoid(int D, gsl_vector * center, gsl_matrix * C, double f, vector<Point *>& ell_pts) {
  int i;
  D_ = D;
  f_ = f;
  center_ = gsl_vector_alloc(D_);
  newcoor_ = gsl_vector_alloc(D_);
  gsl_vector_memcpy(center_,center);
  covMat_ = gsl_matrix_alloc(D_,D_);
  gsl_matrix_memcpy(covMat_,C);

  // find own volume
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

  ell_pts_ = ell_pts;

  // ?? add calculation of Cinv here?? 
}

Ellipsoid::~Ellipsoid() {
  gsl_vector_free(center_);
  gsl_vector_free(newcoor_);
  gsl_matrix_free(covMat_);
}

int Ellipsoid::getD() {
  return D_;
}

gsl_vector * Ellipsoid::getCenter() {
  return center_;
}

gsl_matrix * Ellipsoid::getCovMat() {
  return covMat_;
}

double Ellipsoid::getEnlFac() {
  return f_;
}

void Ellipsoid::setEnlFac(double f) {
  // adjust volume
  vol_ *= pow(f/f_,D_);
  // save new enlargement factor
  f_ = f;
}

double Ellipsoid::getVol() {
  return vol_;
}

void Ellipsoid::printout() {
  int i,j;
  printf("%f\n",f_);
  for(j=0;j<D_;j++){
    printf("%f ",gsl_vector_get(center_,j));
  }
  printf("\n"); 
  for(i=0;i<D_;i++){
    for(j=0;j<D_;j++){
      printf("%f ",gsl_matrix_get(covMat_,i,j));
    }
    printf("\n"); 
  }

  return;
}

double Ellipsoid::mdist(Point *pt) {
  int i;
  gsl_vector * tmpvec = gsl_vector_alloc(D_);
  gsl_vector * tmpvec2 = gsl_vector_alloc(D_);

  for(i=0;i<D_;i++){
    // subtract center from each coor to tmpvec
    gsl_vector_set(tmpvec2,i,pt->get_u(i)-gsl_vector_get(center_,i));
  }
 
  //invert covariance matrix ??save Cinv for each ellipsoid??
  int signum = +1;
  gsl_matrix * tmpmat = gsl_matrix_alloc(D_,D_);
  gsl_matrix * Cinv = gsl_matrix_alloc(D_,D_);
  gsl_permutation *p = gsl_permutation_calloc(D_);
  double dist;

  gsl_matrix_memcpy(tmpmat,covMat_);
  gsl_linalg_LU_decomp (tmpmat, p, &signum);
  gsl_linalg_LU_invert (tmpmat, p, Cinv);

  // multiplication
  gsl_blas_dsymv (CblasUpper, 1/f_, Cinv, tmpvec2, 0.0, tmpvec);
  gsl_blas_ddot (tmpvec2, tmpvec, &dist);

  gsl_vector_free(tmpvec);
  gsl_vector_free(tmpvec2);
  gsl_matrix_free(tmpmat);
  gsl_matrix_free(Cinv);
  gsl_permutation_free(p);

  return dist;
}

void Ellipsoid::SampleEllipsoid()
{
  //returns pseudorandom number coor uniformly distributed in ellipsoid given by the center vector center, 
  //covariance matrix C and enlargement factor f, so that x^T(fC)^-1x<=1
  int i;
 
  // create gsl_vector uniformly sampled from sphere
  float spherical[D_];
  unisphere(&spherical[0],D_);
  gsl_vector * spheresample = gsl_vector_alloc(D_);
  for(i=0;i<D_;i++){
    gsl_vector_set(spheresample,i,spherical[i]);
  }

  gsl_matrix * myC = gsl_matrix_calloc(D_,D_);
  gsl_matrix * T = gsl_matrix_calloc(D_,D_);
  gsl_matrix * X = gsl_matrix_calloc(D_,D_);
  gsl_vector * eval = gsl_vector_alloc(D_);
  gsl_matrix * Dprime = gsl_matrix_calloc(D_,D_);

  // allocate workspace for eigensystem calculation 
  gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc(D_);
  // perform eigensystem calculation
  gsl_matrix_memcpy(myC,covMat_);
  gsl_eigen_symmv(myC, eval, X, w);
  // free workspace
  gsl_eigen_symmv_free(w);

  // fill matrix Dprime
  for(i=0;i<D_;i++){
    gsl_matrix_set(Dprime,i,i,sqrt(gsl_vector_get(eval,i)));
  }
  
  // calculate transfer matrix
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, sqrt(f_), X, Dprime, 0.0, T);

  // final coordinates
  gsl_vector_memcpy(newcoor_,center_);
  gsl_blas_dgemv(CblasNoTrans, 1.0, T, spheresample, 1.0, newcoor_);

  gsl_vector_free(spheresample);
  gsl_matrix_free(myC);
  gsl_matrix_free(T);
  gsl_matrix_free(X);
  gsl_vector_free(eval);  
  gsl_matrix_free(Dprime);
}

bool Ellipsoid::IsMember(gsl_vector * coor)
{
  // returns 1 if point specified by coor lies within ellipsoid specified by center, C, f. returns 0 if not

  gsl_permutation *p = gsl_permutation_calloc(D_);
  //  int sign = +1;
  //  int * signum = &sign;
  int i;
  int myone = 1;
  int * signum;
  signum = &myone;
  gsl_matrix * tmp = gsl_matrix_alloc(D_,D_);
  gsl_matrix * Cinv = gsl_matrix_alloc(D_,D_);
  gsl_vector * tmpvec = gsl_vector_alloc(D_);
  gsl_vector * diff = gsl_vector_alloc(D_);
  double range;
  gsl_matrix_memcpy(tmp,covMat_);
  gsl_linalg_LU_decomp (tmp, p, signum);
  gsl_linalg_LU_invert (tmp, p, Cinv);

  for(i=0;i<D_;i++){
    gsl_vector_set(diff,i,gsl_vector_get(coor, i)-gsl_vector_get(center_,i));
  }

  gsl_blas_dsymv (CblasUpper, 1.0, Cinv, diff, 0.0, tmpvec);
  gsl_blas_ddot (diff, tmpvec, &range);

  gsl_permutation_free(p);
  gsl_matrix_free(tmp);
  gsl_matrix_free(Cinv);
  gsl_vector_free(tmpvec);
  gsl_vector_free(diff);

  if (range > f_) return false;
  else return true;
}
