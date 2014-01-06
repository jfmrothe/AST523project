#include <ellipsoid.hh>
#include <math.h>
//#include <assert.h>

#include <sampleellipsoid.hh>

Ellipsoid::Ellipsoid(int D, gsl_vector * center, gsl_matrix * C, double f) {
  int i;
  D_ = D;
  f_ = f;
  center_ = gsl_vector_alloc(D_);
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

  // ?? add calculation of Cinv here?? 
}

Ellipsoid::~Ellipsoid() {
  gsl_vector_free(center_);
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

double Ellipsoid::mdist(gsl_vector * coor) {
  int i;
  gsl_vector * tmpvec = gsl_vector_alloc(D_);
  gsl_vector * tmpvec2 = gsl_vector_alloc(D_);

  for(i=0;i<D_;i++){
    // subtract center from each coor to tmpvec
    gsl_vector_set(tmpvec2,i,gsl_vector_get(coor,i)-gsl_vector_get(center_,i));
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


int TestMdist() {
  // draw random ellipsoid and 100 points from cube and save everything together with m-distances

  int D = 3;
  int Nsample = 1000;
  int i,j;
  gsl_matrix * C = gsl_matrix_alloc(D,D);
  gsl_vector * center = gsl_vector_alloc(D);
  double f;

  GetRandomEllipsoid(D,center,C,&f);
  Ellipsoid myell(D, center, C, f);

  myell.printout();

  gsl_vector * coors[Nsample];
  double mdists[Nsample];
  for(i=0;i<Nsample;i++) {
    coors[i] = gsl_vector_alloc(D);
    for(j=0;j<D;j++) {
      gsl_vector_set(coors[i],j,uniform(0.0,1.0));
    }
    mdists[i] = myell.mdist(coors[i]);
    printf("%f %f %f %f\n",gsl_vector_get(coors[i],0),gsl_vector_get(coors[i],1),gsl_vector_get(coors[i],2),mdists[i]);
  }

  myell.printout();

  for(i=0;i<Nsample;i++) {
    gsl_vector_free(coors[i]);
  }

  return 0;
}
