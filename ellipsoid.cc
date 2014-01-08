#include <ellipsoid.hh>
#include <math.h>
//#include <assert.h>

#include <sampleellipsoid.hh>

Ellipsoid::Ellipsoid(int D, gsl_vector * center, gsl_matrix * C, double f) {
  int i;
  // copy member variables
  D_ = D;
  f_ = f;
  center_ = gsl_vector_alloc(D_);
  gsl_vector_memcpy(center_,center);
  covMat_ = gsl_matrix_alloc(D_,D_);
  gsl_matrix_memcpy(covMat_,C);

  // save inverse covariance matrix
  Cinv_ = gsl_matrix_alloc(D_,D_);
  
  int signum = +1;
  gsl_matrix * tmpmat = gsl_matrix_alloc(D_,D_);
  gsl_permutation *p = gsl_permutation_calloc(D_);

  gsl_matrix_memcpy(tmpmat,covMat_);
  gsl_linalg_LU_decomp (tmpmat, p, &signum);
  gsl_linalg_LU_invert (tmpmat, p, Cinv_);

  gsl_matrix_free(tmpmat);
  gsl_permutation_free(p);

  // find intersect-check representation matrix
  gsl_matrix * A_ = gsl_matrix_alloc(D+1,D+1);
  gsl_matrix * Tself = gsl_matrix_alloc(D+1,D+1);
  gsl_matrix * Sself = gsl_matrix_alloc(D+1,D+1);
  tmpmat = gsl_matrix_alloc(D+1,D+1);
  for(int i=0;i<D+1;i++) {gsl_matrix_set(Tself,i,i,1.0);}
  for(int i=0;i<D;i++) {gsl_matrix_set(Tself,D,i,-1.0*gsl_vector_get(center_,i));}
  for(int i=0;i<D+1;i++) {gsl_matrix_set(Sself,i,i,-1.0);}
  for(int i=0;i<D;i++) {
    for(int j=0;j<D;j++) {
      gsl_matrix_set(Sself,i,j,gsl_matrix_get(Cinv_,i,j));
    }
  }
  gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1.0, Sself, Tself, 0.0, tmpmat);
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, Tself, tmpmat, 0.0, A_);

  gsl_matrix_free(Tself);
  gsl_matrix_free(Sself);

  // invert check-intersect-representation matrix
  gsl_matrix * Ainv_ = gsl_matrix_alloc(D_+1,D_+1);

  signum = +1;
  p = gsl_permutation_calloc(D_+1);

  gsl_matrix_memcpy(tmpmat,A_);
  gsl_linalg_LU_decomp (tmpmat, p, &signum);
  gsl_linalg_LU_invert (tmpmat, p, Ainv_);

  gsl_matrix_free(tmpmat);
  gsl_permutation_free(p);

  // find own volume
  tmpmat = gsl_matrix_alloc(D,D);
  float pi = 2*acos(0);
  double det;
  signum = +1;
  p = gsl_permutation_calloc(D);
  // use gsl to compute determinant
  gsl_matrix_memcpy(tmpmat,covMat_);
  gsl_linalg_LU_decomp(tmpmat, p, &signum);
  det =  gsl_linalg_LU_det (tmpmat, signum);

  vol_ = 4.0/3.0*pi*sqrt(pow(f_,D)*det);

  gsl_matrix_free(tmpmat); 
  gsl_permutation_free(p);

}

Ellipsoid::~Ellipsoid() {
  gsl_vector_free(center_);
  gsl_matrix_free(covMat_);
  gsl_matrix_free(Cinv_);
  gsl_matrix_free(A_);
  gsl_matrix_free(Ainv_);
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

gsl_matrix * Ellipsoid::getCinv() {
  return Cinv_;
}

gsl_matrix * Ellipsoid::getA() {
  return A_;
}

gsl_matrix * Ellipsoid::getAinv() {
  return Ainv_;
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


bool Ellipsoid::intersect(Ellipsoid& other) {

  if(other.getD() != D_) {return false;}

  gsl_matrix * testMat = gsl_matrix_alloc(D_+1,D_+1);
  // multiply others A with own inverse A
  gsl_matrix_get(testMat,0,0);//ok
  gsl_matrix_get(this->getA(),0,0);
  gsl_matrix_get(A_,0,0);
  gsl_matrix_get(other.getA(),0,0);
  gsl_matrix_get(other.getAinv(),0,0);
  gsl_matrix_get(Ainv_,0,0);
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, other.getA(), Ainv_, 0.0, testMat);
  // obtain eigenvalues
  gsl_vector_complex * eval = gsl_vector_complex_alloc(D_+1);
  gsl_eigen_nonsymm_workspace * w = gsl_eigen_nonsymm_alloc (D_+1);
  gsl_eigen_nonsymm (testMat, eval, w);
  for(int i=0;i<D_+1;i++){ 
    printf("%f+%fi ",GSL_REAL(gsl_vector_complex_get(eval,i)),GSL_IMAG(gsl_vector_complex_get(eval,i)));
  }
  printf("\n");

  bool result;
  // if exactly two are real and negative, result is no

  // else result is yes
  result = true;

  gsl_matrix_free(testMat);
  gsl_eigen_nonsymm_free(w);
  gsl_vector_complex_free(eval);
  return result;
}

int TestIntersect() {
  int D=3;
  double f1;
  gsl_matrix * C1 = gsl_matrix_alloc(D,D);
  gsl_vector * ctr1 = gsl_vector_alloc(D);
  double f2;
  gsl_matrix * C2 = gsl_matrix_alloc(D,D);
  gsl_vector * ctr2 = gsl_vector_alloc(D);
  GetRandomEllipsoid(D,ctr1,C1,&f1);
  Ellipsoid first(D,ctr1,C1,f1);
  GetRandomEllipsoid(D,ctr2,C2,&f2);
  Ellipsoid second(D,ctr2,C2,f2);
  first.printout();
  second.printout();
  printf("%i\n",first.intersect(second));
  return 0;
}
