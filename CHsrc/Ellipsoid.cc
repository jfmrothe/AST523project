#include "Ellipsoid.h"

//float boxmuller()
//{
//  //returns univariate, zero-mean gaussian sample
//  //float u1 = UNIFORM;
//  //float u2 = UNIFORM;
//  float u1 = myrand.doub();
//  float u2 = myrand.doub();
//  while( u1 == 0.0 ){
//    //u1 = UNIFORM;
//    u1 = myrand.doub();
//  }
//  return sqrt(-2.0*log(u1))*cos(2*M_PI*u2);
//}
//
//float quadr()
//{
//  //returns sample from quadratic distribution between 0 and 1
//  //return pow(UNIFORM,1.0/3.0);
//  return pow(myrand.doub(),1.0/3.0);
//}

void unisphere(float * coor, int D)
{
  //returns pseudorandom number uniformly distributed in D-sphere (r=1)
  int i;
  float sample[D];
  float r=0;
  Ran myrand(rand()); 
  for(i=0;i<D;i++){

    double u1 = myrand.doub();
    double u2 = myrand.doub();
    while( u1 == 0.0 ){
    u1 = myrand.doub();
     }
   sample[i]=sqrt(-2.0*log(u1))*cos(2*M_PI*u2);
    r+=pow(sample[i],2);
  }
  r=sqrt(r);
  r/= pow(myrand.doub(),1.0/D);
  for(i=0;i<D;i++){
    coor[i] = sample[i]/r;
  }
}

// constructor for EllipsoidalPartitioning
Ellipsoid::Ellipsoid(int D, gsl_vector * center, gsl_matrix * C, double f, vector<Point *>& ell_pts) {
  D_ = D;
  f_ = f;
  center_ = gsl_vector_alloc(D_);
  newcoor_ = gsl_vector_alloc(D_);
  gsl_vector_memcpy(center_,center);
  covMat_ = gsl_matrix_alloc(D_,D_);
  gsl_matrix_memcpy(covMat_,C);
  Cinv_ = gsl_matrix_alloc(D_,D_);
  A_ = gsl_matrix_alloc(D_+1,D_+1);
  Ainv_ = gsl_matrix_alloc(D_+1,D_+1);

 // save inverse covariance matrix
  int signum = +1;
  gsl_matrix * Cinvtmpmat = gsl_matrix_alloc(D_,D_);
  gsl_permutation *Cinvp = gsl_permutation_calloc(D_);

  gsl_matrix_memcpy(Cinvtmpmat,covMat_);
  gsl_linalg_LU_decomp (Cinvtmpmat, Cinvp, &signum);
  gsl_linalg_LU_invert (Cinvtmpmat, Cinvp, Cinv_);

  gsl_matrix_free(Cinvtmpmat);
  gsl_permutation_free(Cinvp);

  // find intersect-check representation matrix
  gsl_matrix * Tself = gsl_matrix_calloc(D_+1,D_+1);
  gsl_matrix * Sself = gsl_matrix_calloc(D_+1,D_+1);
  gsl_matrix * Atmpmat = gsl_matrix_alloc(D_+1,D_+1);
  for(int i=0;i<D_+1;i++) {gsl_matrix_set(Tself,i,i,1.0);}
  for(int i=0;i<D_;i++) {gsl_matrix_set(Tself,D_,i,-1.0*gsl_vector_get(center_,i));}
  for(int i=0;i<D_+1;i++) {gsl_matrix_set(Sself,i,i,-1.0);}
  for(int i=0;i<D_;i++) {
    for(int j=0;j<D_;j++) {
      gsl_matrix_set(Sself,i,j,gsl_matrix_get(Cinv_,i,j)/f_);
    }
  }
  gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1.0, Sself, Tself, 0.0, Atmpmat);
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, Tself, Atmpmat, 0.0, A_);

  gsl_matrix_free(Atmpmat);//
  gsl_matrix_free(Tself);
  gsl_matrix_free(Sself);

  // invert check-intersect-representation matrix
  signum = +1;
  gsl_permutation * Ainvp = gsl_permutation_calloc(D_+1);
  gsl_matrix * Ainvtmpmat = gsl_matrix_alloc(D_+1,D_+1);//

  gsl_matrix_memcpy(Ainvtmpmat,A_);
  gsl_linalg_LU_decomp (Ainvtmpmat, Ainvp, &signum);
  gsl_linalg_LU_invert (Ainvtmpmat, Ainvp, Ainv_);

  gsl_matrix_free(Ainvtmpmat);
  gsl_permutation_free(Ainvp);

  // find own volume
  gsl_matrix * voltmpmat = gsl_matrix_alloc(D,D);
  double det;
  signum = +1;
  gsl_permutation * volp = gsl_permutation_calloc(D);
  // use gsl to compute determinant
  gsl_matrix_memcpy(voltmpmat,covMat_);
  gsl_linalg_LU_decomp(voltmpmat, volp, &signum);
  det =  gsl_linalg_LU_det (voltmpmat, signum);

  vol_ = 4.0/3.0*M_PI*sqrt(pow(f_,D)*det);


  //  ell_pts_ = ell_pts;
  // Chelsea's suggested problem: as ell_pts goes out of scope, Points referenced within are auto-destroyed
  for(unsigned int i=0;i<ell_pts.size();i++) {
    ell_pts_.push_back(new Point(*ell_pts[i]));
  }

  gsl_matrix_free(voltmpmat); 
  gsl_permutation_free(volp);
}

Ellipsoid::~Ellipsoid() {
  gsl_vector_free(center_);
  gsl_vector_free(newcoor_);
  gsl_matrix_free(covMat_);
  gsl_matrix_free(Cinv_);
  gsl_matrix_free(A_);
  gsl_matrix_free(Ainv_);
  //int size = ell_pts_.size();
  for(unsigned int j=0; j<ell_pts_.size(); j++){delete ell_pts_[j];}
  ell_pts_.clear();
  //delete ell_pts_;
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
  vol_ *= pow(f/f_,1.0*D_);
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
  double dist;
  gsl_vector * tmpvec = gsl_vector_alloc(D_);
  gsl_vector * tmpvec2 = gsl_vector_alloc(D_);

  for(i=0;i<D_;i++){
    // subtract center from each coor to tmpvec
    gsl_vector_set(tmpvec,i,pt->get_u(i)-gsl_vector_get(center_,i));
  }
 
  // multiplication
  gsl_blas_dsymv (CblasUpper, 1/f_, Cinv_, tmpvec, 0.0, tmpvec2);
  gsl_blas_ddot (tmpvec, tmpvec2, &dist);

  gsl_vector_free(tmpvec);
  gsl_vector_free(tmpvec2);

  return dist;
}

void Ellipsoid::SampleEllipsoid()
{
  //creates pseudorandom number coor uniformly distributed in ellipsoid given by the center vector center, 
  //covariance matrix C and enlargement factor f, so that x^T(fC)^-1x<=1
  //stores is in new_coor member variable
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
  //gsl_vector_memcpy(newcoor,newcoor_);

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



bool Ellipsoid::intersect(Ellipsoid& other) {

  if(other.getD() != D_) {return false;}

  gsl_matrix * testMat = gsl_matrix_alloc(D_+1,D_+1);
  // multiply others A with own inverse A into testMat
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, other.getA(), Ainv_, 0.0, testMat);
  // obtain eigenvalues
  gsl_vector_complex * eval = gsl_vector_complex_alloc(D_+1);
  gsl_eigen_nonsymm_workspace * w = gsl_eigen_nonsymm_alloc (D_+1);
  gsl_eigen_nonsymm (testMat, eval, w);
 
  bool result;
  int negrealcount = 0;
  // count negative real eigenvalues
  for(int i=0;i<D_+1;i++){ 
    if(GSL_REAL(gsl_vector_complex_get(eval,i))<0.0 and GSL_IMAG(gsl_vector_complex_get(eval,i))==0.0) {negrealcount++;}
  }
  // if exactly two are real and negative, result is no -- else result is yes
  if (negrealcount == 2) {result=false;}
  else{result = true;}

  gsl_matrix_free(testMat);
  gsl_eigen_nonsymm_free(w);
  gsl_vector_complex_free(eval);
  return result;
}

void Ellipsoid::fetchPoints(Ellipsoid& other) {

  for(int j=0; j<other.ell_pts_.size(); j++) {
    ell_pts_.push_back(other.ell_pts_[j]);    
  }

}


void Ellipsoid::RescaleToCatch() {
    // rescale to catch all points
    double tmp;
    double f = 0.0;
    int Npoints = ell_pts_.size();
    gsl_vector * tmpvec = gsl_vector_alloc(D_);
    gsl_vector * tmpvec2 = gsl_vector_alloc(D_);
    
    for(int j=0;j<Npoints;j++){
      for(int i=0;i<D_;i++){
        // subtract center from each coor to tmpvec
        gsl_vector_set(tmpvec2,i,ell_pts_[j]->get_u(i)-gsl_vector_get(center_,i));
      }
      gsl_blas_dsymv (CblasUpper, 1.0, Cinv_, tmpvec2, 0.0, tmpvec);
      gsl_blas_ddot (tmpvec2, tmpvec, &tmp);
      if (tmp > f) f = tmp;
    }

    if (f>f_) {
      setEnlFac(f);
    }
    gsl_vector_free(tmpvec);
    gsl_vector_free(tmpvec2);
}

