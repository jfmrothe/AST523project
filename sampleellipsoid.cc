#include <sampleellipsoid.hh>
#include <findenclosing.hh>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
float uniform(float min, float max){
  //returns pseudorandom number from uniform distribution
  return (((float) rand())/((float) RAND_MAX))*(max-min) + min;
}

float boxmuller(){
  //returns univariate, zero-mean gaussian sample
  float u1 = uniform(0,1);
  float u2 = uniform(0,1);
  float pi = 2*acos(0);
  while( u1 == 0.0 ){
    u1 = uniform(0,1);
  }
  return sqrt(-2.0*log(u1))*cos(2*pi*u2);
}

float quadr(){
  //returns sample from quadratic distribution between 0 and 1
  return pow(uniform(0,1),1.0/3.0);
}

void unisphere(float * coor, int D){
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

void SampleEllipsoid(int D, gsl_vector * center, gsl_matrix * C, double f, gsl_vector * coor){
  //returns pseudorandom number coor uniformly distributed in ellipsoid given by the center vector center, covariance matrix C and enlargement factor f, so that x^T(fC)^-1x<=1
  int i;

  // create gsl_vector uniformly sampled from sphere
  float spherical[D];
  unisphere(&spherical[0],D);
  gsl_vector * spheresample = gsl_vector_alloc(D);
  for(i=0;i<D;i++){
    gsl_vector_set(spheresample,i,spherical[i]);
  }

  gsl_matrix * myC = gsl_matrix_calloc(D,D);
  gsl_matrix * T = gsl_matrix_calloc(D,D);
  gsl_matrix * X = gsl_matrix_calloc(D,D);
  gsl_vector * eval = gsl_vector_alloc(D);
  gsl_matrix * Dprime = gsl_matrix_calloc(D,D);

  // allocate workspace for eigensystem calculation 
  gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc(D);
  // perform eigensystem calculation
  gsl_matrix_memcpy(myC,C);
  gsl_eigen_symmv(myC, eval, X, w);
  // free workspace
  gsl_eigen_symmv_free(w);

  // fill matrix Dprime
  for(i=0;i<D;i++){
    gsl_matrix_set(Dprime,i,i,sqrt(gsl_vector_get(eval,i)));
  }
  
  // calculate transfer matrix
  // this is where the dog lies buried
  //gsl_blas_dgemm(CblasTrans, CblasNoTrans, f, X, Dprime, 0.0, T);
  //gsl_blas_dgemm(CblasTrans, CblasNoTrans, sqrt(f), X, Dprime, 0.0, T);
gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, sqrt(f), X, Dprime, 0.0, T);

  // final coordinates
  gsl_vector_memcpy(coor,center);
  gsl_blas_dgemv(CblasNoTrans, 1.0, T, spheresample, 1.0, coor);
  
  gsl_vector_free(spheresample);
  gsl_matrix_free(myC);
  gsl_matrix_free(T);
  gsl_matrix_free(X);
  gsl_vector_free(eval);  
  gsl_matrix_free(Dprime);
  gsl_eigen_symmv_free(w);
}

int IsMember(int D, gsl_vector * coor,  gsl_vector * center, gsl_matrix * C, double f){
  // returns 1 if point specified by coor lies within ellipsoid specified by center, C, f. returns 0 if not
  gsl_permutation *p = gsl_permutation_calloc(D);
  //  int sign = +1;
  //  int * signum = &sign;
  int i;
  int myone = 1;
  int * signum;
  signum = &myone;
  gsl_matrix * tmp = gsl_matrix_alloc(D,D);
  gsl_matrix * Cinv = gsl_matrix_alloc(D,D);
  gsl_vector * tmpvec = gsl_vector_alloc(D);
  gsl_vector * diff = gsl_vector_alloc(D);
  double range;
  gsl_matrix_memcpy(tmp,C);
  gsl_linalg_LU_decomp (tmp, p, signum);
  gsl_linalg_LU_invert (tmp, p, Cinv);

  for(i=0;i<D;i++){
    gsl_vector_set(diff,i,gsl_vector_get(coor,i)-gsl_vector_get(center,i));
  }

  gsl_blas_dsymv (CblasUpper, 1.0, Cinv, diff, 0.0, tmpvec);
  gsl_blas_ddot (diff, tmpvec, &range);

  gsl_permutation_free(p);
  gsl_matrix_free(tmp);
  gsl_matrix_free(Cinv);
  gsl_vector_free(tmpvec);
  gsl_vector_free(diff);

  if (range > f) return 0;
  else return 1;
}


double GetRandomCovMat(int D, gsl_matrix * C) {

  int i,j;
  float coor[D];
  float tmp;
  double specrad = 0;

  // generate random covariance matrix
  gsl_matrix * X = gsl_matrix_alloc(D,D);
  gsl_matrix * Diag = gsl_matrix_calloc(D,D);
  gsl_matrix * tmpmat = gsl_matrix_calloc(D,D);
  // make random diagonal matrix  
  for(i=0;i<D;i++){
    tmp = uniform(0.5,3);
    if (tmp>specrad) specrad = tmp;
    gsl_matrix_set(Diag,i,i,tmp);
  }

  //generate random orthogonal matrix by Gram-Schmidt
  float myX[3][3];
  //set random matrix
  for(i=0;i<D;i++){
    for(j=0;j<D;j++){
      myX[i][j] = uniform(-1,1);
    } 
  }
  //normalize first column
  tmp = sqrt(pow(myX[0][0],2)+pow(myX[1][0],2)+pow(myX[2][0],2));
  myX[0][0] /= tmp;
  myX[1][0] /= tmp;
  myX[2][0] /= tmp;
  // project second column
  tmp = myX[0][0]*myX[0][1]+ myX[1][0]*myX[1][1]+ myX[2][0]*myX[2][1];
  myX[0][1] -= tmp*myX[0][0];
  myX[1][1] -= tmp*myX[1][0];
  myX[2][1] -= tmp*myX[2][0];
  //normalize second column
  tmp = sqrt(pow(myX[0][1],2)+pow(myX[1][1],2)+pow(myX[2][1],2));
  myX[0][1] /= tmp;
  myX[1][1] /= tmp;
  myX[2][1] /= tmp;
  // project third column
  tmp = myX[0][0]*myX[0][2]+ myX[1][0]*myX[1][2]+ myX[2][0]*myX[2][2];
  myX[0][2] -= tmp*myX[0][0];
  myX[1][2] -= tmp*myX[1][0];
  myX[2][2] -= tmp*myX[2][0];
  tmp = myX[0][1]*myX[0][2]+ myX[1][1]*myX[1][2]+ myX[2][1]*myX[2][2];
  myX[0][2] -= tmp*myX[0][1];
  myX[1][2] -= tmp*myX[1][1];
  myX[2][2] -= tmp*myX[2][1];
  //normalize third column
  tmp = sqrt(pow(myX[0][2],2)+pow(myX[1][2],2)+pow(myX[2][2],2));
  myX[0][2] /= tmp;
  myX[1][2] /= tmp;
  myX[2][2] /= tmp;
  // copy to X
  for(i=0;i<D;i++){
    for(j=0;j<D;j++){
      gsl_matrix_set(X,i,j,myX[i][j]);
    } 
  }
  // C = X^T.D.X
  gsl_blas_dgemm (CblasTrans, CblasNoTrans, 1.0, X, Diag, 0.0, tmpmat);
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, tmpmat,X, 0.0, C);
  
  gsl_matrix_free(X);
  gsl_matrix_free(Diag);
  gsl_matrix_free(tmpmat);

}


double GetRandomEllipsoid(int D, gsl_vector * center, gsl_matrix * C, double * f) {
  int i;
  double specrad;
  // set center
  for(i=0;i<D;i++){
    gsl_vector_set(center,i,uniform(-10,10));
  }
  // set covariance matrix
  specrad = GetRandomCovMat(D,C);
  //generate random enlargement factor
  *f = uniform(0.1,3);
};

int TestSampleEllipsoid(){
  // randomly creates ellipsoid data, samples 100 times, puts out coordinates and ellipsoid data in format understood by testsampleellipsoid.nb
  srand(time(NULL));

  int D = 3;
  int Nsample = 100;
  int i,j;
  gsl_matrix * C = gsl_matrix_alloc(D,D);
  gsl_vector * center = gsl_vector_alloc(D);
  double f;

  GetRandomEllipsoid(D,center,C,&f);

  //sample this ellipsoid Nsample times
  gsl_vector * coorvec[Nsample];
  for(i=0;i<Nsample;i++){
    coorvec[i] = gsl_vector_alloc(D);
  }

  for(i=0;i<Nsample;i++){
    SampleEllipsoid(D, center, C, f, coorvec[i]);
    for (j=0;j<D;j++){
      printf("%f ",gsl_vector_get(coorvec[i],j));    
    } 
  printf("\n");
  }
  // print ellipsoid data

  printf("%f\n",f);

  printf("\n"); 
  for(j=0;j<D;j++){
    printf("%f ",gsl_vector_get(center,j));
  }
  printf("\n"); 
  for(i=0;i<D;i++){
    for(j=0;j<D;j++){
      printf("%f ",gsl_matrix_get(C,i,j));
    }
    printf("\n"); 
  }

  gsl_matrix_free(C);
  gsl_vector_free(center);
  for(i=0;i<Nsample;i++){
    gsl_vector_free(coorvec[i]);
  }

  return 0;
}




int TestIsMember(){
  // generate random ellipsoid, sample points uniformly from cube surrounding ellipsoid. print out coordinates with IsMember-flag and ellipsoid data in format for testismember.nb

  int i,j;
  int D = 3;
  int Nsamples = 1000;
  float coor[D];
  float tmp;
  float specrad = 0;

  gsl_matrix * C = gsl_matrix_alloc(D,D);
  gsl_vector * center = gsl_vector_alloc(D);
  double f;

  specrad = GetRandomEllipsoid(D, center, C, &f);

  //sample ellipsoid-surrounding box Nsamples times
  gsl_vector * coorvec[Nsamples];

  for(i=0;i<Nsamples;i++){
    coorvec[i] = gsl_vector_alloc(D);
    for (j=0;j<D;j++){
      tmp = uniform(gsl_vector_get(center,j)-specrad,gsl_vector_get(center,j)+specrad);
      printf("%f ",tmp);    
      gsl_vector_set(coorvec[i],j,tmp);
    } 
    printf("%i\n",IsMember(D,coorvec[i],center,C,f));
  }
  // print ellipsoid data

  printf("%f\n",f);

  for(j=0;j<D;j++){
    printf("%f ",gsl_vector_get(center,j));
  }
  printf("\n"); 
  for(i=0;i<D;i++){
    for(j=0;j<D;j++){
      printf("%f ",gsl_matrix_get(C,i,j));
    }
    printf("\n"); 
  }

  return 0;

}
