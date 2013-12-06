#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>

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

void sampleellipsoid(int D, gsl_vector * center, gsl_matrix * C, double f, gsl_vector * coor){
  //returns pseudorandom number coor uniformly distributed in ellipsoid given by the center vector center, covariance matrix C and enlargement factor f, so that x^T(fC)^-1x<=1
  int i,j;

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
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, f, X, Dprime, 0.0, T);


  // final coordinates
  gsl_vector_memcpy(coor,center);
  gsl_blas_dgemv(CblasNoTrans, 1.0, T, spheresample, 1.0, coor);
  
}

int main(){

  srand(time(NULL));

  FILE * fp = fopen("debugout.txt","w");

  int i,j;
  int D=3;
  float coor[D];
  float tmp;
  // initialize covariance matrix
  fprintf(fp,"C:\n");
  gsl_matrix * C = gsl_matrix_alloc(D,D);
  //  float myC[3][3] = {{1.5,0.5,0},{0.5,1.5,0},{0,0,3}};
  //  float myC[3][3] = {{1.5,0,0.5},{0,3,0},{0.5,0,1.5}};
  //  float myC[3][3] = {{3,0,0},{0,1.5,0.5},{0,0.5,1.5}};
  float myC[3][3] = {{9,1,4},{1,9,1},{4,1,9}};
  for(i=0;i<D;i++){
    for(j=0;j<D;j++){
      gsl_matrix_set(C,i,j,myC[i][j]);
      fprintf(fp,"%i,%i: %f\n",i,j,gsl_matrix_get(C,i,j));
    } 
  }

  
  fprintf(fp,"center:\n");
  gsl_vector * center = gsl_vector_alloc(D);
  for(i=0;i<D;i++){
    gsl_vector_set(center,i,10*i);
    fprintf(fp,"%i: %f\n",i,gsl_vector_get(center,i));
  }
  fclose(fp);

  double f = 1.0;

  gsl_vector * coorvec = gsl_vector_alloc(D);

  for(i=0;i<10000;i++){
    sampleellipsoid(D, center, C, f, coorvec);
    // output
    for (j=0;j<D;j++){
      printf("%f ",gsl_vector_get(coorvec,j));    
    } 
    printf("\n");
  }

  return 0;
}
