#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
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
  
}

int IsMember(int D, gsl_vector * coor,  gsl_vector * center, gsl_matrix * C, double f){
  // returns 1 if point specified by coor lies within ellipsoid specified by center, C, f. returns 0 if not
  gsl_permutation *p = gsl_permutation_calloc(D);
  //  int sign = +1;
  //  int * signum = &sign;
  int i;
  int * signum;
  *signum = 1;
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

  if (range > f) return 0;
  else return 1;
}

void FindEnclosingEllipsoid(int D, int N, gsl_vector ** coors,  gsl_vector * center, gsl_matrix * C, double * f){
  // returns enclosing ellipsoid data for a given point cloud (N D-dimensional coordinates)
  int i,j;
  double tmp;
  gsl_vector * tmpvec = gsl_vector_alloc(D);
  gsl_vector * tmpvec2 = gsl_vector_alloc(D);
   
  FILE * fp = fopen("debugout.txt","w");
  
  //print input coordinates -> works
  for(i=0;i<N;i++){
    for(j=0;j<D;j++){
      fprintf(fp,"%f ",gsl_vector_get(coors[i],j));
    }
  fprintf(fp,"\n");
  }


  fprintf(fp,"center:\n");
  //find center -> works
  for(i=0;i<D;i++){
    tmp = 0.0;
    for(j=0;j<N;j++){
      tmp += gsl_vector_get(coors[j],i);
    }
    gsl_vector_set(center,i,tmp/N);
    fprintf(fp,"%f \n",gsl_vector_get(center,i));
  }

  //find covariance matrix
  // set C to zero
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 0.0, C, C, 0.0, C);
  //add covariance
  for(j=0;j<N;j++){
    for(i=0;i<D;i++){
      // subtract center from each coor to tmpvec
      gsl_vector_set(tmpvec,i,gsl_vector_get(coors[j],i)-gsl_vector_get(center,i));
    }
    // symmetric rank-1 update
    gsl_blas_dsyr(CblasUpper, 1.0/N , tmpvec, C);
  }
  // copy upper half into lower half
  for(i=0;i<D;i++){
    for(j=i+1;j<D;j++){
      gsl_matrix_set(C,j,i,gsl_matrix_get(C,i,j));
    }
  }

  // output of C -> works
  fprintf(fp,"\nC:\n");
  for(i=0;i<D;i++){
    for(j=0;j<D;j++){
      fprintf(fp,"%i,%i: %f\n",i,j,gsl_matrix_get(C,i,j));
    } 
  }

  //invert covariance matrix
  int signum = +1;
  gsl_matrix * tmpmat = gsl_matrix_alloc(D,D);
  gsl_matrix * Cinv = gsl_matrix_alloc(D,D);
  gsl_permutation *p = gsl_permutation_calloc(D);
  //  gsl_vector * tmpvec = gsl_vector_alloc(D);

  gsl_matrix_memcpy(tmpmat,C);
  gsl_linalg_LU_decomp (tmpmat, p, &signum);
  gsl_linalg_LU_invert (tmpmat, p, Cinv);

  // output of Cinv -> works
  fprintf(fp,"\nCinv:\n");
  for(i=0;i<D;i++){
    for(j=0;j<D;j++){
      fprintf(fp,"%i,%i: %f\n",i,j,gsl_matrix_get(Cinv,i,j));
    } 
  }
  //find enlargement factor
  *f = 0;
  for(j=0;j<N;j++){
    for(i=0;i<D;i++){
      // subtract center from each coor to tmpvec
      gsl_vector_set(tmpvec2,i,gsl_vector_get(coors[j],i)-gsl_vector_get(center,i));
    }
  gsl_blas_dsymv (CblasUpper, 1.0, Cinv, tmpvec2, 0.0, tmpvec);
  gsl_blas_ddot (tmpvec2, tmpvec, &tmp);
  if (tmp > *f) *f = tmp;
  }
  fprintf(fp,"\nf: %f\n",*f);
  fclose(fp);
}

int main(){

  srand(time(NULL));

  int i,j;
  int D=3;
  float coor[D];
  float tmp;
  // initialize covariance matrix
  gsl_matrix * C = gsl_matrix_alloc(D,D);
  //  float myC[3][3] = {{1.5,0.5,0},{0.5,1.5,0},{0,0,3}};
  //  float myC[3][3] = {{1.5,0,0.5},{0,3,0},{0.5,0,1.5}};
  //  float myC[3][3] = {{3,0,0},{0,1.5,0.5},{0,0.5,1.5}};
  float myC[3][3] = {{9,1,4},{1,9,1},{4,1,9}};
  for(i=0;i<D;i++){
    for(j=0;j<D;j++){
      gsl_matrix_set(C,i,j,myC[i][j]);
    } 
  }

  
  gsl_vector * center = gsl_vector_alloc(D);
  for(i=0;i<D;i++){
    gsl_vector_set(center,i,10.0*i);
  }

  double f = 1.0;

  gsl_vector * coorvec[100];
  for(i=0;i<100;i++){
    coorvec[i] = gsl_vector_alloc(D);
  }

  for(i=0;i<100;i++){
    SampleEllipsoid(D, center, C, f, coorvec[i]);
    for (j=0;j<D;j++){
      printf("%f ",gsl_vector_get(coorvec[i],j));    
    } 
    //    printf("%i\n",IsMember(D, coorvec[i], center, C, f));
  printf("\n");
  }

  gsl_matrix * myCC = gsl_matrix_alloc(D,D);
  gsl_vector * mycenter = gsl_vector_alloc(D);
  double myf;
  
  FindEnclosingEllipsoid(D, 100, &coorvec[0], mycenter, myCC, &myf);


  gsl_vector * testcoor = gsl_vector_alloc(D);
  for(i=0;i<10000;i++){
    SampleEllipsoid(D, mycenter, myCC, myf, testcoor);  
    // output
    for (j=0;j<D;j++){
      printf("%f ",gsl_vector_get(testcoor,j));    
    }
    printf("%i\n",IsMember(D, testcoor, mycenter, myCC, myf));
    //    printf("%i\n",IsMember(D, coorvec, center, C, f));
  }


  return 0;
}











int oldtestcase(){
 srand(time(NULL));

  FILE * fp = fopen("debugout.txt","w");

  int i,j;
  int D=3;
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
    gsl_vector_set(center,i,0.0*i);
    fprintf(fp,"%i: %f\n",i,gsl_vector_get(center,i));
  }
  fclose(fp);

  double f = 1.0;

  gsl_vector * coorvec = gsl_vector_alloc(D);

  for(i=0;i<10000;i++){
    for (j=0;j<D;j++){
      gsl_vector_set(coorvec,j,uniform(-3,3));    
    } 
    //SampleEllipsoid(D, center, C, f, coorvec);
    // output
    for (j=0;j<D;j++){
      printf("%f ",gsl_vector_get(coorvec,j));    
    } 
    printf("%i\n",IsMember(D, coorvec, center, C, f));
  }


  return 0;
}
