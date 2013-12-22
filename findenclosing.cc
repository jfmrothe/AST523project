# include "findenclosing.hh"
# include "sampleellipsoid.hh"

void FindEnclosingEllipsoid(int D, int N, gsl_vector ** coors,  gsl_vector * center, gsl_matrix * C, double * f){
  // returns enclosing ellipsoid data for a given point cloud (N D-dimensional coordinates), enlargement factor f is chosen so that all points are below the ellipsoid surface
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

int TestFindEnclosingEllipsoid(){
  // takes 5 random samples from the unisphere, returns these and data of the enclosing ellipsoid, in a format understood by testfindenclosing.nb
  int myN=5;
  int myD=3;
  int i,j;
  float coor[myD];
  gsl_vector * samples[myN];

  for(i=0;i<myN;++i){
    unisphere(&coor[0],myD);
    samples[i] = gsl_vector_alloc(myD); 
    for(j=0;j<myD;j++){
      printf("%f ",coor[j]);
      gsl_vector_set(samples[i],j,coor[j]);
    }
    printf("\n");
  }
  
  gsl_vector * mycenter = gsl_vector_alloc(myD);
  gsl_matrix * C = gsl_matrix_calloc(myD,myD);
  double f;

  FindEnclosingEllipsoid(myD, myN, &samples[0], mycenter, C, &f);

  printf("\n"); 
  for(j=0;j<myD;j++){
    printf("%f ",gsl_vector_get(mycenter,j));
  }
  printf("\n"); 
  for(i=0;i<myD;i++){
    for(j=0;j<myD;j++){
      printf("%f ",gsl_matrix_get(C,i,j));
    }
    printf("\n"); 
  }
  printf("%f\n",f);

  return 0;
}
