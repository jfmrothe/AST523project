# include <ellipsoid.hh>
# include <findenclosing.hh>
# include <sampleellipsoid.hh>

#include <vector>
#include <math.h>
#include <cluster.hh>

Ellipsoid FindEnclosingEllipsoid(gsl_vector ** coors, int D, int N){
  // returns enclosing ellipsoid data for a given point cloud (N D-dimensional coordinates), enlargement factor f is chosen so that all points are below the ellipsoid surface
  int i,j;
  double tmp;
  gsl_vector * tmpvec = gsl_vector_alloc(D);
  gsl_vector * tmpvec2 = gsl_vector_alloc(D);

  gsl_vector * center = gsl_vector_alloc(D);
  gsl_matrix * C = gsl_matrix_alloc(D,D);
  float f;

  //find center
  for(i=0;i<D;i++){
    tmp = 0.0;
    for(j=0;j<N;j++){
      tmp += gsl_vector_get(coors[j],i);
    }
    gsl_vector_set(center,i,tmp/N);
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

  //invert covariance matrix
  int signum = +1;
  gsl_matrix * tmpmat = gsl_matrix_alloc(D,D);
  gsl_matrix * Cinv = gsl_matrix_alloc(D,D);
  gsl_permutation *p = gsl_permutation_calloc(D);
  //  gsl_vector * tmpvec = gsl_vector_alloc(D);

  gsl_matrix_memcpy(tmpmat,C);
  gsl_linalg_LU_decomp (tmpmat, p, &signum);
  gsl_linalg_LU_invert (tmpmat, p, Cinv);

  //find enlargement factor
  f = 0;
  for(j=0;j<N;j++){
    for(i=0;i<D;i++){
      // subtract center from each coor to tmpvec
      gsl_vector_set(tmpvec2,i,gsl_vector_get(coors[j],i)-gsl_vector_get(center,i));
    }
  gsl_blas_dsymv (CblasUpper, 1.0, Cinv, tmpvec2, 0.0, tmpvec);
  gsl_blas_ddot (tmpvec2, tmpvec, &tmp);
  if (tmp > f) f = tmp;
  }

  gsl_vector_free(tmpvec);
  gsl_vector_free(tmpvec2);
  gsl_matrix_free(tmpmat);
  gsl_matrix_free(Cinv);
  gsl_permutation_free(p);

  Ellipsoid encl(D, center, C, f);

  gsl_vector_free(center);
  gsl_matrix_free(C);

  return encl;
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
  
  Ellipsoid ell = FindEnclosingEllipsoid(&samples[0], myD, myN);

  // retrieve data from ellipsoid object
  gsl_vector * mycenter = ell.getCenter();
  gsl_matrix * C = ell.getCovMat();
  double f = ell.getEnlFac();
  // data print out
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

  for(i=0;i<myN;++i){
    gsl_vector_free(samples[i]);
  }
  return 0;
}

void SelectFromGrouping(gsl_vector ** coor, int D, int N, int * grouping, int index, gsl_vector ** group) {
  // copies the coor-vectors for which the grouping entry equals index into group, length of array group must be pre-arranged and all entries calloced
  int i;
  int j=0;
  for(i=0;i<N;i++) {
    if(grouping[i]==index) {
      gsl_vector_memcpy(group[j],coor[i]);
      j += 1;
    }
  }
  return;
}

int TestSelectFromGrouping() {
  // create a random list of coordinates and a random grouping, print split version
  int D = 3;
  int N = 100;
  gsl_vector * coors[N];
  int grouping[N];
  int i,j;
  double tmp;

  for(i=0;i<N;i++) {
    coors[i] = gsl_vector_alloc(D);
    for(j=0;j<D;j++) {
      tmp = uniform(0.0,1.0);
      gsl_vector_set(coors[i],j,tmp);
      printf("%f ",tmp);
    }
    grouping[i] = rand()%2;
    printf("%i\n",grouping[i]);
  }

  int index,Nsub;
  for(index=0;index<2;index++) {
    // count subsample length
    Nsub = 0;
    for(i=0;i<N;i++) {
      if (grouping[i]==index) {
	Nsub += 1;
      }
    }
    // allocate memory for subsample coordinate array
    gsl_vector * subcoor[Nsub];
    for(i=0;i<Nsub;i++) {
      subcoor[i] = gsl_vector_alloc(D);
    }
    // fill subsample coordinate array
    SelectFromGrouping(coors, D, N, &grouping[0], index, &subcoor[0]);
    // print group
    printf("Group %i, length %i\n",index,Nsub);
    for(i=0;i<Nsub;i++) {
      for(j=0;j<D;j++) {
	printf("%f ",gsl_vector_get(subcoor[i],j));
      }
      printf("%i\n",index);
    }
    // free memory
    for(i=0;i<Nsub;i++) {
      gsl_vector_free(subcoor[i]);
    }
  }
  return 0;
}

Ellipsoid FindSelectedEnclosingEllipsoid(gsl_vector ** coors, int D, int N, int * grouping, int index) {
  // finds enclosing ellipsoid only of the subsample of coor-coordinates which have a grouping entry of index
  // N and D describe the formats of the full coors set, subsample is counted in here
  int i;
  
  // count subsample length
  int Nsub = 0;
  for(i=0;i<N;i++) {
    if (grouping[i]==index) {
      Nsub += 1;
    }
  }
  // allocate memory for subsample coordinate array
  gsl_vector * subcoor[Nsub];
  for(i=0;i<Nsub;i++) {
    subcoor[i] = gsl_vector_alloc(D);
  }
  // fill subsample coordinate array
  SelectFromGrouping(coors, D, N, &grouping[0], index, &subcoor[0]);

  // find enclosing
  Ellipsoid result = FindEnclosingEllipsoid(&subcoor[0], D, Nsub);

  // free memory
  for(i=0;i<Nsub;i++) {
    gsl_vector_free(subcoor[i]);
  }

  return result;
}

int TestFindSelectedEnclosingEllipsoid() {
  // sample a few points from hypercube, group upper/lower half and return two enclosing ellipsoids
  int D = 3;
  int N = 25;

  gsl_vector * coors[N];
  int grouping[N];
  int i,j;
  double tmp;

  for(i=0;i<N;i++){
    coors[i] = gsl_vector_alloc(D);
    for(j=0;j<D;j++){
      tmp = uniform(0,1);
      gsl_vector_set(coors[i],j,tmp);
      printf("%f ",tmp);
    }
    grouping[i] = (gsl_vector_get(coors[i],2)>0.5);
    printf("%i\n",grouping[i]);
  }
  Ellipsoid first = FindSelectedEnclosingEllipsoid(&coors[0], D, N, &grouping[0], 0);
  Ellipsoid second = FindSelectedEnclosingEllipsoid(&coors[0], D, N, &grouping[0], 1);
  first.printout();
  second.printout();
  return 0;
}

int Test2FindSelectedEnclosingEllipsoid() {
  // run kmeans on torus-sample and enclose the parts
  // Sample N points from torus in z=0-plane
  double r[3],ra[3],rc[3];
  double R0 = 0.05;
  double R1 = 0.3;
  rc[0] = 0.5;
  rc[1] = 0.5;
  rc[2] = 0.5;
  int N = 100;
  int D = 3;
  int i = 0;
  int j,k;
  gsl_vector * coors[N];

  while(i<N){
    // sample from cube
    for (j=0;j<D;j++){
      r[j] = uniform(0.0,1.0);
    }
    // check if in torus
    ra[0] = rc[0] + (r[0]-rc[0])*R1/sqrt(pow(r[0]-rc[0],2)+pow(r[1]-rc[1],2));
    ra[1] = rc[1] + (r[1]-rc[1])*R1/sqrt(pow(r[0]-rc[0],2)+pow(r[1]-rc[1],2));
    ra[2] = rc[2];
    if ( pow(r[0]-ra[0],2) + pow(r[1]-ra[1],2) + pow(r[2]-ra[2],2)  < pow(R0,2) ) {
      coors[i] = gsl_vector_alloc(D);
      for(j=0;j<D;j++){
	gsl_vector_set(coors[i],j,r[j]);
	printf("%f ",r[j]);
      }
      printf("\n");
      i++;
    }
    // else just go on
  }
  printf("\n");

  k = 2;
  int grouping[N];
  int groupsize[k];
  KMeans(&coors[0],D,N,k,&grouping[0],&groupsize[0]);

  Ellipsoid first = FindSelectedEnclosingEllipsoid(&coors[0], D, N, &grouping[0], 0);
  Ellipsoid second = FindSelectedEnclosingEllipsoid(&coors[0], D, N, &grouping[0], 1);
  first.printout();
  second.printout();
  return 0;
}
