#include <ellipsoid.hh>

#include <vector>

#include <cluster.hh>
#include <sampleellipsoid.hh>
#include <findenclosing.hh>

#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>



double dist(int D, gsl_vector * pos1, gsl_vector * pos2){
  // returns 2-norm of difference between two given vectors
  double val=0;
  for (int i=0;i<D;i++){
    val += pow(gsl_vector_get(pos1,i)-gsl_vector_get(pos2,i),2);
  }
  return sqrt(val);
}

int KMeans(gsl_vector ** points, int D, int N, int k, int * grouping, int * groupsize)
// Implementation of the K-means algorithm
// following MacKay, D. C. J., 2003, Information Theory, Inference and Learning Algorithms, Cambridge University Press, Cambridge, p. 640
// points gives the coordinates of N points in the D-dimensional [0,1]-hypercube, to be split into k clusters. Association of each point to a cluster is returned in int-array grouping
{
  int i,j;

  gsl_vector * centers[k];
  //  double centers[k][D];
  //  double ** centers = new double * [k];
  //  for(i=0;i<k;i++){
  //    centers[i] = new double [D];
  //  }



  int nomemb[k];
  bool changed = true;
  int current = 0; 
  // assign points randomly to groups
  for(i=0;i<N;i++){
    grouping[i] = rand()%k;
  }

  while(changed){
  // calculate group centers
    for(i=0;i<k;i++){
      nomemb[i] = 0;
      centers[i] = gsl_vector_calloc(D);
    }
    for(i=0;i<N;i++){
      nomemb[grouping[i]]++;
      for(j=0;j<D;j++){
	gsl_vector_set(centers[grouping[i]],j, gsl_vector_get(centers[grouping[i]],j)+gsl_vector_get(points[i],j));
      }
    }
    for(i=0;i<k;i++){
      for(j=0;j<D;j++){
	gsl_vector_set(centers[i],j,gsl_vector_get(centers[i],j)/nomemb[i]);
      }
    }
  // assign points to closest group center
    changed = false;
    for(i=0;i<N;i++){
      current = grouping[i];
      double mindist = dist(D, points[i], centers[current]);
      int mink = current;
      for(j=0;j<k;j++){
	if( dist(D, points[i], centers[j])<mindist ){
	  mindist = dist(D, points[i], centers[j]);
	  mink = j;
	}
      }
      if( mink != current ){
	changed = true;
	grouping[i] = mink;
      }
    }
  }
  //count group sizes
  for(i=0;i<k;i++){
    groupsize[i]=0;
  }
  for(i=0;i<N;i++){
    groupsize[grouping[i]] += 1;
  }
  return 0;
}

int TestKMeans()
{
  // Clusters N points from D-dimensional [0,1]-hypercube into k groups, prints list of groups for testkmeans.nb
  srand(time(NULL));

  int D = 3;
  int N = 1000;
  int k = 9;
  gsl_vector * points[N];
  int grouping[N];
  int i,j,m;

  for(i=0;i<N;i++){
    points[i] = gsl_vector_alloc(D);
    for(j=0;j<D;j++) {
      gsl_vector_set(points[i],j,uniform(0,1));
    }
  }
  int groupsize[k];
  KMeans(&points[0],D,N,k,&grouping[0],&groupsize[k]);

  for(i=0;i<k;i++){
    printf("\n");
    for(j=0;j<N;j++){
      if(grouping[j]==i){
	for(m=0;m<D;m++){
  	  printf("%f\t",gsl_vector_get(points[j],m));
	}
	printf("\n");
      }
    }
  }
  return 0;
}


double ClusteringQuality(std::vector<Ellipsoid*>& clustering, double Xtot) {
  // computes the ratio of remaining prior volume Xtot to the sum of all ellipsoid volumes given in clustering (vector of pointers to ellipsoid objects)
  // no accounting for overlapping ellipsoids -> overlap means bad quality.
  int i;
  double ellvol = 0;

  for (i=0;i<clustering.size();i++){
    ellvol += clustering[i]->getVol();
  }

  return ellvol / Xtot;
}

int TestClusteringQuality() {
  // make some spheres and call ClusteringQuality with Niter = 0
  int i;
  int D = 3;
  int N = 1000;
  int Niter = 0;
  int Nell = 4;
  double rad[Nell];
  for (i=0;i<Nell;i++) {
    rad[i] = i;
  }

  gsl_vector * center = gsl_vector_calloc(D);
  gsl_matrix * C = gsl_matrix_calloc(D,D);
  for(i=0;i<D;i++){
    gsl_matrix_set(C,i,i,1.0);
  }

  std::vector<Ellipsoid*> clustering;
  clustering.push_back(new Ellipsoid(D,center,C,pow(rad[0],2)));
  clustering.push_back(new Ellipsoid(D,center,C,pow(rad[1],2)));
  clustering.push_back(new Ellipsoid(D,center,C,pow(rad[2],2)));
  clustering.push_back(new Ellipsoid(D,center,C,pow(rad[3],2)));
  
  printf("%f\n",ClusteringQuality(clustering, 1));

  for(i=0;i++;i<clustering.size()){
    delete clustering[i];
  }
  gsl_vector_free(center);
  gsl_matrix_free(C);
  return 0;
}


void EllipsoidalPartitioning(gsl_vector ** coors, int D, int N, double Xtot, std::vector<Ellipsoid*>& clustering) {
  // this is the variant currently in development, because segfaults could be avoided. will check if this implementation is "good" in some sense...
  // vector of ellipsoid pointers is returned, these are deleted after use in the function calling EllipsoidalPartitioningVec

  // performs Algorithm I from Feroz, Hobson and Bridges (2009) on N points in D-dimensional [0,1]-hypercube given in coors. Returns resulting number of ellipsoids in Nell, uses new to create array of ellipsoids, first address is returned in clustering

  /////////////////
  // ALGORITHM I //
  /////////////////

  // create mainEllipsoid which may be split further
  // ?? in recursion depth, this first ellipsoid can be passed by parent??
  Ellipsoid mainEll = FindEnclosingEllipsoid(coors,D,N);

  //enlarge if neccessary
  if(Xtot>mainEll.getVol()) {
    mainEll.setEnlFac(mainEll.getEnlFac()*pow(Xtot/mainEll.getVol(),1.0/D));
  }
  // initialize splitting using kmeans
  int i;
  int grouping[N];
  int k=2;
  int groupsize[k];
  KMeans(&coors[0],D,N,k,&grouping[0],&groupsize[0]);

  // return main ellipsoid immediately if partitioning would create singular (flat) ellipsoid
  if(groupsize[0]<D+1 or groupsize[1]<D+1) {
    clustering.push_back (new Ellipsoid(D, mainEll.getCenter(), mainEll.getCovMat(), mainEll.getEnlFac()) );
    return;
  }

  bool changed = true;
  double X1;
  double X2;
  double h1;
  double h2;
  double tmp1,tmp2;
  double vol1,vol2;

  while(changed) {

    // find new sub-Ellipsoids
    //"locality" of these variables removes object overwriting trouble
    Ellipsoid subEll1 = FindSelectedEnclosingEllipsoid(&coors[0],D,N,&grouping[0],0);
    vol1 = subEll1.getVol();
    X1 = ((double)groupsize[0]/N)*Xtot;
    Ellipsoid subEll2 = FindSelectedEnclosingEllipsoid(&coors[0],D,N,&grouping[0],1);
    vol2 = subEll2.getVol();
    X2 = ((double)groupsize[1]/N)*Xtot;

    //enlarge if neccessary
    if(X1>vol1) {
      subEll1.setEnlFac(subEll1.getEnlFac()*pow(X1/vol1,(double)(1.0/D)));
      vol1 = subEll1.getVol();
    }
    if(X2>vol2) {
    subEll2.setEnlFac(subEll2.getEnlFac()*pow(X2/vol2,(double)(1.0/D)));
    vol2 = subEll2.getVol();
    }

    // optimize splitting using mahalanobis-h-metric
    changed = false;    
    for(i=0;i<N;i++) {

      tmp1 = subEll1.mdist(coors[i]);
      tmp2 = subEll2.mdist(coors[i]);
   
      h1 = vol1 / X1 * tmp1;
      h2 = vol2 / X2 * tmp2;

      if (h2<h1 and grouping[i]==0) {
	changed = true;
	grouping[i] = 1;
	groupsize[0] -= 1;
	groupsize[1] += 1;
      }
      if (h1<h2 and grouping[i]==1) {
	changed = true;
	grouping[i] = 0;
	groupsize[0] += 1;
	groupsize[1] -= 1;
      }
    }
  }
  // judgement if splitting should be continued
  if( vol1+vol2<mainEll.getVol() or mainEll.getVol()>2*Xtot) {

    // recursively start the splittings of subEll1 and subEll2
    gsl_vector * subcoors1[groupsize[0]];
    for(i=0;i<groupsize[0];i++){
      subcoors1[i] = gsl_vector_alloc(D);
    }
    SelectFromGrouping(coors, D, N, grouping, 0, &subcoors1[0]);
    EllipsoidalPartitioning( &subcoors1[0], D, groupsize[0], X1, clustering);
    // second one
    gsl_vector * subcoors2[groupsize[1]];
    for(i=0;i<groupsize[1];i++){
      subcoors2[i] = gsl_vector_alloc(D);
    }
    SelectFromGrouping(coors, D, N, grouping, 1, &subcoors2[0]);
    EllipsoidalPartitioning( &subcoors2[0], D, groupsize[1], X2, clustering);

  }
  else{
    // allocate memory for mainEll and push it to the vector
    clustering.push_back (new Ellipsoid(D, mainEll.getCenter(), mainEll.getCovMat(), mainEll.getEnlFac()) );
  }
  return;
}


int TestEllipsoidalPartitioning(){ 
  // Sample N points from torus in z=0-plane
  double pi = 2*acos(0);

  double r[3],ra[3],rc[3];
  double R0 = 0.05;
  double R1 = 0.3;
  rc[0] = 0.5;
  rc[1] = 0.5;
  rc[2] = 0.5;
  int N = 1000;
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
  
  std::vector<Ellipsoid*> clustering;
  // inherent volume of data cloud: torus
  double Xtot = pow(R0,2)*pi*2*pi*R1;
  EllipsoidalPartitioning(&coors[0], D, N, Xtot, clustering);
  int Nell = clustering.size();
  // output
  gsl_vector * mycenter;
  gsl_matrix * C;
  double f;
  for(i=0;i<Nell;i++){
    (*clustering[i]).printout();
  }
  printf("%f\n",ClusteringQuality(clustering,Xtot));

  // free all stuff
  for(i=0;i<Nell;i++){
    delete clustering[i];
  }
  for(i=0;i<N;i++){
    gsl_vector_free(coors[i]);
  }
  return 0;
}
