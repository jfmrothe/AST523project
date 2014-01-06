#include <ellipsoid.hh>
#include <cluster.hh>
#include <sampleellipsoid.hh>

#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>

float dist(int D, float * pos1, float * pos2){
  // returns 2-norm of difference between two given vectors
  float val=0;
  for (int i=0;i<D;i++){
    val += pow(pos1[i]-pos2[i],2);
  }
  return sqrt(val);
}

int indexofSmallestElement(float * array, int size)
// returns the index of the smallest element from a given array
{
  int index = 0;

  if (size != 1){
    int n = array[0];
    for (int i = 1; i < size; i++){
      if (array[i] < n){
	n = array[i];
        index = i;
      }
    }
  }
return index;
}

int KMeans(int D, int N, int k, float * points, int * grouping)
// Implementation of the K-means algorithm
// following MacKay, D. C. J., 2003, Information Theory, Inference and Learning Algorithms, Cambridge University Press, Cambridge, p. 640
// points gives the coordinates of N points in the D-dimensional [0,1]-hypercube, to be split into k clusters. Association of each point to a cluster is returned in int-array grouping
{
  int i,j;

  float centers[k][D];
  //  float ** centers = new float * [k];
  //  for(i=0;i<k;i++){
  //    centers[i] = new float [D];
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
      for(j=0;j<D;j++){
	centers[i][j] = 0;
      }
    }
    for(i=0;i<N;i++){
      nomemb[grouping[i]]++;
      for(j=0;j<D;j++){
	centers[grouping[i]][j] += points[D*i+j];
      }
    }
        for(i=0;i<k;i++){
          for(j=0;j<D;j++){
	    centers[i][j] = centers[i][j] / nomemb[i];
          }
        }
  // assign points to closest group center
    changed = false;
    for(i=0;i<N;i++){
      current = grouping[i];
      float mindist = dist(D,&points[D*i],&centers[current][0]);
      int mink = current;
      for(j=0;j<k;j++){
	if( dist(D,&points[D*i],&centers[j][0])<mindist ){
	  mindist = dist(D,&points[D*i],&centers[j][0]);
	  mink = j;
	}
      }
      if( mink != current ){
	changed = true;
	grouping[i] = mink;
      }
    }
  }
  // de-allocation
  //  for(i=0;i<k;i++){
  //  delete [] centers[i];
  //}
  //delete [] centers;

  return 0;
}

int TestKMeans()
{
  // Clusters N points from D-dimensional [0,1]-hypercube into k groups, prints list of groups for testkmeans.nb
  srand(time(NULL));

  int D = 3;
  int N = 1000;
  int k = 5;
  float points[N*D];
  int grouping[N];
  int i,j,m;

  for(i=0;i<N*D;i++){
    points[i] = uniform(0,1);
  }
 
  KMeans(D,N,k,&points[0],&grouping[0]);

  for(i=0;i<k;i++){
    printf("\n");
    for(j=0;j<N;j++){
      if(grouping[j]==i){
	for(m=D*j;m<D*j+D;m++){
  	  printf("%f\t",points[m]);
	}
	printf("\n");
      }
    }
  }
  return 0;
}


double ClusteringQuality(int D, int N, int Niter, Ellipsoid * clustering, int Nell) {
  // computes the ratio of remaining prior volume to the sum of all ellipsoid volumes given in clustering (adress of vector of ellipsoid objects, length Nell)
  // no accounting for overlapping ellipsoids -> overlap means bad quality.
  int i;
  float pi = 2*acos(0);
  double ellvol = 0;

  for (i=0;i<Nell;i++){
    ellvol += clustering[i].GetVol();
  }

  return ellvol / exp(-1.0*((double)Niter)/((double)N));
}

int TestEllipsoidalPartitioning(){ 
  // Sample 1000 points from torus in z=0-plane
  double r[3],ra[3],rc[3];
  double R0 = 0.05;
  double R1 = 0.3;
  rc[0] = 0.5;
  rc[1] = 0.5;
  rc[2] = 0.5;
  int N = 1000;
  int D = 3;
  int i = 0;
  int j;
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

  // free all stuff
  for(i=0;i<N;i++){
    gsl_vector_free(coors[i]);
  }
  return 0;
}
