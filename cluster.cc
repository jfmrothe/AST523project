#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>

float uniform(float min, float max){
  //returns pseudorandom number from uniform distribution
  return (((float) rand())/((float) RAND_MAX))*(max-min) + min;
}

float dist(int D, float * pos1, float * pos2){
  float val=0;
  for (int i=0;i<D;i++){
    val += pow(pos1[i]-pos2[i],2);
  }
  return sqrt(val);
}

int indexofSmallestElement(float * array, int size)
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

int k_means(int D, int N, int k, float * points, int * grouping)
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

int main()
{
  srand(time(NULL));

  int D = 3;
  int N = 20;
  int k = 3;
  float points[N*D];
  int grouping[N];
  int i,j,m;

  for(i=0;i<N*D;i++){
    points[i] = uniform(0,1);
  }
 
  k_means(D,N,k,&points[0],&grouping[0]);

  for(i=0;i<k;i++){
    printf("Group %i\n",i);
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
