#include "Samplers.h"

double dist(int D, Point *pt, gsl_vector * pos2)
{
  // returns 2-norm of difference between two given vectors
  double val=0;
  for (int i=0;i<D;i++){
    val += pow(pt->get_u(i)-gsl_vector_get(pos2,i),2);
  }
  return sqrt(val);
}

int KMeans(vector<Point *>& pts, int D, int k, int * grouping)
{
 // Implementation of the K-means algorithm following MacKay, D. C. J., 2003,
 // Information Theory, Inference and Learning Algorithms, Cambridge University
 // Press, Cambridge, p. 640 points gives the coordinates of N points in the
 // D-dimensional [0,1]-hypercube, to be split into k clusters. Association of
 // each point to a cluster is returned in int-array grouping
  int i,j;
  int N = pts.size();

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
    gsl_vector_set(centers[grouping[i]],j, gsl_vector_get(centers[grouping[i]],j)+pts[i]->get_u(j));
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
      double mindist = dist(D, pts[i], centers[current]);
      int mink = current;
      for(j=0;j<k;j++){
    if( dist(D, pts[i], centers[j])<mindist ){
      mindist = dist(D, pts[i], centers[j]);
      mink = j;
    }
      }
      if( mink != current ){
    changed = true;
    grouping[i] = mink;
      }
    }
  }

  return 0;
}

Ellipsoid FindEnclosingEllipsoid(vector<Point *>& pts, int D)
{
  // enclosing ellipsoid data for a given point cloud (N D-dimensional
  // coordinates), enlargement factor f is chosen so that all points are below
  // the ellipsoid surface

  int i,j;
  double tmp;
  int N = pts.size();
  gsl_vector * tmpvec = gsl_vector_alloc(D);
  gsl_vector * tmpvec2 = gsl_vector_alloc(D);

  gsl_vector * center = gsl_vector_alloc(D);
  gsl_matrix * C = gsl_matrix_alloc(D,D);
  float f;

  //find center
  for(i=0;i<D;i++){
    tmp = 0.0;
    for(j=0;j<N;j++){
      tmp += pts[j]->get_u(i);
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
      gsl_vector_set(tmpvec,i,pts[j]->get_u(i)-gsl_vector_get(center,i));
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
      gsl_vector_set(tmpvec2,i,pts[j]->get_u(i)-gsl_vector_get(center,i));
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

  Ellipsoid encl(D, center, C, f, pts);

  gsl_vector_free(center);
  gsl_matrix_free(C);

  return encl;
}

void SelectFromGrouping(vector<Point *>& pts, int D, int * grouping, int index, vector<Point *>& pts_subset) 
{
  // copies the coor-vectors for which the grouping entry equals index into group, 
  // length of array group must be pre-arranged and all entries calloced
  int i;
  for(i=0;i<pts.size();i++) {
    if(grouping[i]==index) {
      pts_subset.push_back(pts[i]);
    }
  }
}

bool u_in_hypercube(gsl_vector * newcoor, int D)
{
    double x_i;

    for(int i = 0; i<D; i++)
    {
        x_i = gsl_vector_get(newcoor, i);
        if(x_i < 0 || x_i > 1) {return false;}
    }
    return true;
}

