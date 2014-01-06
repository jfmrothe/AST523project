#include "Samplers.h"
#include "Point.h"

float boxmuller()
{
  //returns univariate, zero-mean gaussian sample
  float u1 = UNIFORM;
  float u2 = UNIFORM;
  float pi = 2*acos(0);
  while( u1 == 0.0 ){
    u1 = UNIFORM;
  }
  return sqrt(-2.0*log(u1))*cos(2*pi*u2);
}

float quadr()
{
  //returns sample from quadratic distribution between 0 and 1
  return pow(UNIFORM,1.0/3.0);
}

void unisphere(float * coor, int D)
{
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

// ************************************** Samplers methods start here

Samplers::Samplers(int Dim)
{
    D = Dim;
    center = gsl_vector_alloc(D);
    coor = gsl_vector_alloc(D);
    C = gsl_matrix_calloc(D,D);
}

Samplers::~Samplers()
{
    gsl_vector_free(center);
    gsl_vector_free(coor);
    gsl_matrix_free(C);
}

void Samplers::set_vectors_zero()
{
    gsl_vector_set_zero(center);
    gsl_vector_set_zero(coor);
}

bool Samplers::u_in_hypercube()
{
    double x_i;

    for(int i = 0; i<D; i++)
    {
        x_i = gsl_vector_get(coor, i);
        if(x_i < 0 || x_i > 1) {return false;}
    }
    return true;
}

void Samplers::FindEnclosingEllipsoid(int N, Point *pts[])
{
  /* returns enclosing ellipsoid data for a given point cloud (N D-dimensional coordinates), 
     enlargement factor f is chosen so that all points are below the ellipsoid surface */
  int i,j;
  double tmp;
  gsl_vector * tmpvec = gsl_vector_alloc(D);
  gsl_vector * tmpvec2 = gsl_vector_alloc(D);
   
  //find center -> works
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
}

void Samplers::SampleEllipsoid()
{
  /*returns pseudorandom number coor uniformly distributed in ellipsoid given by the center vector center,
    covariance matrix C and enlargement factor f, so that x^T(fC)^-1x<=1 */
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
}

void Samplers::mcmc(Point* pt, Data data_obj, double logLmin)
{
    vector<double> new_coords(D);
    double step = 0.1;
    int accept = 0;
    int reject = 0;
    Point *trial;
    trial = new Point(D);
    *trial = *pt;

    for(int j = 20; j > 0; j--)
    {
        for(int i = 0; i<D; i++)
        {
            new_coords[i] = pt->get_u(i) + step*(2.0*UNIFORM - 1.0);
            new_coords[i] -= floor(new_coords[i]);
            trial->set_u(i, new_coords[i]);
        }

        trial->transform_prior();
        data_obj.lighthouse_logL(trial);

        if(trial->get_logL() > logLmin){*pt = *trial; accept++;}
        else reject++;

        if(accept > reject) step *= exp(1.0 / accept);
        else if(accept < reject) step /= exp(1.0 / reject);
    }

    delete trial;
}
