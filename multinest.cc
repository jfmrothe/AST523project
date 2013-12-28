<<<<<<< HEAD
#include "sampleellipsoid.hh"
#include "findenclosing.hh"
=======
#include <ellipsoid.hh>
#include <sampleellipsoid.hh>
#include <findenclosing.hh>
#include <cluster.hh>       
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>

>>>>>>> 51485a1352dbd1bcd890ed4bc324c8451cc36d01

int main(){
  srand(time(NULL));
  int D = 3;
  double f;
  gsl_vector * center = gsl_vector_alloc(D);
  gsl_matrix * C = gsl_matrix_alloc(D,D);
  GetRandomEllipsoid(D, center, C, &f);

  Ellipsoid mfell (D, center, C, f);
  printf("%f\n",mfell.GetEnlFac());
  printf("%f\n",mfell.GetVol());

  gsl_vector_free(center);
  gsl_matrix_free(C);
  return 0;
}

int TestMain() {
  return TestEllipsoidalPartitioning();
}
