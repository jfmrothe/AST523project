# ifndef _AST523_GSL_LIB_
# define _AST523_GSL_LIB_
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
# endif

# ifndef _AST523_SAMPLE_ELLIPSOID_
# define _AST523_SAMPLE_ELLIPSOID_

float uniform(float min, float max);
float boxmuller();
float quadr();
void unisphere(float * coor, int D);
void SampleEllipsoid(Ellipsoid &ell, gsl_vector * coor, int D);
int IsMember(Ellipsoid &ell, gsl_vector * coor, int D);
double GetRandomCovMat(int D, gsl_matrix * C);
double GetRandomEllipsoid(int D, gsl_vector * center, gsl_matrix * C, double * f);
int TestSampleEllipsoid();
int TestIsMember();
# endif
