# ifndef _AST523_GSL_LIB_
# define _AST523_GSL_LIB_
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
# endif

# ifndef _AST523_FIND_ENCLOSING_
# define _AST523_FIND_ENCLOSING_
void FindEnclosingEllipsoid(int D, int N, gsl_vector ** coors,  gsl_vector * center, gsl_matrix * C, double * f);
int TestFindEnclosingEllipsoid();
# endif
