# ifndef _AST523_GSL_LIB_
# define _AST523_GSL_LIB_
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
# endif

# ifndef _AST523_FIND_ENCLOSING_
# define _AST523_FIND_ENCLOSING_
Ellipsoid FindEnclosingEllipsoid(gsl_vector ** coors, int D, int N);
int TestFindEnclosingEllipsoid();
void SelectFromGrouping(gsl_vector ** coor, int D, int N, int * grouping, int index, gsl_vector ** group);
int TestSelectFromGrouping();
Ellipsoid FindSelectedEnclosingEllipsoid(gsl_vector ** coors, int D, int N, int * grouping, int index);
int TestFindSelectedEnclosingEllipsoid();
int Test2FindSelectedEnclosingEllipsoid();
# endif
