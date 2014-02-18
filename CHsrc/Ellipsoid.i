%module Ellipsoid
%{
#define SWIG_FILE_WITH_INIT
%}

%include "numpy.i"
%init %{
import_array();
%}

%{
#include <list>
#include <fstream>
#include <string>
#include <float.h>
#include <math.h>
#include <vector>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include "Point.h"
#include "Ellipsoid.h"
%}

%include "Ellipsoid.h"
