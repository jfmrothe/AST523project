%module Samplers
%{
#define SWIG_FILE_WITH_INIT
%}
%include "std_string.i"
%{
#define SWIG_FILE_WITH_INIT
%}
%include "numpy.i"
%init %{
import_array();
%}

%apply(double * INPLACE_ARRAY1, int DIM1){(double* min_vals, int nmin), (double* max_vals, int nmax)};
%apply(double * INPLACE_ARRAY2, int DIM1, int DIM2){(double *Alltheta, int nx, int ny)};
%apply(double * INPLACE_ARRAY1, int DIM1){(double * logL,int nL)}; 
%{
#include "Samplers.h"
#include "string"
%}

%include "Samplers.h"
