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
%apply(double * INPLACE_ARRAY2, int DIM1, int DIM2){(double * posterior, int nx, int ny)};
%apply(double * INPLACE_ARRAY1, int DIM1){(double * prob, int np)};
%apply(double * INPLACE_ARRAY1, int DIM1){(double * logL,int nL)}; 
%apply(double * INPLACE_ARRAY1, int DIM1){(double * theta,int nt)}; 
%apply(double * INPLACE_ARRAY1, int DIM1){(double * logzinfo,int nz)};

%{
#include "Samplers.h"
#include "Ellipsoid.h"
#include "string"
%}

%include "Ellipsoid.h"
%include "Samplers.h"
