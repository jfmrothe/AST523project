%module Ellipsoid
%{
#define SWIG_FILE_WITH_INIT
%}

%include "numpy.i"
%init %{
import_array();
%}

%{
#include "Ellipsoid.h"
%}

%include "Ellipsoid.h"
