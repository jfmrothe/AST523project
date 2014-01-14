%module Point
%{
#define SWIG_FILE_WITH_INIT
%}

%include "Point.i"
%init %{
import_array();
%}
%{
#include "Point.h"
%}

%include "Point.h"
