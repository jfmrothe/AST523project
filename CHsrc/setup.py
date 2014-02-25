#! /usr/bin/env python

from distutils.core import *
from distutils import sysconfig
from distutils.extension import Extension

import numpy

numpy_include = numpy.get_include()
_Point = Extension("_Point",
        ["Point_wrap.cxx",
            "Point.cc"],
        include_dirs = [numpy_include],
        )
_Ellipsoid = Extension("_Ellipsoid",
        ["Ellipsoid_wrap.cxx",
            "Ellipsoid.cc"],
        include_dirs = [numpy_include],
        )
_Samplers = Extension("_Samplers",
        ["Samplers_wrap.cxx",
            "Samplers.cc"],
        include_dirs = [numpy_include],
        )

setup(name="MutiNest",
        description = "Bayesian tools for model fitting",
        author = "Johnny Greco, Johannes Rothe for the c++ part, Xu Huang for the python part and modeling",
        author_email = "",
        url = "",
        version = "0.0.0",
        py_modules = ["Point","Ellipsoid, Samplers"],
        ext_modules = [_Point,_Ellipsoid,_Samplers])
