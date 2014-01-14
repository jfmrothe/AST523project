#! /usr/bin/env python

from distutils.core import *
from distutils import sysconfig
from distutils.extension import Extension

import numpy

numpy_include = numpy.get_include()
_Points = Extension("_Points",
        ["Points_wrap.cxx",
            "Points.cc"],
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
        description = "Baysen tools for model fitting",
        author = "Johnny Greco, Johannes Rothe for the c++ part, Xu Huang for the python part and modeling",
        author_email = "",
        url = "",
        version = "0.0.0",
        py_modules = ["oblateness","elliptic","oblatenessfast"],
        ext_modules = [_oblateness,_elliptic,_oblatenessfast])
