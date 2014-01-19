swig -c++ -python Point.i
swig -c++ -python Ellipsoid.i
swig -c++ -python Samplers.i
./setup.py build_ext -i
running build_ext
building '_Point' extension
creating build
creating build/temp.macosx-10.6-universal-2.6
gcc-4.2 -fno-strict-aliasing -fno-common -dynamic -DNDEBUG -g -fwrapv -Os -Wall -Wstrict-prototypes -DENABLE_DTRACE -arch x86_64 -pipe -I/Library/Python/2.6/site-packages/numpy/core/include -I/System/Library/Frameworks/Python.framework/Versions/2.6/include/python2.6 -I/usr/local/include/ -c Point_wrap.cxx -o build/temp.macosx-10.6-universal-2.6/Point_wrap.o -L/usr/local/lib/ 
gcc-4.2 -fno-strict-aliasing -fno-common -dynamic -DNDEBUG -g -fwrapv -Os -Wall -Wstrict-prototypes -DENABLE_DTRACE -arch x86_64 -pipe -I/Library/Python/2.6/site-packages/numpy/core/include -I/System/Library/Frameworks/Python.framework/Versions/2.6/include/python2.6 -I/usr/local/include/ -c Point.cc -o build/temp.macosx-10.6-universal-2.6/Point.o -L/usr/local/lib/
g++-4.2 -Wl,-F. -bundle -undefined dynamic_lookup -arch x86_64 build/temp.macosx-10.6-universal-2.6/Point_wrap.o build/temp.macosx-10.6-universal-2.6/Point.o -o _Point.so -lgsl -lgslcblas -lm
gcc-4.2 -fno-strict-aliasing -fno-common -dynamic -DNDEBUG -g -fwrapv -Os -Wall -Wstrict-prototypes -DENABLE_DTRACE -arch x86_64 -pipe -I/Library/Python/2.6/site-packages/numpy/core/include -I/System/Library/Frameworks/Python.framework/Versions/2.6/include/python2.6 -I/usr/local/include -c Ellipsoid_wrap.cxx -o build/temp.macosx-10.6-universal-2.6/Ellipsoid_wrap.o -L/usr/local/lib
gcc-4.2 -fno-strict-aliasing -fno-common -dynamic -DNDEBUG -g -fwrapv -Os -Wall -Wstrict-prototypes -DENABLE_DTRACE -arch x86_64 -pipe -I/Library/Python/2.6/site-packages/numpy/core/include -I/System/Library/Frameworks/Python.framework/Versions/2.6/include/python2.6 -I/usr/local/include -c Ellipsoid.cc -o build/temp.macosx-10.6-universal-2.6/Ellipsoid.o -L/usr/local/lib
g++-4.2 -Wl,-F. -bundle -undefined dynamic_lookup -arch x86_64 build/temp.macosx-10.6-universal-2.6/Ellipsoid_wrap.o build/temp.macosx-10.6-universal-2.6/Ellipsoid.o -o _Ellipsoid.so -lgsl -lgslcblas -lm
gcc-4.2 -fno-strict-aliasing -fno-common -dynamic -DNDEBUG -g -fwrapv -Os -Wall -Wstrict-prototypes -DENABLE_DTRACE -arch x86_64 -pipe -I/Library/Python/2.6/site-packages/numpy/core/include -I/System/Library/Frameworks/Python.framework/Versions/2.6/include/python2.6 -I/usr/local/include -c Samplers_wrap.cxx -o build/temp.macosx-10.6-universal-2.6/Samplers_wrap.o -L/usr/local/lib
gcc-4.2 -fno-strict-aliasing -fno-common -dynamic -DNDEBUG -g -fwrapv -Os -Wall -Wstrict-prototypes -DENABLE_DTRACE -arch x86_64 -pipe -I/Library/Python/2.6/site-packages/numpy/core/include -I/System/Library/Frameworks/Python.framework/Versions/2.6/include/python2.6 -I/usr/local/include -c Samplers.cc -o build/temp.macosx-10.6-universal-2.6/Samplers.o -L/usr/local/lib
g++-4.2 -Wl,-F. -bundle -undefined dynamic_lookup -arch x86_64 build/temp.macosx-10.6-universal-2.6/Samplers_wrap.o build/temp.macosx-10.6-universal-2.6/Samplers.o -o _Samplers.so -lgsl -lgslcblas -lm
