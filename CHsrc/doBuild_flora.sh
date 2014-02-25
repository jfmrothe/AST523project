swig -c++ -python Point.i
swig -c++ -python Samplers.i
swig -c++ -python Ellipsoid.i
mkdir build
mkdir build/temp.linux-x86_64-2.7
gcc -fno-strict-aliasing -g -O2 -DNDEBUG -g -fwrapv -O3 -Wall -Wstrict-prototypes -fPIC -I/usr/peyton/utils/Python-2.7.2_rwx/lib/python2.7/site-packages/numpy/core/include -I/usr/peyton/utils/Python-2.7.2_rwx/include/python2.7 -I/u/chelsea/Install/include/ -c Point_wrap.cxx -o build/temp.linux-x86_64-2.7/Point_wrap.o -L/u/chelsea/Install/lib 
gcc -fno-strict-aliasing -g -O2  -DNDEBUG -g -fwrapv -O3 -Wall -Wstrict-prototypes -fPIC -I/usr/peyton/utils/Python-2.7.2_rwx/lib/python2.7/site-packages/numpy/core/include -I/usr/peyton/utils/Python-2.7.2_rwx/include/python2.7 -I/u/chelsea/Install/include/ -c Point.cc -o build/temp.linux-x86_64-2.7/Point.o -L/u/chelsea/Install/lib/ 
g++ -pthread -shared build/temp.linux-x86_64-2.7/Point_wrap.o build/temp.linux-x86_64-2.7/Point.o -o _Point.so -lgsl -lgslcblas -lm
#building '_Ellipsoid' extension
gcc -fno-strict-aliasing -g -O2  -DNDEBUG -g -fwrapv -O3 -Wall -Wstrict-prototypes -fPIC -I/usr/peyton/utils/Python-2.7.2_rwx/lib/python2.7/site-packages/numpy/core/include -I/usr/peyton/utils/Python-2.7.2_rwx/include/python2.7 -c Ellipsoid_wrap.cxx -o build/temp.linux-x86_64-2.7/Ellipsoid_wrap.o -L/usr/lib/ -L/u/chelsea/Install/lib/ -I/usr/include/
gcc -fno-strict-aliasing -g -O2  -DNDEBUG -g -fwrapv -O3 -Wall -Wstrict-prototypes -fPIC -I/usr/peyton/utils/Python-2.7.2_rwx/lib/python2.7/site-packages/numpy/core/include -I/usr/peyton/utils/Python-2.7.2_rwx/include/python2.7 -c Ellipsoid.cc -o build/temp.linux-x86_64-2.7/Ellipsoid.o -L/u/chelsea/Install/lib/ -I/usr/include/
g++ -pthread -shared build/temp.linux-x86_64-2.7/Ellipsoid_wrap.o build/temp.linux-x86_64-2.7/Ellipsoid.o build/temp.linux-x86_64-2.7/Point.o -o _Ellipsoid.so -lgsl -lgslcblas -lm
#building '_Samplers' extension
gcc -fno-strict-aliasing -g -O2 -DNDEBUG -g -fwrapv -O3 -Wall -Wstrict-prototypes -fPIC -I/usr/peyton/utils/Python-2.7.2_rwx/lib/python2.7/site-packages/numpy/core/include -I/usr/peyton/utils/Python-2.7.2_rwx/include/python2.7 -c Samplers_wrap.cxx -o build/temp.linux-x86_64-2.7/Samplers_wrap.o -L/u/chelsea/Install/lib/ -I/u/chelsea/Install/include/
gcc -fno-strict-aliasing -g -O2 -DNDEBUG -g -fwrapv -O3 -Wall -Wstrict-prototypes -fPIC -I/usr/peyton/utils/Python-2.7.2_rwx/lib/python2.7/site-packages/numpy/core/include -I/usr/peyton/utils/Python-2.7.2_rwx/include/python2.7 -c Samplers.cc -o build/temp.linux-x86_64-2.7/Samplers.o -L/u/chelsea/Install/lib/ -I/u/chelsea/Install/include/
g++ -pthread -shared build/temp.linux-x86_64-2.7/Samplers_wrap.o build/temp.linux-x86_64-2.7/Ellipsoid.o build/temp.linux-x86_64-2.7/Point.o build/temp.linux-x86_64-2.7/Samplers.o -o _Samplers.so -lgsl -lgslcblas -lm
