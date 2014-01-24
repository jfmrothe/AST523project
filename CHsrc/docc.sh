c++    -g -c Point.cc -o Point.o -I/u/chelsea/Install/include/ -L/u/chelsea/Install/lib
c++    -g -c Samplers.cc -o Samplers.o -I/u/chelsea/Install/include/ -L/u/chelsea/Install/lib
c++    -g -c Ellipsoid.cc -o Ellipsoid.o  -I/u/chelsea/Install/include/ -L/u/chelsea/Install/lib
c++    -g -c main.cc -o main.o -I/u/chelsea/Install/include/ -L/u/chelsea/Install/lib
c++ -Wall -g -pg -o testmulti Point.o Samplers.o Ellipsoid.o main.o -lgsl -lgslcblas -lm


