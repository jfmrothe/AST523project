c++    -c -o Point.o Point.cc
c++    -c -o Samplers.o Samplers.cc
c++    -c -o Ellipsoid.o Ellipsoid.cc
c++    -c -o Data.o Data.cc
c++    -c -o main.o main.cc
c++ -Wall -pg -o testmulti Point.o Samplers.o Ellipsoid.o Data.o main.o -lgsl -lgslcblas -lm
