multinest:multinest.o findenclosing.o sampleellipsoid.o cluster.o
	g++ -g -o multinest multinest.o findenclosing.o sampleellipsoid.o cluster.o -lgsl -lgslcblas -lm
multinest.o:multinest.cc
	g++ -I /home/johannes/Documents/USA/courses/algo/finalproject/src_develop -g -c multinest.cc
sampleellipsoid:sampleellipsoid.o
	g++ -g -o sampleellipsoid sampleellipsoid.o -lgsl -lgslcblas -lm
sampleellipsoid.o:sampleellipsoid.cc
	g++ -I /home/johannes/Documents/USA/courses/algo/finalproject/src_develop -g -c sampleellipsoid.cc
findenclosing:findenclosing.o
	g++ -g -o findenclosing findenclosing.o -lgsl -lgslcblas -lm
findenclosing.o:findenclosing.cc
	g++ -I /home/johannes/Documents/USA/courses/algo/finalproject/src_develop -g -c findenclosing.cc
cluster:cluster.o
	g++ -g -o cluster cluster.o -lgsl -lgslcblas -lm
cluster.o:cluster.cc
	g++ -I /home/johannes/Documents/USA/courses/algo/finalproject/src_develop -g -c cluster.cc

clean:
	rm multinest.o multinest
	rm sampleellipsoid.o sampleellipsoid
	rm findenclosing.o findenclosing