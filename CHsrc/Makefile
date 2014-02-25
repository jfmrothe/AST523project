OBJECT = Point.o Samplers.o Ellipsoid.o Data.o main.o 
PROJECT = testmulti
MULTI_INTERFACES = Point.i Samplers.i Ellipsoid.i
MULTI_SOURCES= Point.cc Point.h Samplers.cc Samplers.h Ellipsoid.cc Ellipsoid.h 
SOURCES = $(MULTI_SOURCES) 
INTERFACES = $(MULTI_INTERFACES) 
WRAPPERS = $(INTERFACES:.i=_wrap.cxx)
PROXIES = $(INTERFACES:.i=.py)
CXX = c++
LDFLAGS = -pg 
CFLAGS = -g -Wall   
LIBS = -lgsl -lgslcblas -lm
HEAD = Point.h Samplers.h Ellipsoid.h
INCLUDES = 
CCSOURCES = $(MULTI_SOURCES) main.cc Data.cc Data.h
.PHONY: all
all: $(WRAPPERS) $(SOURCES)
	./setup.py build_ext -i

%_wrap.cxx: %.i %.h ./numpy.i
	swig -c++ -python $<

.PHONY: cc
cc: $(PROJECT)

$(PROJECT) :$(OBJECT)
	$(CXX) $(CFLAGS) -o $@ $^ $(LIBS)

$(OBJECT): $(CCSOURCES)
	$(CXX) $(INCLUDES) $(LDFLAGS) $(CFLAGS) -c $^

.PHONY: clean
clean:
	$(RM) $(PROJECT)
	$(RM) $(PROJECT_DEBUG)
	$(RM) *.so *.pyc *_wrap.h
	$(RM) -r *.dSYM
	$(RM) -r build
	$(RM) $(WRAPPERS)
	$(RM) $(PROXIES)
	$(RM) *.o
	$(RM) *.gch
	$(RM) .depend
	$(RM) *.out
	$(RM) *~
