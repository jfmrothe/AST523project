OBJECT = Point.o Samplers.o Ellipsoid.o main.o 
PROJECT = testmulti
MULTI_INTERFACES = Point.i Samplers.i Ellipsoid.i
MULTI_SOURCES= Point.cc Point.h Samplers.cc Samplers.h Ellipsoid.cc Ellipsoid.h 
SOURCES = $(MULTI_SOURCES) 
INTERFACES = $(MULTI_INTERFACES) 
WRAPPERS = $(INTERFACES:.i=_wrap.cxx)
PROXIES = $(INTERFACES:.i=.py)
CXX = c++
.PHONY: all
all: $(WRAPPERS) $(SOURCES)
	./setup.py build_ext -i

%_wrap.cxx: %.i %.h ./numpy.i
	swig -c++ -python $<

.PHONY: cc
cc: $(PROJECT)

$(PROJECT) :$(OBJECT)
	$(CXX) -Wall -pg -o $@ $^



.PHONY: clean
clean:
	$(RM) $(PROJECT)
	$(RM) $(PROJECT_DEBUG)
	$(RM) *.so *.pyc *_wrap.h
	$(RM) -r *.dSYM
	$(RM) build/temp.linux-x86_64-2.7/*
	$(RM) $(WRAPPERS)
	$(RM) $(PROXIES)
	$(RM) *.o
	$(RM) .depend
	$(RM) *.out
	$(RM) *~
