OBJECT = oblatenessfast.o main.o
PROJECT = oblateness
MULTI_INTERFACES = Points.i Samplers.i Ellipsoid.i
MULTI_SOURCES= Points.cc Points.h Samplers.cc Samplers.h Ellipsoid.cc Ellipsoid.h 
SOURCES = $(MULTI_SOURCES) 
INTERFACES = $(MULTI_INTERFACES) 
WRAPPERS = $(INTERFACES:.i=_wrap.cxx)
PROXIES = $(INTERFACES:.i=.py)

.PHONY: all
all: $(WRAPPERS) $(SOURCES)
	./setup.py build_ext -i

%_wrap.cxx: %.i %.h ./numpy.i
	swig -c++ -python $<


.PHONY: cc
cc: $(PROJECT)

$(PROJECT) :$(OBJECT)
	$(CXX) -pg -o $@ $^


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
	$(RM) .depend
	$(RM) *.out
	$(RM) *~
