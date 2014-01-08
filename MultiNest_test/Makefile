#!/bin/shell

OBJS= \
main.o \
Point.o \
Data.o \
Samplers.o \
SamplersUtilities.o \
Ellipsoid.o

HDRS= \
Point.h \
Data.h \
Samplers.h \
Ellipsoid.h

CC=g++
CFLAGS= -g -O3
WARNFLAGS= -Werror -Wall -W -Wshadow -fno-common
LDFLAGS= -lgsl -lgslcblas -lm

all: multinest 

multinest : $(OBJS)
	$(CC) $(CFLAGS) $(LDFLAGS) $(OBJS) -o multinest

main.o : main.cc $(HDRS) 
	$(CC) $(CFLAGS) -c main.cc

Point.o : Point.cc $(HRDS)
	$(CC) $(CFLAGS) -c Point.cc

Data.o : Data.cc $(HDRS)
	$(CC) $(CFLAGS) -c Data.cc

Samplers.o : Samplers.cc $(HDRS)
	$(CC) $(CFLAGS) -c Samplers.cc

SamplersUtilities.o : SamplersUtilities.cc $(HDRS)
	$(CC) $(CFLAGS) -c SamplersUtilities.cc

Elipsoid.o : Ellipsoid.cc $(HDRS)
	$(CC) $(CFLAGS) -c Ellipsoid.cc

clean:
	rm *.o multinest 
