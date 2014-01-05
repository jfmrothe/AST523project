#!/bin/shell

OBJS= \
main.o \
Point.o \
Data.o \
Samplers.o

HDRS= \
Point.h \
Data.h \
Samplers.h

CC=g++
CFLAGS= -g -O3
WARNFLAGS= -Werror -Wall -W -Wshadow -fno-common
LDFLAGS= -lgsl -lgslcblas -lm

all: run

run : $(OBJS)
	$(CC) $(CFLAGS) $(LDFLAGS) $(OBJS) -o run

main.o : main.cc $(HDRS) 
	$(CC) $(CFLAGS) -c main.cc

Point.o : Point.cc $(HRDS)
	$(CC) $(CFLAGS) -c Point.cc

Data.o : Data.cc $(HDRS)
	$(CC) $(CFLAGS) -c Data.cc

Samplers.o : Samplers.cc $(HDRS)
	$(CC) $(CFLAGS) -c Samplers.cc

clean:
	rm *.o run
