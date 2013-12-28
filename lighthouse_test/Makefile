#!/bin/shell
CC=g++
CFLAGS= -g -O3
WARNFLAGS= -Werror -Wall -W -Wshadow -fno-common
LDFLAGS= -lgsl -lgslcblas -lm

all: drive

drive : driver.o EllSamp.o 
	$(CC) $(CFLAGS) $(LDFLAGS) driver.o EllSamp.o -o drive

drive.o : driver.cc Point.h
	$(CC) $(CFLAGS) -c driver.cc

EllSamp.o : EllSamp.cc
	$(CC) $(CFLAGS) -c EllSamp.cc
 
#findenclosing.o : findenclosing.cc
#	$(CC) $(CFLAGS) -c findenclosing.cc

clean:
	rm *.o drive
