SHELL=/bin/sh

# This file contains a set of rules used by the "make"
#   command.  This makefile $(MAKEFILE) tells "make" how the
#   executable Point $(COMMAND) should be created
#   from the source files $(SRCS) via the object
#   files $(OBJS) and the header files $(HDRS); type the command:
#        "make -f make_Point"
#   where make_Point should be replaced by the name of
#   the makefile.
# 
# This file also tells make what to do if arguments
#   follow the "make" command.
#
# To remove the OBJS files; type the command:
#        "make -f make_Point clean"
#
# To create a gzipped tar file with name $(COMMAND).tar.gz
#  containing this makefile and the SRCS files,
#  type the command:
#        "make -f make_plane_parallel tarz"
 
SRCS= \
driver.cc \
  
OBJS= \
driver.o \

HDRS= \

MAKEFILE= Makefile

COMMAND= drive
 
CC= g++
CFLAGS= -g -O3
WARNFLAGS= -Werror -Wall -W -Wshadow -fno-common
MOREFLAGS= -ansi -pedantic -Wpointer-arith -Wcast-qual -Wcast-align \
           -Wwrite-strings -fshort-enums
LDFLAGS= -lgsl -lgslcblas # -lgomp
 
$(COMMAND): $(OBJS) 
	$(CC) -o $(COMMAND) $(OBJS) $(LDFLAGS) $(LIBS)

driver.o : driver.cc $(HDRS)
	$(CC) $(CFLAGS) $(WARNFLAGS)  -c driver.cc -o driver.o 
                 

clean:
	rm -f $(OBJS)
 
tarz:
	tar zcf - $(MAKEFILE) $(SRCS) $(HDRS) > $(COMMAND).tar.gz
