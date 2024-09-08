# Generated automatically from Makefile.in by configure.
#
#   This file is part of CONCORDE
#
#   (c) Copyright 1995--1999 by David Applegate, Robert Bixby,
#   Vasek Chvatal, and William Cook
#
#   Permission is granted for academic research use.  For other uses,
#   contact the authors for licensing options.
#
#   Use at your own risk.  We make no guarantees about the
#   correctness or usefulness of this code.
#


SHELL = /bin/sh
SRCROOT = ..
BLDROOT = ..
CCINCDIR=$(SRCROOT)/INCLUDE

srcdir = .

CC = gcc
CFLAGS = -O3  -I$(BLDROOT)/INCLUDE -I$(CCINCDIR)
LDFLAGS = -O3 
LIBFLAGS = -lm 
RANLIB = ranlib

OBJ_SUFFIX = o
o = $(OBJ_SUFFIX)

THISLIB=linkern.a
LIBSRCS=linkern.c linkern_fixed.c linkern_path.c flip_two.c
ALLSRCS=lk_main.c $(LIBSRCS)

LIBS=$(BLDROOT)/EDGEGEN/edgegen.a $(BLDROOT)/KDTREE/kdtree.a \
     $(BLDROOT)/FMATCH/fmatch.a   $(BLDROOT)/UTIL/util.a

all: linkern LKSR $(THISLIB)

everything: all

linkern: lk_main.$o $(THISLIB) $(LIBS)
	$(CC) $(LDFLAGS) -o $@ lk_main.$o $(THISLIB) $(LIBS) $(LIBFLAGS)

LKSR: 
	gcc ttp_util.c ttp_heuristic.c ttp_final_ga.c ../UTIL/allocrus.c logistic_regression.c math_util.c ttp_optimisation.c ../KDTREE/kdbuild.c ../KDTREE/kdnear.c ../KDTREE/kdspan.c ../UTIL/zeit.c ../UTIL/getdata.c flip_two.c ../UTIL/edgelen.c ../UTIL/dheaps_i.c ../UTIL/urandom.c ../UTIL/safe_io.c ../UTIL/fastread.c ../UTIL/sortrus.c ../UTIL/edgeutil.c ../UTIL/edgemap.c ../UTIL/util.c ttp.c ttp_linkern.c linkern.c mbfs.c cvector.c ../EDGEGEN/delaunay.c -I../INCLUDE -o LKSR -lm -O3

clean:
	-rm -f *.$o $(THISLIB) linkern

OBJS=$(LIBSRCS:.c=.o)

$(THISLIB): $(OBJS)
	$(AR) $(ARFLAGS) $(THISLIB) $(OBJS)
	$(RANLIB) $(THISLIB)

.PHONY: $(BLDROOT)/concorde.a
$(BLDROOT)/concorde.a: $(OBJS)
	$(AR) $(ARFLAGS) $(BLDROOT)/concorde.a $(OBJS)
	$(RANLIB) $(BLDROOT)/concorde.a

include ../INCLUDE/Makefile.common

# DO NOT DELETE THIS LINE -- make depend depends on it.

I=$(CCINCDIR)
I2=$(BLDROOT)/INCLUDE

flip_two.$o: flip_two.c $(I)/machdefs.h $(I2)/config.h  $(I)/util.h     \
        $(I)/linkern.h  
linkern.$o:  linkern.c  $(I)/machdefs.h $(I2)/config.h  $(I)/linkern.h  \
        $(I)/util.h     $(I)/kdtree.h   $(I)/macrorus.h 
lk_main.$o:  lk_main.c  $(I)/machdefs.h $(I2)/config.h  $(I)/linkern.h  \
        $(I)/util.h     $(I)/kdtree.h   $(I)/edgegen.h  $(I)/macrorus.h 
ttp_final_ga.$o: ttp_final_ga.c $(I)/

#include "ttp_util.h"
#include "kdtree.h"
#include "ttp_heuristic.h"
#include "ttp_optimisation.h"
#include "linkern.h"
#include "ttp.h"
#include "ttp_linkern.h"
#include "mbfs.h"

# gcc ttp_util.c ttp_heuristic.c test.c ../UTIL/allocrus.c logistic_regression.c math_util.c ttp_optimisation.c ../KDTREE/kdbuild.c ../KDTREE/kdnear.c ../KDTREE/kdspan.c ../UTIL/zeit.c ../UTIL/getdata.c flip_two.c ../UTIL/edgelen.c ../UTIL/dheaps_i.c ../UTIL/urandom.c ../UTIL/safe_io.c ../UTIL/fastread.c ../UTIL/sortrus.c ../UTIL/edgeutil.c ../UTIL/edgemap.c ../UTIL/util.c ttp.c ttp_linkern.c linkern.c mbfs.c cvector.c ../EDGEGEN/delaunay.c -I../INCLUDE -o testpgch2 -lm -O3