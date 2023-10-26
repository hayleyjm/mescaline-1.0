#--------------------------------------------------------------
# Mescaline: a post-simulation analysis tool for cactus data
# (c) Daniel Price and Hayley Macpherson 2017-
#
#     Makefile for averaging part only
#--------------------------------------------------------------

# -------------------------------------
# Directories containing Mescaline source code
# -------------------------------------
TYDIR=../src/tyrosine
AVGSDIR=../src/avgs
VPATH=$(AVGSDIR) $(TYDIR) $(SPLASH_DIR)/src
SHELL=/bin/bash

BINDIR=../bin

# -------------------------------------
# Compiler and HDF5 library flags
# -------------------------------------

# GNU compilers
ifeq ($(SYSTEM), gfortran)
    FC=gfortran
    CC=gcc
    FFLAGS= -g -O3 -Wall -fbacktrace -fdefault-real-8 -fdefault-double-8 -fopenmp #-fcheck=bounds -fcheck-all
    CFLAGS= -g
endif

# Intel compilers
ifeq ($(SYSTEM), ifort)
    FC=ifort
    CC=icc
    FFLAGS= -g -O3 -wrap-margin-
    CFLAGS= -g
endif

HDF5LIBS    = -L$(HDF5ROOT)/lib -lhdf5 -lhdf5_fortran
HDF5INCLUDE = -I$(HDF5ROOT)/include

# Some External Libraries needed for HDF5 Fortran interface
EXTLIBZ = -lm

# -------------------------------------
# Modules borrowed from splash: all should be located in $SPLASH_DIR/src
# -------------------------------------
SPLASH_FILES=asciiutils.f90 sort.f90 utils_vectors.f90 read_data_cactus_hdf5_futils.f90

# -------------------------------------
# Mescaline source files
# -------------------------------------
#

# If using ifort, need some special statement for output of averaging when using lots of spheres
ifeq ($(SYSTEM), ifort)
    TRIP_FILE=tripreports_ifort.f90
endif
ifeq ($(SYSTEM), gfortran)
    TRIP_FILE=tripreports.f90
endif

# Modules that are used in all make cases
MOD_FILES=params.f90 $(SPLASH_FILES) random.f90 periodic.f90 \
manipulations.f90 prints.f90 derivatives.f90 analytic_solns.f90 gauge.f90 \
spines.f90 $(TRIP_FILE) violation.f90 riemann.f90 weyl.f90

#
# Files used in the default mescaline case
AVGS_FILES=options_avgs.f90 $(MOD_FILES) init.f90 roots.f90 ricci.f90 cactoid.f90 analysis.f90 mescaline.f90
OBJ_AVGS=read_data_cactus_hdf5_utils.o $(AVGS_FILES:.f90=.o)

# -------------------------------------
# Now compile...
# -------------------------------------

%.o: %.f90
	$(FC) $(FFLAGS) -o $@ -c $< $(HDF5INCLUDE)

## Default - ``regular'' Mescaline case, using the new Buchert+2019 formalism
.PHONY: mescaline
mescaline: $(OBJ_AVGS)
	$(FC) $(FFLAGS) -o $(BINDIR)/mescaline $(OBJ_AVGS) $(LDFLAGS) $(HDF5LIBS) $(EXTLIBZ)


## Dependencies for splash HDF5 read
read_data_cactus_hdf5_utils.o: $(SPLASH_DIR)/src/read_data_cactus_hdf5_utils.c
	$(CC) $(CFLAGS) $(HDF5INCLUDE) $(H5PART_INCLUDE) -c $< -o $@

read_data_cactus_hdf5_utils.c: checkhdf5 checksplash

## Clean up some things 
.PHONY: clean cleanhard
clean:
	rm *.o *.mod

cleanhard:
	rm *.o *.mod $(BINDIR)/mescaline*

checkhdf5:
   ifeq (X${HDF5ROOT}, X)
	@echo; echo "ERROR: Mescaline needs HDF5ROOT to be set before compiling"; echo; echo "try: export HDF5ROOT=/usr/local"; echo; ${MAKE} err;
   else
	@if [ -d $$HDF5ROOT ]; then echo; echo "HDF5ROOT=$$HDF5ROOT"; echo; else echo; echo "ERROR: Directory given by HDF5ROOT=$$HDF5ROOT does not exist"; echo; ${MAKE} err; fi;
   endif

checksplash:
   ifeq (X${SPLASH_DIR}, X)
	@echo; echo "ERROR: SPLASH_DIR should be set before compiling"; echo; ${MAKE} err;
   else
	@if [ -d $$SPLASH_DIR ]; then echo; echo "SPLASH_DIR=$$SPLASH_DIR"; echo; else echo; echo "ERROR: Directory given by SPLASH_DIR=$$SPLASH_DIR does not exist"; echo; ${MAKE} err; fi;
   endif

err:
	$(error aborting);
