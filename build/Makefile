#--------------------------------------------------------------
#
# Mescaline 1.0: a post-simulation analysis tool for cactus data
#
# (c) Daniel Price and Hayley Macpherson 2017-
#
#--------------------------------------------------------------

# -------------------------------------
# Directories containing Mescaline source code
# -------------------------------------
SRCDIR=../src
VPATH=$(SRCDIR) $(SPLASH_DIR)/src
SHELL=/bin/bash
BINDIR=../bin

# -------------------------------------
# Compiler and HDF5 library flags
# -------------------------------------

# GNU compilers
ifeq ($(SYSTEM), gfortran)
    FC=gfortran
    CC=gcc
    FFLAGS= -g -O3 -Wall -fdefault-real-8 -fopenmp -fcheck=bounds
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

# Location of some External Libraries needed for HDF5 Fortran interface
#LIBZ    = /Users/hayleymac/anaconda3/lib/libz.a
#LIBSZ   = /usr/local/lib/libsz.a
EXTLIBZ = -lm

# -------------------------------------
# Modules borrowed from splash
# -------------------------------------
SPLASH_FILES=asciiutils.f90 sort.f90 read_data_cactus_hdf5_futils.f90

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
MOD_FILES=params.f90 $(SPLASH_FILES) options.f90 random.f90 periodic.f90 manipulations.f90 prints.f90 derivatives.f90 gauge.f90 spines.f90 $(TRIP_FILE) init.f90 roots.f90 violation.f90

#
# Files used in the default mescaline case
MESC_FILES=$(MOD_FILES) ricci.f90 cactoid.f90 analysis.f90 mescaline.f90
OBJ_MESC=read_data_cactus_hdf5_utils.o $(MESC_FILES:.f90=.o)


# -------------------------------------
# Now compile...
# -------------------------------------

%.o: %.f90
	$(FC) $(FFLAGS) -o $@ -c $< $(HDF5INCLUDE)

mescaline: $(OBJ_MESC)
	$(FC) $(FFLAGS) -o $(BINDIR)/mescaline $(OBJ_MESC) $(LDFLAGS) $(HDF5LIBS) $(EXTLIBZ)

## Dependencies for splash HDF5 read
read_data_cactus_hdf5_utils.o: $(SPLASH_DIR)/src/read_data_cactus_hdf5_utils.c
	$(CC) $(CFLAGS) $(HDF5INCLUDE) $(H5PART_INCLUDE) -c $< -o $@

read_data_cactus_hdf5_utils.c: checkhdf5 checksplash

## Clean up some things
.PHONY: clean cleanhard
clean:
	rm *.o *.mod

cleanhard:
	rm -f *.o *.mod $(BINDIR)/mescaline $(BINDIR)/mescaline_ics $(BINDIR)/mescaline_test $(BINDIR)/mescaline_larena $(BINDIR)/mescaline_raytrace

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
