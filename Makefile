###########
#
# This file was originally part of the GADGET3 code developed by
#   Volker Springel. The code has been modified
#   slighty by Phil Hopkins (phopkins@caltech.edu) for GIZMO (mostly 
#   dealing with new files and filename conventions)
#
#############

CONFIG   =  Config.sh
PERL     =  /usr/bin/perl

RESULT     := $(shell CONFIG=$(CONFIG) PERL=$(PERL) make -f config-makefile)
CONFIGVARS := $(shell cat GIZMO_config.h)

HG_COMMIT := $(shell git rev-parse --short HEAD 2>/dev/null)
HG_REPO := $(shell git config --get remote.origin.url)
HG_BRANCH := $(shell git rev-parse --abbrev-ref HEAD 2>/dev/null)
HOST := $(shell hostname)
BUILDINFO = "Built on $(HOST) by $(USER) from $(HG_BRANCH):$(HG_COMMIT) at $(HG_REPO)"
OPT += -DBUILDINFO='$(BUILDINFO)'


# initialize some default flags -- these will all get re-written below
CC	= mpicc		# sets the C-compiler (default, will be set for machine below)
CXX	= mpiCC		# sets the C++-compiler (default, will be set for machine below)
FC	= mpif90	# sets the fortran compiler (default, will be set for machine below)
OPTIMIZE = -Wall  -g   # optimization and warning flags (default)


# one annoying thing here is the FFTW libraries, since they are named differently depending on
#  whether they are compiled in different precision levels, or with different parallelization options, so we
#  have to have a big block here 'sorting them out'.
#
ifeq (NOTYPEPREFIX_FFTW,$(findstring NOTYPEPREFIX_FFTW,$(CONFIGVARS)))  # fftw installed without type prefix?
    FFTW_LIBNAMES =  #-lrfftw_mpi -lfftw_mpi -lrfftw -lfftw
else
ifeq (DOUBLEPRECISION_FFTW,$(findstring DOUBLEPRECISION_FFTW,$(CONFIGVARS)))  # test for double precision libraries
    FFTW_LIBNAMES =  #-ldrfftw_mpi -ldfftw_mpi -ldrfftw -ldfftw
else
    FFTW_LIBNAMES =  #-lsrfftw_mpi -lsfftw_mpi -lsrfftw -lsfftw
endif
endif
# we only need fftw if PMGRID is turned on
ifneq (USE_FFTW3, $(findstring USE_FFTW3, $(CONFIGVARS)))
ifeq (PMGRID, $(findstring PMGRID, $(CONFIGVARS)))
ifeq (NOTYPEPREFIX_FFTW,$(findstring NOTYPEPREFIX_FFTW,$(CONFIGVARS)))  # fftw installed without type prefix?
  FFTW_LIBNAMES = -lrfftw_mpi -lfftw_mpi -lrfftw -lfftw
else
ifeq (DOUBLEPRECISION_FFTW,$(findstring DOUBLEPRECISION_FFTW,$(CONFIGVARS)))  # test for double precision libraries
  FFTW_LIBNAMES = -ldrfftw_mpi -ldfftw_mpi -ldrfftw -ldfftw
else
  FFTW_LIBNAMES = -lsrfftw_mpi -lsfftw_mpi -lsrfftw -lsfftw
endif
endif
else
# or if TURB_DRIVING_SPECTRUMGRID is activated
ifeq (TURB_DRIVING_SPECTRUMGRID, $(findstring TURB_DRIVING_SPECTRUMGRID, $(CONFIGVARS)))
ifeq (NOTYPEPREFIX_FFTW,$(findstring NOTYPEPREFIX_FFTW,$(CONFIGVARS)))  # fftw installed without type prefix?
  FFTW_LIBNAMES = -lrfftw_mpi -lfftw_mpi -lrfftw -lfftw
else
ifeq (DOUBLEPRECISION_FFTW,$(findstring DOUBLEPRECISION_FFTW,$(CONFIGVARS)))  # test for double precision libraries
  FFTW_LIBNAMES = -ldrfftw_mpi -ldfftw_mpi -ldrfftw -ldfftw
else
  FFTW_LIBNAMES = -lsrfftw_mpi -lsfftw_mpi -lsrfftw -lsfftw
endif
endif
else
  FFTW_LIBNAMES = #
endif
endif
else # use FFTW3 instead of FFTW2.?
ifeq (PMGRID, $(findstring PMGRID, $(CONFIGVARS)))
ifeq (DOUBLEPRECISION_FFTW,$(findstring DOUBLEPRECISION_FFTW,$(CONFIGVARS)))  # test for double precision libraries
  FFTW_LIBNAMES = -lfftw3_mpi -lfftw3
else #single precision 
  FFTW_LIBNAMES = -lfftw3f_mpi -lfftw3f
endif
else 
# or if TURB_DRIVING_SPECTRUMGRID is activated
ifeq (TURB_DRIVING_SPECTRUMGRID, $(findstring TURB_DRIVING_SPECTRUMGRID, $(CONFIGVARS)))
ifeq (DOUBLEPRECISION_FFTW,$(findstring DOUBLEPRECISION_FFTW,$(CONFIGVARS)))  # test for double precision libraries
  FFTW_LIBNAMES = -lfftw3_mpi -lfftw3
else #single precision  
  FFTW_LIBNAMES = -lfftw3f_mpi -lfftw3f
endif
else 
  FFTW_LIBNAMES = #
endif
endif
endif


## read the systype information to use the blocks below for different machines
ifdef SYSTYPE
SYSTYPE := "$(SYSTYPE)"
-include Makefile.systype
else
include Makefile.systype
endif

ifeq ($(wildcard Makefile.systype), Makefile.systype)
INCL = Makefile.systype
else
INCL =
endif



#----------------------------------------------------------------------------------------------
ifeq ($(SYSTYPE),"macOS")
export OMPI_CC = clang
export OMPI_CXX = clang++
export OMPI_FC = clang
CC       = mpicc
CXX      = mpicxx
FC       = mpif90
LDFLAGS	 = -L/usr/local/Cellar/gcc/11.1.0/lib/gcc/11 -lgfortran
# WARNING: do *NOT* run with -ffast-math !!
OPTIMIZE += -O2 -march=native -ffp-contract=off -fstandalone-debug #-fsanitize=address # clang options
ifeq (OPENMP,$(findstring OPENMP,$(CONFIGVARS)))
OPTIMIZE += -fopenmp # openmp required compiler flags
endif
GRACKLE_HOME = $(HOME)/grackle_install
GRACKLEINCL = -I$(GRACKLE_HOME)/include
GRACKLELIBS = -L$(GRACKLE_HOME)/lib -Wl,-rpath,$(GRACKLE_HOME)/lib
HDF5INCL = -I/usr/local/Cellar/hdf5/1.12.0_4/include -DH5_USE_16_API
HDF5LIB  = -L/usr/local/Cellar/hdf5/1.12.0_4/lib -lhdf5 -lz
GSLINCL  = -I/usr/local/Cellar/gsl/2.6/include
GSLLIB   = -L/usr/local/Cellar/gsl/2.6/lib -lgsl -lgslcblas
MPICHLIB = 
OPT     += -DUSE_MPI_IN_PLACE -DNO_ISEND_IRECV_IN_DOMAIN
endif
#----------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------
ifeq ($(SYSTYPE),"Ubuntu")
export OMPI_CC = clang
export OMPI_CXX = clang++
export OMPI_FC = clang
CC       = mpicc
CXX      = mpicc
FC       = mpif90
LDFLAGS	 = -lgfortran
# WARNING: do *NOT* run with -ffast-math !!
OPTIMIZE += -O2 -march=native -ffp-contract=off -fstandalone-debug #-fsanitize=address # clang options
ifeq (OPENMP,$(findstring OPENMP,$(CONFIGVARS)))
OPTIMIZE += -fopenmp # openmp required compiler flags
endif
GRACKLE_HOME = $(HOME)/grackle_install
GRACKLEINCL = -I$(GRACKLE_HOME)/include
GRACKLELIBS = -L$(GRACKLE_HOME)/lib -Wl,-rpath=$(GRACKLE_HOME)/lib
HDF5INCL = -I/usr/include/hdf5/serial -DH5_USE_16_API
HDF5LIB  = -L/usr/lib/x86_64-linux-gnu -lhdf5_serial -lz
GSLINCL  = -I/usr/include/gsl
GSLLIB   = -I/usr/lib/x86_64-linux-gnu -lgsl
MPICHLIB = 
OPT     += -DUSE_MPI_IN_PLACE -DNO_ISEND_IRECV_IN_DOMAIN
endif
#----------------------------------------------------------------------------------------------


# different code groups that need to be compiled. the groupings below are
# arbitrary (they will all be added to OBJS and compiled, and if they are
# un-used it should be through use of macro flags in the source code). But
# they are grouped below for the sake of clarity when adding/removing code
# blocks in the future
#
CORE_OBJS =	allvars.o grackle.o cooling.o eos.o main.o

## name of executable and optimizations
EXEC   = cooling_curve
OPTIONS = $(OPTIMIZE) $(OPT)

## combine all the objects above
OBJS  = $(CORE_OBJS)

## include files needed at compile time for the above objects
INCL    += 	allvars.h \
      proto.h \
			cooling.h \
			Makefile


## now we add special cases dependent on compiler flags. normally we would
##  include the files always, and simply use the in-file compiler variables
##  to determine whether certain code is compiled [this allows us to take
##  advantage of compiler logic, and makes it easier for the user to
##  always specify what they want]. However special cases can arise, if e.g.
##  there are certain special libraries needed, or external compilers, for
##  certain features

# if grackle libraries are installed they must be a shared library as defined here
ifeq (COOL_GRACKLE,$(findstring COOL_GRACKLE,$(CONFIGVARS)))
OPTIONS += -DCONFIG_BFLOAT_8
GRACKLEINCL +=
GRACKLELIBS += -lgrackle
else
GRACKLEINCL =
GRACKLELIBS =
endif

# linking libraries (includes machine-dependent options above)
CFLAGS = $(OPTIONS) $(GSLINCL) $(HDF5INCL) $(GRACKLEINCL)

LIBS = $(GSLLIB) $(HDF5LIB) -g -lm $(GRACKLELIBS)

ifeq (PTHREADS_NUM_THREADS,$(findstring PTHREADS_NUM_THREADS,$(CONFIGVARS))) 
LIBS += -lpthread
endif

$(EXEC): $(OBJS) $(FOBJS)  
	$(CXX) $(LDFLAGS) $(OPTIMIZE) $(OBJS) $(LIBS) $(RLIBS) -o $(EXEC)

$(OBJS): $(INCL) $(CONFIG)

clean:
	rm -f $(OBJS) $(FOBJS) $(EXEC) *.oo *.c~


