#----------------------------------------------------------------
# Revised Makefile (19/02/2014 PS)
# How to call:
# calling as normal
#     > make
# will trigger the default compiler options.
#
# You can now also call using
#     > make CO=debug
# for extra warnings, gprof and gdb output, exception trapping
# at runtime, and bounds-checking.
# This option is slow (about 2x slower than make all)
#
# For working code, call
#     > make CO=fast
# to enable optimisation. This is about 2x faster than normal
# (using gfortran) but has little error trapping.
#
# Finally,
#     > make new
# simply calls clean then all to force a re-build.
#
# I have also included similar options for ifort. Since I have
# the compiler here, and it is potentially significantly faster,
# it may be useful when it comes time to do science.
#
#----------------------------------------------------------------

FC=gfortran
LD=gfortran
FFLAGS=-ffree-line-length-0 -Jsource/ -fopenmp
CUBEFLAGS=-L/software/cfitsio/ -lcfitsio -lm

ifeq ($(FC),gfortran)
  ifeq ($(CO),debug)
    FFLAGS += -fbounds-check -Wall -Wuninitialized #-ffpe-trap=zero,overflow,invalid,underflow,denormal
  endif
  ifeq ($(CO),debug2)
    FFLAGS += -g -pg -fbounds-check -Wall -Wuninitialized #-ffpe-trap=zero,overflow,invalid,underflow,denormal
  endif
  ifeq ($(CO),debug3)
    FFLAGS += -g -pg -fbounds-check -Wall -Wuninitialized -ffpe-trap=zero,overflow,invalid,underflow,denormal
  endif
  ifeq ($(CO),pedantic)
    FFLAGS += -g -pg -fbounds-check -Wall -Wuninitialized -Werror -pedantic -ffpe-trap=zero,overflow,invalid,underflow,denormal
  endif
  ifeq ($(CO),fast)
    FFLAGS += -O3 -fno-backtrace
  endif
endif

ifeq ($(FC),ifort)
  FFLAGS = -O0 #-warn all -warn errors
  LD=ifort
  ifeq ($(CO),debug)
    FFLAGS = -pg -g -check bounds -check uninit -warn all -warn nodeclarations -WB -zero -traceback # -std
  endif
  ifeq ($(CO),pedantic)
    FFLAGS = -pg -g -check bounds -check uninit -warn all -warn nodeclarations -WB -zero -traceback -std
  endif
  ifeq ($(CO),fast)
    FFLAGS = -axavx -msse3 -O3 -ip -ipo # for today's CPUs
#    FFLAGS = -fast -tune pn4 # for older pentium 4
  endif
endif

.PHONY: all clean install

new: clean all

all: alfa alfacube

%.o: %.f95
	$(FC) $(FFLAGS) $< -c -o $@

alfa: source/types.o source/functions.o source/readfiles.o source/quicksort.o source/continuum.o source/fit.o source/uncertainties.o source/alfa.o
	$(LD) $(LDFLAGS) $(FFLAGS) -o $@ $^

alfacube: source/types.o source/functions.o source/readfiles.o source/quicksort.o source/continuum.o source/fit.o source/uncertainties.o source/alfa_cube.o
	$(LD) $(LDFLAGS) $(FFLAGS) -o $@ $^ ${CUBEFLAGS}

clean:
	rm -f alfa alfacube source/*.o source/*.mod

install:

	test -e /etc/alfa || mkdir /etc/alfa
	install -m 644 linelists/* /etc/alfa
	install alfa /usr/bin
	install alfacube /usr/bin
