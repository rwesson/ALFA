#----------------------------------------------------------------
# Makefile based on the one written for NEAT by Peter Scicluna
# How to call:
# calling as normal
#     > make
# will trigger the default compiler options which include optimisation.
#
# You can now also call using
#     > make CO=debug
# for extra warnings, gprof and gdb output, exception trapping
# at runtime, and bounds-checking.
# This option is slow (about 2x slower than make all)
#
#     > make CO=debug2, debug2, pedantic
# offer further levels of checks in case of problems
#
#     > make new
# simply calls clean then all to force a re-build.
#
#     > (sudo) make install
# places the files in the standard UNIX directories.
#----------------------------------------------------------------

FC=gfortran
LD=gfortran
FFLAGS=-ffree-line-length-0 -Jsource/ -fopenmp
CUBEFLAGS=-L/usr/lib/x86_64-linux-gnu/ -lcfitsio -lm
MANDIR=${DESTDIR}/usr/share/man/man1

ifeq ($(FC),gfortran)
  ifeq ($(CO),debug)
    FFLAGS += -fbounds-check -Wall -Wuninitialized #-ffpe-trap=zero,overflow,invalid,underflow,denormal
  else ifeq ($(CO),debug2)
    FFLAGS += -g -pg -fbounds-check -Wall -Wuninitialized #-ffpe-trap=zero,overflow,invalid,underflow,denormal
  else ifeq ($(CO),debug3)
    FFLAGS += -g -pg -fbounds-check -Wall -Wuninitialized -ffpe-trap=zero,overflow,invalid,underflow,denormal
  else ifeq ($(CO),pedantic)
    FFLAGS += -g -pg -fbounds-check -Wall -Wuninitialized -Werror -pedantic -ffpe-trap=zero,overflow,invalid,underflow,denormal
  else
    FFLAGS += -O3 -fno-backtrace
  endif
endif

ifeq ($(FC),ifort)
  FFLAGS = -O0 #-warn all -warn errors
  LD=ifort
  ifeq ($(CO),debug)
    FFLAGS = -pg -g -check bounds -check uninit -warn all -warn nodeclarations -WB -zero -traceback # -std
  else ifeq ($(CO),pedantic)
    FFLAGS = -pg -g -check bounds -check uninit -warn all -warn nodeclarations -WB -zero -traceback -std
  else
    FFLAGS = -axavx -msse3 -O3 -ip -ipo # for today's CPUs
#    FFLAGS = -fast -tune pn4 # for older pentium 4
  endif
endif

.PHONY: all clean install

all: alfa alfacube

new: clean all

%.o: %.f95
	$(FC) $(FFLAGS) $< -c -o $@

alfa: source/types.o source/functions.o source/readfiles.o source/quicksort.o source/continuum.o source/fit.o source/uncertainties.o source/alfa.o
	$(LD) $(LDFLAGS) $(FFLAGS) -o $@ $^

alfacube: source/types.o source/functions.o source/readfiles.o source/quicksort.o source/continuum.o source/fit.o source/uncertainties.o source/alfa_cube.o
	$(LD) $(LDFLAGS) $(FFLAGS) -o $@ $^ ${CUBEFLAGS}

clean:
	rm -f alfa alfacube source/*.o source/*.mod

install:
	test -e ${DESTDIR}/usr/share/alfa || mkdir -p ${DESTDIR}/usr/share/alfa
	test -e ${DESTDIR}/usr/bin || mkdir -p ${DESTDIR}/usr/bin
	test -e ${MANDIR} || mkdir -p ${MANDIR}
	install -m 644 linelists/* ${DESTDIR}/usr/share/alfa
	install alfa ${DESTDIR}/usr/bin
	install alfacube ${DESTDIR}/usr/bin
	install -g 0 -o 0 -m 644 man/alfa.1 ${MANDIR}
	gzip -f ${MANDIR}/alfa.1
	ln -s -f ${MANDIR}/alfa.1.gz ${MANDIR}/alfacube.1.gz

uninstall:
	rm -rf ${DESTDIR}/usr/share/alfa
	rm -f ${DESTDIR}/usr/bin/alfa ${DESTDIR}/usr/bin/alfacube
	rm -f ${MANDIR}/alfa.1.gz ${MANDIR}/alfacube.1.gz
