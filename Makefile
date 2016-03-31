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
PREFIX=/usr
MANDIR=${DESTDIR}${PREFIX}/share/man/man1

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

all: alfa alfacube alfarss

new: clean all

%.o: %.f95
	$(FC) $(FFLAGS) $< -c -o $@

alfa: source/types.o source/functions.o source/readfiles.o source/quicksort.o source/continuum.o source/fit.o source/uncertainties.o source/alfa.o
	$(LD) $(LDFLAGS) $(FFLAGS) -o $@ $^

alfacube: source/types.o source/functions.o source/readfiles.o source/quicksort.o source/continuum.o source/fit.o source/uncertainties.o source/alfa_cube.o
	$(LD) $(LDFLAGS) $(FFLAGS) -o $@ $^ ${CUBEFLAGS}

alfarss: source/types.o source/functions.o source/readfiles.o source/quicksort.o source/continuum.o source/fit.o source/uncertainties.o source/alfa_rss.o
	$(LD) $(LDFLAGS) $(FFLAGS) -o $@ $^ ${CUBEFLAGS}

clean:
	rm -f alfa alfacube source/*.o source/*.mod

install:
	test -e ${DESTDIR}${PREFIX}/share/alfa || mkdir -p ${DESTDIR}${PREFIX}/share/alfa
	test -e ${DESTDIR}${PREFIX}/bin || mkdir -p ${DESTDIR}${PREFIX}/bin
	test -e ${MANDIR} || mkdir -p ${MANDIR}
	install -m 644 linelists/* ${DESTDIR}${PREFIX}/share/alfa
	install alfa ${DESTDIR}${PREFIX}/bin
	install alfacube ${DESTDIR}${PREFIX}/bin
	install alfarss ${DESTDIR}${PREFIX}/bin
	install -g 0 -o 0 -m 644 man/alfa.1 ${MANDIR}
	gzip -f ${MANDIR}/alfa.1
	ln -s -f ${MANDIR}/alfa.1.gz ${MANDIR}/alfacube.1.gz

uninstall:
	rm -rf ${DESTDIR}${PREFIX}/share/alfa
	rm -f ${DESTDIR}${PREFIX}/bin/alfa ${DESTDIR}${PREFIX}/bin/alfacube ${DESTDIR}${PREFIX}/bin/alfarss
	rm -f ${MANDIR}/alfa.1.gz ${MANDIR}/alfacube.1.gz
