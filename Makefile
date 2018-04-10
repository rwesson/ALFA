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
#     > make CO=debug2, debug3, pedantic
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

# set prefix depending on OS
OS := $(shell uname)
ifeq ($(OS),Darwin)
  PREFIX=/usr/local
else
  PREFIX=/usr
endif

# get version from changelog if debian package, or git log otherwise
VERSION := $(shell if [ -e debian/ ]; then dpkg-parsechangelog -S version; elif [ -e .git/ ]; then git describe --always --tags --dirty; fi)

FFLAGS+=-cpp -DPREFIX=\"$(PREFIX)\" -DVERSION=\"$(VERSION)\"
LDFLAGS+=
CFITSIOFLAGS=-lcfitsio -lm #-L/usr/lib/x86_64-linux-gnu/
MANDIR=$(DESTDIR)$(PREFIX)/share/man/man1

ifeq ($(FC),gfortran)
  FFLAGS += -ffree-line-length-0 -Jsource/ -fopenmp
  ifeq ($(CO),debug)
    FFLAGS += -fbounds-check -Wall -Wuninitialized -DCO=\"$(CO)\" #-ffpe-trap=zero,overflow,invalid,underflow,denormal
  else ifeq ($(CO),debug2)
    FFLAGS += -g -pg -fbounds-check -Wall -Wuninitialized -DCO=\"$(CO)\" #-ffpe-trap=zero,overflow,invalid,underflow,denormal
  else ifeq ($(CO),debug3)
    FFLAGS += -g -pg -fbounds-check -Wall -Wuninitialized -ffpe-trap=zero,overflow,invalid,underflow,denormal -DCO=\"$(CO)\"
  else ifeq ($(CO),pedantic)
    FFLAGS += -g -pg -fbounds-check -Wall -Wuninitialized -Werror -pedantic -ffpe-trap=zero,overflow,invalid,underflow,denormal -DCO=\"$(CO)\"
  else
    FFLAGS += -O3 -fno-backtrace
  endif
endif

ifeq ($(FC),ifort)
  FFLAGS += -module source/ -openmp
  LD=ifort
  ifeq ($(CO),debug)
    FFLAGS += -pg -g -check bounds -check uninit -warn all -warn nodeclarations -WB -zero -traceback -DCO=\"$(CO)\" # -std
  else ifeq ($(CO),pedantic)
    FFLAGS += -pg -g -check bounds -check uninit -warn all -warn nodeclarations -WB -zero -traceback -std -DCO=\"$(CO)\"
  else
    FFLAGS += -O3 -ip -ipo # for today's CPUs
#    FFLAGS = -fast -tune pn4 # for older pentium 4
  endif
endif

.PHONY: all clean install new

all: alfa

new: clean all

%.o: %.f90
	$(FC) $(FFLAGS) $< -c -o $@

alfa: source/rnglib.o source/types.o source/functions.o source/commandline.o source/readfiles.o source/quicksort.o source/continuum.o source/linefit.o source/uncertainties.o source/alfa.o
	$(LD) $(LDFLAGS) $(FFLAGS) -o $@ $^ $(CFITSIOFLAGS)
	@echo "Compilation complete"

clean:
	rm -f alfa source/*.o source/*.mod man/alfa.html

install: alfa
	test -e $(DESTDIR)$(PREFIX)/share/alfa || mkdir -p $(DESTDIR)$(PREFIX)/share/alfa
	test -e $(DESTDIR)$(PREFIX)/bin || mkdir -p $(DESTDIR)$(PREFIX)/bin
	test -e $(MANDIR) || mkdir -p $(MANDIR)
	install -m 644 linelists/* $(DESTDIR)$(PREFIX)/share/alfa
	install alfa $(DESTDIR)$(PREFIX)/bin
	install -m 644 man/alfa.1 $(MANDIR)
	test -e $(DESTDIR)$(PREFIX)/share/bash-completion/completions || mkdir -p $(DESTDIR)$(PREFIX)/share/bash-completion/completions
	install -m 644 source/bashcompletion $(DESTDIR)$(PREFIX)/share/bash-completion/completions/alfa
	gzip -f $(MANDIR)/alfa.1
	@echo "Installation complete"

test: alfa
	@mkdir -p testoutput
	@printf "testing reading and fitting of data formats:\n"
	@printf " ..1d ascii"
	@./alfa test/1d.ascii --output-dir testoutput && printf "....success!\n"
	@printf " ..1d fits image"
	@./alfa test/1d.fits --output-dir testoutput && printf "....success!\n"
	@printf " ..2d fits image (also tests image section handling)"
	@./alfa test/2d_small.fits[*,1:2] --output-dir testoutput && printf "....success!\n"
	@printf " ..3d fits image"
	@./alfa test/3dspec_2x2.fits --output-dir testoutput && printf "....success!\n"
	@printf " ..3d fits image collapsed"
	@./alfa test/3dspec_2x2.fits --collapse --output-dir testoutput && printf "....success!\n"
	@printf "all tests ran successfully!\n"
	@rm -rf testoutput

uninstall:
	rm -rf $(DESTDIR)$(PREFIX)/share/alfa
	rm -f $(DESTDIR)$(PREFIX)/bin/alfa
	rm -f $(DESTDIR)$(PREFIX)/share/bash-completion/completions/alfa
	rm -f $(MANDIR)/alfa.1.gz

htmlmanual:
	groff -m mandoc -Thtml man/alfa.1 > man/alfa.html
