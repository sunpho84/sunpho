#include m4 dir
ACLOCAL_AMFLAGS=-I m4

AUTOMAKE_OPTIONS=subdir-objects

#include fftw3
LIBS=-lfftw3

#parser
AM_YFLAGS=-d

bin_PROGRAMS=$(top_builddir)/bin/main
__top_builddir__bin_main_SOURCES= \
	macros.hpp \
	debug.cpp \
	debug.hpp \
	main.cpp \
        parse.ypp \
        tokenize.lpp \
        parse.hpp \
        driver.hpp \
        redefine_yy.hpp
