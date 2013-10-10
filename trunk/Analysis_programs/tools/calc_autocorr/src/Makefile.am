#include m4 dir
ACLOCAL_AMFLAGS=-I m4

AUTOMAKE_OPTIONS=subdir-objects

#include fftw3
LIBS=-lfftw3

#parser
AM_YFLAGS=-d -r all

bin_PROGRAMS=$(top_builddir)/bin/main
__top_builddir__bin_main_SOURCES= \
	debug.cpp \
	debug.hpp \
        driver.cpp \
        driver.hpp \
        load.cpp \
        load.hpp \
	macros.hpp \
	main.cpp \
	new_types.cpp \
	new_types.hpp \
        parse.ypp \
        tokenize.lpp \
        redefine_yy.hpp

CLEANFILES=parse.cpp parse.hpp tokenize.cpp
