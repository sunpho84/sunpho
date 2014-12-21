ACLOCAL_AMFLAGS=-I m4

LDADD=-lfftw3
LDADD+=-lfftw3_threads

on_source_files= \
	src/geometry.cpp \
	src/init.cpp \
	src/parameters.cpp \
	src/random.cpp \
	src/tools.cpp

bin_PROGRAMS= \
	$(top_builddir)/bin/main

__top_builddir__bin_main_SOURCES=${on_source_files} src/main.cpp
