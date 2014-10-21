ACLOCAL_AMFLAGS=-I m4

LDADD=-lfftw3
LDADD+=-lfftw3_threads

cpn_source_files= \
	src/action.cpp \
        src/close.cpp \
        src/corr.cpp \
        src/data.cpp \
        src/debug.cpp \
        src/geometry.cpp \
        src/hmc.cpp \
        src/init.cpp \
        src/lambda.cpp \
        src/metro.cpp \
        src/micro.cpp \
        src/overheat.cpp \
        src/parameters.cpp \
        src/random.cpp \
        src/staples.cpp \
        src/stout.cpp \
        src/tools.cpp \
        src/topology.cpp \
        src/zeta.cpp


bin_PROGRAMS= \
	$(top_builddir)/bin/main \
	$(top_builddir)/bin/reversibility_check

__top_builddir__bin_main_SOURCES=${cpn_source_files} src/main.cpp
__top_builddir__bin_reversibility_check_SOURCES=${cpn_source_files} src/reversibility_check.cpp
