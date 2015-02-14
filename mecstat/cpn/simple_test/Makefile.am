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
        src/polyakov.cpp \
        src/random.cpp \
        src/staples.cpp \
        src/stout.cpp \
        src/tools.cpp \
        src/topology.cpp \
        src/zeta.cpp


include_HEADERS= \
	src/action.hpp \
        src/close.hpp \
        src/corr.hpp \
        src/cpn.hpp \
        src/data.hpp \
        src/debug.hpp \
        src/geometry.hpp \
        src/hmc.hpp \
        src/init.hpp \
        src/lambda.hpp \
        src/macros.hpp \
        src/metro.hpp \
        src/micro.hpp \
        src/overheat.hpp \
        src/parameters.hpp \
        src/polyakov.hpp \
        src/random.hpp \
        src/routines.hpp \
        src/sitmo_prng_engine.hpp \
        src/staples.hpp \
        src/stout.hpp \
        src/tools.hpp \
        src/types.hpp \
        src/topology.hpp \
        src/zeta.hpp


bin_PROGRAMS= \
	$(top_builddir)/bin/main \
	$(top_builddir)/bin/reversibility_check

__top_builddir__bin_main_SOURCES=${cpn_source_files} src/main.cpp
__top_builddir__bin_reversibility_check_SOURCES=${cpn_source_files} src/reversibility_check.cpp
