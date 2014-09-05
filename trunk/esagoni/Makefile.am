ACLOCAL_AMFLAGS=-I m4

bin_PROGRAMS=$(top_builddir)/bin/esa
__top_builddir__bin_esa_SOURCES= \
        src/esa.cpp \
        src/pos.cpp \
        src/tile.cpp

if PLOT
__top_builddir__bin_esa_SOURCES+=src/plot.cpp
endif
