ACLOCAL_AMFLAGS=-I m4

AM_CPPFLAGS=-I$(top_srcdir)/src

include_HEADERS=src/bissa.hpp
lib_LIBRARIES=$(top_builddir)/lib/libbissa.a

############################################# programs #####################################

LDADD=$(top_builddir)/lib/libbissa.a

bin_PROGRAMS= \
        $(top_builddir)/bin/test

__top_builddir__bin_test_SOURCES=projects/test/main.cpp

########################################### base sources ##################################

__top_builddir__lib_libbissa_a_SOURCES= \
	src/base/close.cpp \
	src/base/debug.cpp \
	src/base/global_variables.cpp \
	src/base/init.cpp \
	src/base/linalgs.cpp \
	src/base/random.cpp \
	src/base/vectors.cpp
include_HEADERS+= \
	src/base/close.hpp \
	src/base/debug.hpp \
	src/base/global_variables.hpp \
	src/base/init.hpp \
	src/base/linalgs.hpp \
	src/base/macros.hpp \
	src/base/random.hpp \
	src/base/thread_macros.hpp \
	src/base/vectors.hpp

#svn versioning: specify to build sources and include script
BUILT_SOURCES=update_svnversion svnversion.hpp
update_svnversion:
	$(srcdir)/src/base/update_svnversion.hpp.sh $(top_srcdir)
.PHONY: update_svnversion
EXTRA_DIST=svnversion.hpp src/base/update_svnversion.hpp.sh

#include global_variables.cpp in the installation
#needed because included in global_variables.hpp
global_variabledir=$(prefix)/include/base
global_variable_DATA=src/base/global_variables.cpp

#clean svnversion.hpp
CLEANFILES=base/svnversion.hpp

######################################### communicate ######################################

__top_builddir__lib_libbissa_a_SOURCES+= \
	src/communicate/all_to_all.cpp \
	src/communicate/borders.cpp

include_HEADERS+= \
	src/communicate/all_to_all.hpp \
	src/communicate/borders.hpp \
	src/communicate/communicate.hpp

################################### geometry ###################################

__top_builddir__lib_libbissa_a_SOURCES+= \
	src/geometry/geometry_eo.cpp \
	src/geometry/geometry_lx.cpp

include_HEADERS+= \
	src/geometry/geometry_eo.hpp \
	src/geometry/geometry_lx.hpp

############################################### io #################################

__top_builddir__lib_libbissa_a_SOURCES+= \
	src/io/checksum.cpp \
	src/io/endianness.cpp \
	src/io/input.cpp \
	src/io/ILDG_File.cpp \
	src/io/reader.cpp \
	src/io/writer.cpp

include_HEADERS+= \
	src/io/checksum.hpp \
	src/io/endianness.hpp \
	src/io/input.hpp \
	src/io/ILDG_File.hpp \
	src/io/reader.hpp \
	src/io/writer.hpp

########################################## new types #################################

#__top_builddir__lib_libbissa_a_SOURCES+=

include_HEADERS+= \
	src/new_types/new_types_definitions.hpp

########################################## operations ###################################

__top_builddir__lib_libbissa_a_SOURCES+= \
	src/operations/remap_vector.cpp

include_HEADERS+= \
	src/operations/remap_vector.hpp

########################################## routines ###################################

__top_builddir__lib_libbissa_a_SOURCES+= \
	src/routines/ios.cpp \
	src/routines/math_routines.cpp \
	src/routines/mpi_routines.cpp

include_HEADERS+= \
	src/routines/ios.hpp \
	src/routines/math_routines.hpp \
	src/routines/mpi_routines.hpp

#compile or not thread support
if USE_THREADS
__top_builddir__lib_libbissa_a_SOURCES+=src/routines/thread.cpp
include_HEADERS+=src/routines/thread.hpp
endif
