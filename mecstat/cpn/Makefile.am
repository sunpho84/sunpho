ACLOCAL_AMFLAGS=-I m4

#svn versioning: specify to build sources and include script
BUILT_SOURCES=update_svnversion src/svnversion.hpp
update_svnversion:
	$(top_srcdir)/src/update_svnversion.hpp.sh $(top_srcdir)
.PHONY: update_svnversion
EXTRA_DIST=src/svnversion.hpp $(top_srcdir)/src/update_svnversion.hpp.sh
CLEANFILES=src/svnversion.hpp

bin_PROGRAMS=$(top_builddir)/bin/cpn
__top_builddir__bin_cpn_SOURCES= \
	src/debug.cpp \
	src/endianness.cpp \
	src/field.cpp \
	src/geometry.cpp \
	src/main.cpp \
	src/neighs.cpp \
	src/random.cpp \
	src/simul.cpp \
	src/threads.cpp \
	src/utils.cpp \
	src/vectors.cpp

# \
#	src/system.cpp \
#	src/utils.cpp
