ACLOCAL_AMFLAGS=-I m4

#svn versioning: specify to build sources and include script
BUILT_SOURCES=update_svnversion svnversion.hpp
update_svnversion:
	$(srcdir)/src/update_svnversion.hpp.sh $(top_srcdir)
.PHONY: update_svnversion
EXTRA_DIST=svnversion.hpp update_svnversion.hpp.sh
CLEANFILES=svnversion.hpp

bin_PROGRAMS=$(top_builddir)/bin/cpn
__top_builddir__bin_cpn_SOURCES= \
	src/debug.cpp \
	src/fields.cpp \
	src/geometry.cpp \
	src/main.cpp \
	src/simul.cpp \
	src/threads.cpp \
	src/utils.cpp \
	src/vectors.cpp

# \
#	src/random.cpp \
#	src/system.cpp \
#	src/utils.cpp
