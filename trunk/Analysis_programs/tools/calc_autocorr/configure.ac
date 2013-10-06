AC_PREREQ([2.69])
AC_INIT([calc_autocorr],[1.0],[francesco.sanfilippo@th.u-psud.fr])
AC_CONFIG_SRCDIR([src])
AC_CONFIG_AUX_DIR([config])
AM_CONFIG_HEADER([src/config.hpp])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([1.0 no-define foreign subdir-objects])

# Checks for programs.
AC_PROG_CXX
AC_PROG_INSTALL
AC_PROG_MAKE_SET

# Checks for libraries.
AC_CHECK_LIB([fftw3], [fftw_export_wisdom_to_string],[],AC_MSG_ERROR([unable to find fftw3 library]))

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h])

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([sqrt])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT