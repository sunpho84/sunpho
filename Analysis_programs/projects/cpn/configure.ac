AC_INIT([analysis],[1.0],[fr.sanfilippo@gmail.com])

AC_CONFIG_SRCDIR([src])
AC_CONFIG_AUX_DIR([config])
AM_INIT_AUTOMAKE([1.0 no-define foreign subdir-objects])
AM_CONFIG_HEADER(src/config.hpp)
AC_CONFIG_MACRO_DIR([m4])
AC_CANONICAL_HOST

#silent automake
AM_SILENT_RULES([yes])

AC_PROG_CC
AC_PROG_CXX

#Checks for headers
AC_CHECK_HEADERS([fftw3.h],[fftw_found_headers=yes])

# Checks for libraries
AC_SEARCH_LIBS(fftw_plan_dft_r2c_2d,fftw3)

AX_CXX_COMPILE_STDCXX_11(noext,mandatory)

AC_CONFIG_FILES(Makefile)


AC_OUTPUT
