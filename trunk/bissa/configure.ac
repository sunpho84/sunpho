AC_INIT([bissa],[1.0],[fr.sanfilippo@gmail.com])

#include all config flags
config_flags="$*"
AC_DEFINE_UNQUOTED([CONFIG_FLAGS],["$config_flags"],[Flags passed to configure])

#get configuration time
AC_MSG_CHECKING([configuration date and time])
CONFIG_DATE=`date`
AC_DEFINE_UNQUOTED([CONFIG_TIME], "$CONFIG_TIME", [time when configured])
AC_MSG_RESULT([$CONFIG_DATE])

AC_CONFIG_SRCDIR([src])
AC_CONFIG_AUX_DIR([config])
AM_INIT_AUTOMAKE([1.0 no-define foreign subdir-objects])
AM_CONFIG_HEADER(src/config.hpp)
AC_CONFIG_MACRO_DIR([m4])
AC_CANONICAL_HOST

AC_PROG_CC
AC_PROG_CXX
AC_PROG_RANLIB

#MPI I/O
AC_ARG_ENABLE(mpi-io,
        AS_HELP_STRING([--disable-mpi-io],[Enable MPI I/O]),
        enable_mpi_io="${enableval}",
        enable_mpi_io="yes")
if test "$enable_mpi_io" == "yes";then AC_DEFINE([USE_MPI_IO],1,[Enable MPI I/O]) fi
AM_CONDITIONAL([USE_MPI_IO],[test "$enable_mpi_io" == "yes" ],[true],[false])
AC_MSG_RESULT([enabling MPI I/O... $enable_mpi_io])

#threads
AC_ARG_ENABLE(threads,
        AS_HELP_STRING([--disable-threads],[Enable threads]),
        enable_threads="${enableval}",
        enable_threads="yes")
if test "$enable_threads" == "yes";then AC_DEFINE([USE_THREADS],1,[Enable threads]) fi
AM_CONDITIONAL([USE_THREADS],[test "$enable_threads" == "yes" ],[true],[false])
AC_MSG_RESULT([enabling threads... $enable_threads])

#MPI
AC_ARG_ENABLE(MPI,
        AS_HELP_STRING([--disable-MPI],[Enable MPI]),
        enable_MPI="${enableval}",
        enable_MPI="yes")
if test "$enable_MPI" == "yes";then AC_DEFINE([USE_MPI],1,[Enable MPI]) fi
AM_CONDITIONAL([USE_MPI],[test "$enable_MPI" == "yes" ],[true],[false])
AC_MSG_RESULT([enabling MPI... $enable_MPI])

#thread debug
AC_ARG_ENABLE(thread-debug,
        AS_HELP_STRING([--enable-thread-debug],[Enable thread debugging]),
        enable_thread_debug="${enableval}",
        enable_thread_debug="no")
if test "$enable_thread_debug" = "yes";then AC_DEFINE([THREAD_DEBUG],1,[Enable thread_debug]) fi
if test "$enable_thread_debug" = "2";then AC_DEFINE([THREAD_DEBUG],2,[Enable thread_debug]) fi
AM_CONDITIONAL([THREAD_DEBUG],[test "$enable_thread_debug" == "yes" ],[true],[false])
AC_MSG_RESULT([enabling thread debug... $enable_thread_debug])

#max verbosity lv
AC_ARG_ENABLE(max-verbosity-lv,
        AS_HELP_STRING([--enable-max-verbosity-lv=lv],[Set max verbosity lv]),
        max_verbosity_lv="${enableval}",
        max_verbosity_lv=2)
if ! printf "%f" "$max_verbosity_lv" >/dev/null 2>&1; then AC_MSG_ERROR([String $max_verbosity_lv is not a number]); fi
AC_MSG_RESULT([enabling max verbosity lv... $max_verbosity_lv])
AC_DEFINE_UNQUOTED([MAX_VERBOSITY_LV],$max_verbosity_lv,[Max_verbosity_lv])

AC_CONFIG_FILES(Makefile)

AC_OUTPUT
