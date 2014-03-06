AC_INIT([mecstat],[1.0],[fr.sanfilippo@gmail.com])

#include all config flags
config_flags="$*"
AC_DEFINE_UNQUOTED([CONFIG_FLAGS],["$config_flags"],[Flags passed to configure])

#get configuration time
AC_MSG_CHECKING([configuration date and time])
CONFIG_DATE=`date`
AC_DEFINE_UNQUOTED([CONFIG_TIME], "$CONFIG_TIME", [Time when configured])
AC_MSG_RESULT([$CONFIG_DATE])

AC_CONFIG_SRCDIR([src])
AC_CONFIG_AUX_DIR([config])
AM_INIT_AUTOMAKE([1.0 no-define foreign subdir-objects])
AM_CONFIG_HEADER(src/config.hpp)
AC_CONFIG_MACRO_DIR([m4])
AC_CANONICAL_HOST

AC_PROG_CXX
AC_CHECK_HEADERS([stdlib.h])
AC_HEADER_STDBOOL
AC_CONFIG_FILES(Makefile)

AC_OUTPUT
