AC_INIT([cpn],[1.0],[fr.sanfilippo@gmail.com])

AC_CONFIG_SRCDIR([src])
AC_CONFIG_AUX_DIR([config])
AM_INIT_AUTOMAKE([1.0 no-define foreign subdir-objects])
AM_CONFIG_HEADER(src/config.hpp)
AC_CONFIG_MACRO_DIR([m4])
AC_CANONICAL_HOST

AC_PROG_CXX
AC_PROG_CC
AC_CONFIG_FILES(Makefile)

#debug hmc
AC_ARG_ENABLE(plot,
        AS_HELP_STRING([--enable-plot],[Enable plotting (require root)]),
        enable_plot="${enableval}",
        enable_plot="no")
if test "$enable_plot" == "yes";then AC_DEFINE([PLOT],1,[Enable plotting]) fi
AM_CONDITIONAL([PLOT],[test "$enable_plot" == "yes" ],[true],[false])
AC_MSG_RESULT([enabling plotting... $enable_plot])

AC_OUTPUT
