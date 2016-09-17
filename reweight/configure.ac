AC_INIT([reweight],[1.0],[fr.sanfilippo@gmail.com])

AC_CONFIG_SRCDIR([src])
AC_CONFIG_AUX_DIR([config])
AM_INIT_AUTOMAKE([1.0 no-define foreign subdir-objects])
AM_CONFIG_HEADER(src/config.hpp)
AC_CONFIG_MACRO_DIR([m4])
AC_CANONICAL_HOST

AC_PROG_CXX
AC_CONFIG_FILES(Makefile)

#good-generator
AC_ARG_ENABLE(good-generator,
        AS_HELP_STRING([--enable-good-generator],[Enable costly generation of random variables]),
        enable_good_generator="${enableval}",
        enable_good_generator="no")
if test "$enable_good_generator" == "yes";then AC_DEFINE([GOOD_GENERATOR],1,[Enable good generator]) fi
AC_MSG_RESULT([enabling good generator... $enable_good_generator])

AC_OUTPUT
