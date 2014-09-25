AC_INIT([cpn],[1.0],[fr.sanfilippo@gmail.com])

AC_CONFIG_SRCDIR([src])
AC_CONFIG_AUX_DIR([config])
AM_INIT_AUTOMAKE([1.0 no-define foreign subdir-objects])
AM_CONFIG_HEADER(src/config.hpp)
AC_CONFIG_MACRO_DIR([m4])
AC_CANONICAL_HOST

AC_PROG_CXX
AC_CONFIG_FILES(Makefile)

PKG_CHECK_MODULES([pkg_fftw3], [fftw3])

#good-generator
AC_ARG_ENABLE(good-generator,
        AS_HELP_STRING([--enable-good-generator],[Enable costly generation of random variables]),
        enable_good_generator="${enableval}",
        enable_good_generator="no")
if test "$enable_good_generator" == "yes";then AC_DEFINE([GOOD_GENERATOR],1,[Enable good generator]) fi
AC_MSG_RESULT([enabling good generator... $enable_good_generator])

#smart-extraction
AC_ARG_ENABLE(smart-extraction,
        AS_HELP_STRING([--enable-smart-extraction],[Enable smart extraction of random variables]),
        enable_smart_extraction="${enableval}",
        enable_smart_extraction="yes")
if test "$enable_smart_extraction" == "yes";then AC_DEFINE([SMART_EXTRACTION],1,[Enable smart extraction]) fi
AC_MSG_RESULT([enabling smart extraction... $enable_smart_extraction])

#debug hmc
AC_ARG_ENABLE(debug-hmc,
        AS_HELP_STRING([--enable-debug-hmc],[Enable debugging of hmc]),
        enable_debug_hmc="${enableval}",
        enable_debug_hmc="no")
if test "$enable_debug_hmc" == "yes";then AC_DEFINE([DEBUG_HMC],1,[Enable debug hmc]) fi
AM_CONDITIONAL([DEBUG_HMC],[test "$enable_debug_hmc" == "yes" ],[true],[false])
AC_MSG_RESULT([enabling debug hmc... $enable_debug_hmc])

AC_OUTPUT
