AC_INIT([cpn],[1.0],[fr.sanfilippo@gmail.com])

AC_CONFIG_SRCDIR([src])
AC_CONFIG_AUX_DIR([config])
AM_INIT_AUTOMAKE([1.0 no-define foreign subdir-objects])
AM_CONFIG_HEADER(src/config.hpp)
AC_CONFIG_MACRO_DIR([m4])
AC_CANONICAL_HOST

#silent automake
AM_SILENT_RULES([yes])

AC_PROG_CXX

AX_CXX_COMPILE_STDCXX_11(noext,mandatory)

AC_CONFIG_FILES(Makefile)

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
