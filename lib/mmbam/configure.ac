# Perlude.

AC_INIT([mmbam], [0.1.0], [atrusqiao@gmail.com])
AC_PREREQ([2.69])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([m4])

## some common options
AX_CHECK_ENABLE_DEBUG()
AM_CONDITIONAL([DEBUG], [test x"$ax_enable_debug" = x"yes"])
# set log level, meant to be used with YiCppLib::Log
AM_COND_IF([DEBUG], [AC_DEFINE([LOGLVL], [LOGLV_DEBUG], [Compile time log level])])

AC_PROG_CC
AC_PROG_CXX

# C++14
AX_CXX_COMPILE_STDCXX_14

# Automake
AM_INIT_AUTOMAKE([-Wall -Werror foreign subdir-objects])

AM_PROG_AR

# Libtool
LT_PREREQ([2.4])
LT_INIT()

AX_PTHREAD

# PKG_CHECK_MODULES([GTEST], [gtest_main])
AC_CHECK_LIB([deflate], [libdeflate_gzip_decompress_ex], [], [AC_MSG_ERROR([libdeflate not found])])
AC_CHECK_LIB([tbb], [main], [], [AC_MSG_ERROR([libtbb not found])])

# Config files to output
# Makefile: top level makefile
# src/Makefile: actual targets for the project
# test/Makefile: any test targets for the project
# contrib/Makefile: how to handle 3rd party projects
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT
