#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([MetaSMC], [1.0], [fjliao@gmail.com])
AM_INIT_AUTOMAKE(MetaSMC, 1.0)

AC_CONFIG_SRCDIR([src/])
AC_CONFIG_HEADERS([config.h])

AC_ENABLE_SHARED
AC_DISABLE_STATIC
LT_INIT

# Checks for programs.
AC_PROG_CC
PKG_CHECK_MODULES(LAPACK, lapack)

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lgfortran':
AC_CHECK_LIB([gfortran], [main], [], [AC_MSG_ERROR([-lgfortran not found])])
# FIXME: Replace `main' with a function in `-llapack':
AC_CHECK_LIB([lapack], [main], [], [AC_MSG_ERROR([-llapack not found])])
# FIXME: Replace `main' with a function in `-llapacke':
AC_CHECK_LIB([lapacke], [LAPACKE_dgeev], [], [AC_MSG_ERROR([LAPACKE_dgeev of -llapacke not found])])
AC_CHECK_LIB([lapacke], [LAPACKE_dgetrf], [], [AC_MSG_ERROR([LAPACKE_dgetrf of -llapacke not found])])
AC_CHECK_LIB([lapacke], [LAPACKE_dgetri], [], [AC_MSG_ERROR([LAPACKE_dgetri of -llapacke not found])])

# FIXME: Replace `main' with a function in `-lm':
AC_CHECK_LIB([m], [main], [], [AC_MSG_ERROR([-lm not found])])
# FIXME: Replace `main' with a function in `-lrt':
AC_CHECK_LIB([rt], [main])

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h string.h unistd.h])
#AC_CHECK_HEADERS(LAPACK_FLAGS/lapacke.h, [], [AC_MSG_ERROR(lapacke.h not found)])

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_FUNC_STRTOD
AC_CHECK_FUNCS([clock_gettime memset sqrt strtol])

AC_CONFIG_FILES([Makefile
                 src/Makefile
                 src/mt19937-64/Makefile])
AC_OUTPUT
