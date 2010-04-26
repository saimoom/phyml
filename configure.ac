#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
dnl AC_INIT([PhyML(beta)],m4_esyscmd([sh -c "head -n 1 version | sed -e 's|^|v3.0-|g' | sed -e 's|$||g' | tr -d '\n'"]), [guindon@stat.auckland.ac.nz])
AC_INIT([PhyML],[beta],[guindon@stat.auckland.ac.nz])
dnl AC_CONFIG_SRCDIR([simu.c])
AC_CONFIG_SRCDIR([src/simu.c])
AC_CONFIG_HEADER([config.h])

AM_INIT_AUTOMAKE
AC_CANONICAL_HOST

# Checks for programs.
AC_PROG_CC

# Checks for libraries.
AC_CHECK_LIB([m], [log])

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([float.h stdlib.h string.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_FUNC_SETVBUF_REVERSED
AC_FUNC_VPRINTF
AC_CHECK_FUNCS([floor pow rint sqrt strchr strstr])

# Thanks to Shane Sturrock from Biomatters
dnl case "${host}" in
dnl *darwin*)
dnl AM_INIT_AUTOMAKE([no-dependencies])
dnl ARCH_flag="-arch i386 -arch ppc -arch x86_64 -mmacosx-version-min=10.4"
dnl ;;
dnl esac


CFLAGS="-O3 -fomit-frame-pointer -funroll-loops ${ARCH_flag}"
dnl CFLAGS="-O3 -Wfloat-equal -fomit-frame-pointer -funroll-loops"
dnl CFLAGS="-O3 -Wfloat-equal"
dnl CFLAGS="-g"

LDFLAGS="${ARGC_flag}"

AC_PROG_LIBTOOL

AC_ARG_ENABLE([tiporder],[  --enable-tiporder    Compile tiporder],[tiporder=yes])
AM_CONDITIONAL([WANT_TIPORDER], [test "$tiporder" = yes]) 

AC_ARG_ENABLE([times],[  --enable-times    Compile PhyTime],[times=yes])
AM_CONDITIONAL([WANT_TIMES], [test "$times" = yes]) 

AC_ARG_ENABLE([part],[  --enable-part    Compile Part],[part=yes])
AM_CONDITIONAL([WANT_PART], [test "$part" = yes]) 

AC_ARG_ENABLE([rwrap],[  --enable-rwrap    Compile Rwrap],[rwrap=yes])
AM_CONDITIONAL([WANT_RWRAP], [test "$rwrap" = yes]) 

AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT

