# configure.in for the exonerate package
#
# Guy St.C. Slater.  16:11:2000
#
# ( Borrowed some config tricks from gtk+, gnet, xscreensaver, mutt )

AC_INIT([exonerate], [2.4.0])
AC_CONFIG_SRCDIR([src/sequence/sequence.c])
AM_INIT_AUTOMAKE
AC_PROG_CC
AC_PROG_INSTALL
AC_PROG_MAKE_SET
# AM_PROG_LIBTOOL

# OTHER STUFF
AC_DEFINE([_XOPEN_SOURCE],[600])

# Linux lacks _XOPEN_PATH_MAX in limits.h...
AC_CHECK_DECLS([_XOPEN_PATH_MAX],[],[AC_DEFINE([_XOPEN_PATH_MAX],[1024])],\
               [[#include <limits.h>]])
AC_C_CONST
AC_C_INLINE
AC_HEADER_STDC # CHECKS FOR stdlib.h stdarg.h string.h and float.h
AC_CHECK_HEADERS(limits.h errno.h math.h strings.h ctype.h \
                 sys/stat.h sys/types.h time.h unistd.h)

# SET host (REQUIRES config.guess AND config.sub)
AC_CANONICAL_HOST
AC_SUBST(host)

echo "Set host to $host"

# AVOID autoconf WARNINGS ABOUT --datarootdir
AC_SUBST(datarootdir)

custom_guint64_format="llu"

# TRY TO GET TO BUILD ON ALL PLATFORMS BY DEFAULT
echo ""$host
case $host in
   *-apple-darwin*)
   echo "Customising build for OSX"
   C4_AR_INIT="csq"
   C4_AR_APPEND="sq"
   ;;
   *-dec-osf*)
   echo "Customising build for OSF"
   C4_AR="/bin/ar"
   C4_AR_INIT="czr"
   C4_AR_APPEND="zr"
   # use the TRU64 C compiler
   CC="cc"
   GCC= # No longer true (prevents GCC specific CFLAGS being set)
   CFLAGS="-O3 -arch ev6 -fast"
   ;;
   *-irix*)
   echo "Customising build for IRIX"
   C4_AR_INIT="cq"
   C4_AR_APPEND="q"
   ;;
   *-solaris*)
   echo "Customising build for Solaris"
   C4_AR_INIT="cq"
   C4_AR_APPEND="q"
   ;;
   x86_64-*)
   echo "Customising build for X86_64"
   custom_guint64_format="lu"
   ;;
   *)
   echo "*** Unknown build platform: $host ***"
   ;;
esac

AC_SUBST(custom_guint64_format)

# SET SOURCE_ROOT_DIR (REQUIRED FOR viterbi.h)
source_root_dir=`pwd`
AC_SUBST(source_root_dir)
echo "Set source_root_dir to $source_root_dir"

echo "----------------------------------------------------------------"

# CHECK FOR socklen_t
# THIS IS TO FIX PROBLEMS WHERE socklen_t ISN'T DEFINED (eg. ON OSF)
#
AC_MSG_CHECKING(for socklen_t)
AC_TRY_COMPILE([#include <sys/types.h>
#include <sys/socket.h>
socklen_t x;], [], [AC_MSG_RESULT(yes)],
[AC_TRY_COMPILE([#include <sys/types.h>
#include <sys/socket.h>
int accept(int, struct sockaddr *, size_t *);],
[], [AC_MSG_RESULT(size_t)
AC_DEFINE(socklen_t, size_t)], [
AC_MSG_RESULT(int)
AC_DEFINE(socklen_t, int)])])

# ALLOW EITHER GLIB-1 OR GLIB-2 TO BE USED
#
# IT SEEMS TRICKY TO ALLOW COMPILATION AGAINST BOTH GLIB-1 AND GLIB-2
# THAT WILL WORK FOR ONE LIBRARY WHEN THE OTHER IS ABSENT ...
#
# ... WHAT FOLLOWS IS FAR FROM IDEAL, BUT SEEMS TO WORK.
# ANY SUGGESTIONS FOR A BETTER WAY OF DOING THIS WILL BE WELCOME.
#
AC_ARG_ENABLE(glib2,
[
  --enable-glib2  Use glib2 library
  --disable-glib2 Do not glib2 (use glib1 instead)],
[enable_glib2="$enableval"],[enable_glib2=yes])
if test "$enable_glib2" = yes; then
#    AM_PATH_GLIB_2_0(2.0.0,
#            [LIBS="$LIBS $GLIB_LIBS" CFLAGS="$CFLAGS $GLIB_CFLAGS"],
#            AC_MSG_ERROR(Cannot find GLIB2: Is pkg-config in path?))
#    PKG_CHECK_MODULES(GLIB, [glib-2.0], [:], [:])
    AC_PATH_PROG(PKG_CONFIG, pkg-config, no)
    if test "$PKG_CONFIG" = no; then
        echo "ERROR: Could not find pkg-config ... is glib-2 installed ???"
        exit 1
    fi
    echo "Using GLIB-2"
    glib_cflags=`pkg-config --cflags glib-2.0`
    glib_libs=`pkg-config --libs glib-2.0`
    CFLAGS="$CFLAGS $glib_cflags"
    LIBS="$LIBS $glib_libs"
elif test "$enable_glib2" = no; then
    echo "Using GLIB-1"
#    AM_PATH_GLIB([1.2.0],
#            [LIBS="$LIBS $GLIB_LIBS" CFLAGS="$CFLAGS $GLIB_CFLAGS"],
#            AC_MSG_ERROR(Cannot find GLIB1: Is glib-config in path?))
    AC_PATH_PROG(GLIB_CONFIG, glib-config, no)
    if test "$GLIB_CONFIG" = no; then
        echo "ERROR: Could not find glib-config ... is glib-1 installed ???"
        exit 1
    fi
    glib_cflags=`glib-config --cflags glib`
    glib_libs=`glib-config --libs glib`
    CFLAGS="$CFLAGS $glib_cflags"
    LIBS="$LIBS $glib_libs"
else
    echo "error: must be yes or no: --enable-glib2:[$enable_glib2]"
    exit 1
fi

AC_SUBST(glib_libs)
AC_SUBST(glib_cflags)

# ALLOW ASSERTIONS TO BE TURNED ON OR OFF
# THIS ALSO ACTIVATES DEBUGGING OR OPTIMISATION COMPILE OPTIONS
AC_ARG_ENABLE(assert,
[
  --enable-assert        Use code for assertion tests
  --disable-assert       Do not use assertion code],
[enable_assert="$enableval"],[enable_assert=no])
if test "$enable_assert" = yes; then
    echo "Turning ASSERTIONS on"
    CFLAGS="$CFLAGS -g"
elif test "$enable_assert" = no; then
    CFLAGS="$CFLAGS -DG_DISABLE_ASSERT"
    echo "Turning assertions off"
    if test "$GCC" = "yes"; then
        # Not currently using -fomit-frame-pointer as clashes with -pg
        # CFLAGS="$CFLAGS -O3 -fomit-frame-pointer -finline-functions"
        CFLAGS="$CFLAGS -O3 -finline-functions"
    fi
else
    echo "error: must be yes or no: --enable-assert:[$enable_assert]"
    exit 1
fi

# ALLOW PARANOID CODE CHECKING TO BE TURNED ON OR OFF
AC_ARG_ENABLE(paranoia,
[
  --enable-paranoia        Use paranoid compilation options
  --disable-paranoia       Do not use paranoid compilation options],
[enable_paranoia="$enableval"],[enable_paranoia=no])
if test "$enable_paranoia" = yes; then
    echo "Turning PARANOID options on"
    if test "$GCC" = "yes"; then
        CFLAGS="$CFLAGS -Wall -Werror \
                -Wstrict-prototypes -Wmissing-prototypes\
                -Wmissing-declarations"
    fi
elif test "$enable_paranoia" = no; then
    echo "Turning paranoid options off"
else
    echo "error: must be yes or no: --enable-paranoia:[$enable_paranoia]"
    exit 1
fi

# ALLOW PROFILING WITH GPROF
AC_ARG_ENABLE(gprof,
[
  --enable-gprof    Generate executables suitable for gprof profiling
  --disable-gprof   Do not create executables for gprof profiling],
[enable_gprof="$enableval"],[enable_gprof=no])
if test "$enable_gprof" = yes; then
    echo "Building executables for profiling with GPROF"
    echo "Statically linking for gprof"
    echo "    Run: <program> then run: gprof <program> gmon.out"
    CFLAGS="$CFLAGS -O0 -static -pg"
    LDFLAGS="$LDFLAGS -g -pg -nodefaultlibs"
    LIBS="$LIBS -lc_p -lgcc -lm"
elif test "$enable_gprof" = no; then
    echo "Not using options for profiling with gprof"
else
    echo "error: must be yes or no: --enable-gprof:[$enable_gprof]"
    exit 1
fi

# ALLOW COVERAGE TESTING WITH GCOV
AC_ARG_ENABLE(gcov,
[
  --enable-gcov Generate executables suitable for gcov analysis
  --disable-gcov Do not create executables for gcov analysis],
[enable_gcof="$enableval"],[enable_gcov=no])
if test "$enable_gcov" = yes; then
    echo "Building executables for analysis with GCOV"
    echo "    Run <program> then run: gcov <source_file>.c"
    echo "    then view: <source_file>.c.gcov"
    CFLAGS="$CFLAGS -fprofile-arcs -ftest-coverage"
elif test "$enable_gcov" = no; then
    echo "Not using options for analysis with gcov"
else
    echo "error: must be yes or no: --enable-gcov:[$enable_gcov]"
    exit 1
fi

# FORCE LARGEFILE SUPPORT FOR 32 BIT SYSTEMS
AC_ARG_ENABLE(largefile,
[
  --enable-largefile Enable largefile support on 32 bit systems
  --disable-largefile Do not enable largefile support on 32 bit systems],
[enable_largefile="$enableval"],[enable_largefile=yes])
if test "$enable_largefile" = yes; then
    echo "Adding support for LARGE FILES (>2Gb) for a 32 bit system"
    CFLAGS="$CFLAGS -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64"
elif test "$enable_largefile" = no; then
    echo "Not forcing largefile support for 32 bit systems"
else
    echo "error: must be yes or no:"
    echo "    --enable-largefile:[$enable_largefile]"
    exit 1
fi

# ALLOW USE OF COMPILED MODELS
AC_ARG_ENABLE(compiledmodels,
[
  --enable-compiledmodels Use compiled C4 models
  --disable-compiledmodels Do not use compiled C4 models],
[enable_compiledmodels="$enableval"],[enable_compiledmodels=yes])
if test "$enable_compiledmodels" = yes; then
    echo "Using COMPILED C4 MODELS"
    codegen_extra_ldadd="viterbi.o scheduler.o c4_model_archive.a"
    codegen_extra_sources="c4_model_archive.h"
elif test "$enable_compiledmodels" = no; then
    echo "Not using compiled C4 models"
# model_extra_programs="bootstrapper"
    codegen_extra_sources=""
    codegen_extra_ldadd="`pwd`/src/c4/viterbi.o  \
                         `pwd`/src/sdp/scheduler.o"
else
    echo "error: must be yes or no:"
    echo "    --enable-compiledmodels:[$enable_compiledmodels]"
    exit 1
fi
# AC_SUBST(model_extra_programs)
AC_SUBST(codegen_extra_sources)
AC_SUBST(codegen_extra_ldadd)

# ALLOW INSTALLATION OF UTILITIES
AC_ARG_ENABLE(utilities,
[
  --enable-utilities Install all utilities
  --disable-utilities Do not install utilities],
[enable_utilities="$enableval"],[enable_utilities=yes])
if test "$enable_utilities" = yes; then
    echo "Installing UTILITIES"
    installed_util_list="                                      \
         fastaclean fastaclip fastachecksum fastacomposition   \
         fastadiff fastaexplode fastafetch fastahardmask       \
         fastaindex fastalength fastanrdb fastaoverlap         \
         fastareformat fastaremove fastarevcomp fastasoftmask  \
         fastasort fastasplit fastasubseq fastatranslate       \
         fastavalidcds fasta2esd fastaannotatecdna esd2esi"
elif test "$enable_utilities" = no; then
    echo "Not installing utilities"
    installed_util_list=""
else
    echo "error: must be yes or no:"
    echo "    --enable-utilities:[$enable_utilities]"
    exit 1
fi
AC_SUBST(installed_util_list)

# 2012-01-20 (nathanweeks): Disable pthreads code due to reliability issues &
# limited utility.
## ALLOW DISABLING OF PTHREADS
#AC_ARG_ENABLE(pthreads,
#[
#  --enable-pthreads Use pthreads for exonerate-server
#  --disable-pthreads Do not use pthreads for exonerate-server],
#[enable_pthreads="$enableval"],[enable_pthreads=yes])
#if test "$enable_pthreads" = yes; then
#    echo "Using PTHREADS"
#    CFLAGS="$CFLAGS -DUSE_PTHREADS"
#    # for g_thread_init()
#    g_thread_init_ldflags=`pkg-config --libs gthread-2.0`
#    LDFLAGS="$LDFLAGS $g_thread_init_ldflags -lpthread"
#elif test "$enable_pthreads" = no; then
#    echo "Not using pthreads"
#else
#    echo "error: must be yes or no:"
#    echo "    --enable-pthreads:[$enable_pthreads]"
#    exit 1
#fi
#AC_SUBST(use_pthreads)

echo "----------------------------------------------------------------"

AC_OUTPUT(Makefile                   \
          doc/Makefile               \
          doc/man/Makefile           \
          doc/man/man1/Makefile      \
          src/Makefile               \
          src/general/Makefile       \
          src/sequence/Makefile      \
          src/database/Makefile      \
          src/struct/Makefile        \
          src/comparison/Makefile    \
          src/c4/Makefile            \
          src/bsdp/Makefile          \
          src/sdp/Makefile           \
          src/hub/Makefile           \
          src/model/Makefile         \
          src/program/Makefile       \
          src/util/Makefile          \
          test/Makefile              \
          test/util/Makefile         \
          test/ipcress/Makefile      \
          test/exonerate/Makefile    \
          test/data/Makefile         \
          test/data/cdna/Makefile    \
          test/data/protein/Makefile )

