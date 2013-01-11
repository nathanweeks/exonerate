
TESTS = region.test c4.test alignment.test codegen.test optimal.test  \
        layout.test viterbi.test opair.test subopt.test cgutil.test

noinst_PROGRAMS = $(TESTS)

AM_CPPFLAGS = -I$(top_srcdir)/src/struct                \
           -I$(top_srcdir)/src/comparison            \
           -I$(top_srcdir)/src/sequence              \
           -I$(top_srcdir)/src/general               \
           -DSOURCE_ROOT_DIR="\"@source_root_dir@\"" \
           -DGLIB_CFLAGS="\"@glib_cflags@\""         \
           -DHOSTTYPE="\"@host@\""

noinst_HEADERS = region.h c4.h alignment.h codegen.h optimal.h  \
                 layout.h viterbi.h opair.h subopt.h \
                 cgutil.h

region_test_SOURCES = region.test.c region.c

c4_test_SOURCES = c4.test.c c4.c region.c
c4_test_LDADD = $(top_srcdir)/src/struct/matrix.o      \
                $(top_srcdir)/src/general/threadref.o

SEQUENCE_OBJ = $(top_srcdir)/src/struct/sparsecache.o  \
               $(top_srcdir)/src/struct/matrix.o       \
               $(top_srcdir)/src/sequence/sequence.o   \
               $(top_srcdir)/src/sequence/alphabet.o   \
               $(top_srcdir)/src/sequence/splice.o     \
               $(top_srcdir)/src/sequence/translate.o  \
               $(top_srcdir)/src/general/argument.o    \
               $(top_srcdir)/src/general/lineparse.o   \
               -lm

ALIGNMENT_OBJ = $(top_srcdir)/src/comparison/match.o     \
                $(top_srcdir)/src/sequence/codonsubmat.o \
                $(top_srcdir)/src/sequence/submat.o      \
                $(top_srcdir)/src/general/threadref.o    \
                $(SEQUENCE_OBJ)
alignment_test_SOURCES = alignment.test.c alignment.c c4.c region.c
alignment_test_LDADD = $(ALIGNMENT_OBJ)

codegen_test_SOURCES = codegen.test.c codegen.c
codegen_test_LDADD = $(top_srcdir)/src/general/argument.o

cgutil_test_SOURCES = cgutil.test.c cgutil.c codegen.c c4.c region.c
cgutil_test_LDADD = $(top_srcdir)/src/general/argument.o  \
                    $(top_srcdir)/src/struct/matrix.o     \
                    $(top_srcdir)/src/general/threadref.o

viterbi_test_SOURCES = viterbi.test.c viterbi.c c4.c region.c \
                       alignment.c codegen.c layout.c subopt.c \
                       cgutil.c
viterbi_test_LDADD = $(top_srcdir)/src/struct/slist.o      \
                     $(top_srcdir)/src/struct/recyclebin.o \
                     $(top_srcdir)/src/struct/rangetree.o  \
                     $(ALIGNMENT_OBJ)

optimal_test_SOURCES = optimal.test.c optimal.c c4.c alignment.c \
                       codegen.c region.c layout.c viterbi.c     \
                       subopt.c cgutil.c
optimal_test_LDADD = $(top_srcdir)/src/struct/slist.o      \
                     $(top_srcdir)/src/struct/recyclebin.o \
                     $(top_srcdir)/src/struct/rangetree.o  \
                     $(ALIGNMENT_OBJ)

layout_test_SOURCES = layout.test.c layout.c
layout_test_LDADD = $(top_srcdir)/src/struct/matrix.o

opair_test_SOURCES = opair.test.c opair.c optimal.c c4.c alignment.c \
                     codegen.c region.c layout.c viterbi.c subopt.c  \
                     cgutil.c
opair_test_LDADD = $(top_srcdir)/src/struct/slist.o      \
                   $(top_srcdir)/src/struct/recyclebin.o \
                   $(top_srcdir)/src/struct/rangetree.o  \
                   $(ALIGNMENT_OBJ)

subopt_test_SOURCES = subopt.test.c subopt.c region.c
subopt_test_LDADD = $(top_srcdir)/src/struct/rangetree.o \
                    $(top_srcdir)/src/struct/recyclebin.o


# Files to clear away

MAINTAINERCLEANFILES = Makefile.in

