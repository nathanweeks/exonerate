
TESTS = heuristic.test bsdp.test hpair.test sar.test

noinst_PROGRAMS = $(TESTS)

AM_CPPFLAGS = -I$(top_srcdir)/src/struct                \
           -I$(top_srcdir)/src/comparison            \
           -I$(top_srcdir)/src/sequence              \
           -I$(top_srcdir)/src/general               \
           -I$(top_srcdir)/src/c4                    \
           -DSOURCE_ROOT_DIR="\"@source_root_dir@\"" \
           -DGLIB_CFLAGS="\"@glib_cflags@\""         \
           -DHOSTTYPE="\"@host@\""

noinst_HEADERS = heuristic.h bsdp.h hpair.h sar.h

SEQUENCE_OBJ = $(top_srcdir)/src/struct/sparsecache.o  \
               $(top_srcdir)/src/struct/matrix.o       \
               $(top_srcdir)/src/sequence/sequence.o   \
               $(top_srcdir)/src/sequence/alphabet.o   \
               $(top_srcdir)/src/sequence/splice.o     \
               $(top_srcdir)/src/sequence/translate.o  \
               $(top_srcdir)/src/general/argument.o    \
               $(top_srcdir)/src/general/lineparse.o   \
               $(top_srcdir)/src/general/threadref.o   \
               -lm

C4_OBJ = $(top_srcdir)/src/c4/c4.o        \
         $(top_srcdir)/src/c4/alignment.o \
         $(top_srcdir)/src/c4/optimal.o   \
         $(top_srcdir)/src/c4/codegen.o   \
         $(top_srcdir)/src/c4/region.o    \
         $(top_srcdir)/src/c4/layout.o    \
         $(top_srcdir)/src/c4/viterbi.o   \
         $(top_srcdir)/src/c4/subopt.o    \
         $(top_srcdir)/src/c4/cgutil.o

ALIGNMENT_OBJ = $(top_srcdir)/src/comparison/match.o     \
                $(top_srcdir)/src/sequence/codonsubmat.o \
                $(top_srcdir)/src/sequence/submat.o      \
                $(C4_OBJ) $(SEQUENCE_OBJ)

heuristic_test_SOURCES = heuristic.test.c heuristic.c
heuristic_test_LDADD = $(top_srcdir)/src/struct/slist.o      \
                       $(top_srcdir)/src/struct/recyclebin.o \
                       $(top_srcdir)/src/struct/rangetree.o  \
                       $(ALIGNMENT_OBJ)

bsdp_test_SOURCES = bsdp.test.c bsdp.c
bsdp_test_LDADD = $(top_srcdir)/src/general/argument.o  \
                  $(top_srcdir)/src/struct/pqueue.o     \
                  $(top_srcdir)/src/struct/recyclebin.o

hpair_test_SOURCES = hpair.test.c hpair.c heuristic.c bsdp.c sar.c
hpair_test_LDADD = $(top_srcdir)/src/struct/slist.o         \
                   $(top_srcdir)/src/struct/pqueue.o        \
                   $(top_srcdir)/src/struct/recyclebin.o    \
                   $(top_srcdir)/src/struct/rangetree.o     \
                   $(top_srcdir)/src/comparison/hspset.o    \
                   $(top_srcdir)/src/comparison/wordhood.o  \
                   $(ALIGNMENT_OBJ)

sar_test_SOURCES = sar.test.c sar.c heuristic.c
sar_test_LDADD = $(top_srcdir)/src/struct/slist.o      \
                 $(top_srcdir)/src/struct/recyclebin.o \
                 $(top_srcdir)/src/struct/rangetree.o  \
                 $(ALIGNMENT_OBJ)

# Files to clear away

MAINTAINERCLEANFILES = Makefile.in

