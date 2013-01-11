
TESTS = fastadb.test fastapipe.test dataset.test index.test
# seqfmi.test
noinst_PROGRAMS = $(TESTS)

AM_CPPFLAGS = -I$(top_srcdir)/src/sequence   \
           -I$(top_srcdir)/src/struct     \
           -I$(top_srcdir)/src/general    \
           -I$(top_srcdir)/src/comparison \
           -DCUSTOM_GUINT64_FORMAT="\"@custom_guint64_format@\""

noinst_HEADERS = fastadb.h fastapipe.h dataset.h index.h
# seqfmi.h

SEQUENCE_OBJ = $(top_srcdir)/src/struct/sparsecache.o  \
               $(top_srcdir)/src/struct/matrix.o       \
               $(top_srcdir)/src/sequence/sequence.o   \
               $(top_srcdir)/src/sequence/alphabet.o   \
               $(top_srcdir)/src/sequence/splice.o     \
               $(top_srcdir)/src/sequence/translate.o  \
               $(top_srcdir)/src/general/argument.o    \
               $(top_srcdir)/src/general/lineparse.o   \
               -lm
fastadb_test_SOURCES = fastadb.test.c fastadb.c
fastadb_test_LDADD = $(top_srcdir)/src/general/compoundfile.o \
                     $(SEQUENCE_OBJ)

fastapipe_test_SOURCES = fastapipe.test.c fastapipe.c
fastapipe_test_LDADD = fastadb.o                              \
                     $(top_srcdir)/src/general/compoundfile.o \
                     $(SEQUENCE_OBJ)

dataset_test_SOURCES = dataset.test.c dataset.c
dataset_test_LDADD = fastadb.o                                \
                     $(top_srcdir)/src/general/compoundfile.o \
                     $(top_srcdir)/src/struct/bitarray.o      \
                     $(SEQUENCE_OBJ)

index_test_SOURCES = index.test.c index.c
index_test_LDADD = fastadb.o dataset.o                        \
                   $(top_srcdir)/src/general/compoundfile.o   \
                   $(top_srcdir)/src/general/threadref.o      \
                   $(top_srcdir)/src/struct/bitarray.o        \
                   $(top_srcdir)/src/struct/vfsm.o            \
                   $(top_srcdir)/src/struct/pqueue.o          \
                   $(top_srcdir)/src/struct/recyclebin.o      \
                   $(top_srcdir)/src/struct/rangetree.o       \
                   $(top_srcdir)/src/struct/splaytree.o       \
                   $(top_srcdir)/src/struct/noitree.o         \
                   $(top_srcdir)/src/comparison/wordhood.o    \
                   $(top_srcdir)/src/comparison/hspset.o      \
                   $(top_srcdir)/src/comparison/match.o       \
                   $(top_srcdir)/src/sequence/submat.o        \
                   $(top_srcdir)/src/sequence/codonsubmat.o   \
                   $(SEQUENCE_OBJ)

# seqfmi_test_SOURCES = seqfmi.test.c seqfmi.c
# seqfmi_test_LDADD = fastadb.o                                  \
#                     $(top_srcdir)/src/general/compoundfile.o   \
#                     $(top_srcdir)/src/struct/fmindex.o         \
#                     $(top_srcdir)/src/struct/bitarray.o        \
#                     $(top_srcdir)/src/struct/pqueue.o          \
#                     $(top_srcdir)/src/struct/recyclebin.o      \
#                     $(top_srcdir)/src/struct/huffman.o         \
#                     $(SEQUENCE_OBJ)

# Files to clear away

MAINTAINERCLEANFILES = Makefile.in

