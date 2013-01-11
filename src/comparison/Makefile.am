
TESTS = wordhood.test hspset.test pcr.test \
        comparison.test match.test seeder.test

noinst_PROGRAMS = $(TESTS)

AM_CPPFLAGS = -I$(top_srcdir)/src/sequence \
           -I$(top_srcdir)/src/struct   \
           -I$(top_srcdir)/src/general

noinst_HEADERS = wordhood.h hspset.h pcr.h \
                 comparison.h match.h seeder.h

SEQUENCE_OBJ = $(top_srcdir)/src/struct/sparsecache.o  \
               $(top_srcdir)/src/struct/matrix.o       \
               $(top_srcdir)/src/sequence/sequence.o   \
               $(top_srcdir)/src/sequence/alphabet.o   \
               $(top_srcdir)/src/sequence/splice.o     \
               $(top_srcdir)/src/sequence/translate.o  \
               $(top_srcdir)/src/general/argument.o    \
               $(top_srcdir)/src/general/lineparse.o   \
               -lm
wordhood_test_SOURCES = wordhood.test.c wordhood.c
wordhood_test_LDADD = $(top_srcdir)/src/sequence/submat.o      \
                      $(top_srcdir)/src/sequence/codonsubmat.o \
                      $(SEQUENCE_OBJ)
hspset_test_SOURCES = hspset.test.c hspset.c match.c wordhood.c
hspset_test_LDADD = $(top_srcdir)/src/sequence/submat.o      \
                    $(top_srcdir)/src/sequence/codonsubmat.o \
                    $(top_srcdir)/src/struct/pqueue.o        \
                    $(top_srcdir)/src/struct/recyclebin.o    \
                    $(top_srcdir)/src/general/threadref.o    \
                    $(SEQUENCE_OBJ)

pcr_test_SOURCES = pcr.test.c pcr.c wordhood.c
pcr_test_LDADD = $(top_srcdir)/src/sequence/submat.o      \
                 $(top_srcdir)/src/sequence/codonsubmat.o \
                 $(top_srcdir)/src/struct/fsm.o           \
                 $(top_srcdir)/src/struct/slist.o         \
                 $(top_srcdir)/src/struct/recyclebin.o    \
                 $(SEQUENCE_OBJ)

comparison_test_SOURCES = comparison.test.c comparison.c hspset.c \
                          match.c wordhood.c
comparison_test_LDADD = $(top_srcdir)/src/sequence/submat.o       \
                        $(top_srcdir)/src/sequence/codonsubmat.o  \
                        $(top_srcdir)/src/struct/recyclebin.o     \
                        $(top_srcdir)/src/struct/pqueue.o         \
                        $(top_srcdir)/src/general/threadref.o     \
                        $(SEQUENCE_OBJ)

match_test_SOURCES = match.test.c match.c
match_test_LDADD = $(top_srcdir)/src/sequence/submat.o      \
                   $(top_srcdir)/src/sequence/codonsubmat.o \
                   $(SEQUENCE_OBJ)

seeder_test_SOURCES = seeder.test.c seeder.c wordhood.c hspset.c \
                      match.c comparison.c
seeder_test_LDADD = $(top_srcdir)/src/sequence/submat.o      \
                    $(top_srcdir)/src/sequence/codonsubmat.o \
                    $(top_srcdir)/src/struct/recyclebin.o    \
                    $(top_srcdir)/src/struct/fsm.o           \
                    $(top_srcdir)/src/struct/vfsm.o          \
                    $(top_srcdir)/src/struct/pqueue.o        \
                    $(top_srcdir)/src/general/threadref.o    \
                    $(SEQUENCE_OBJ)

# Files to clear away

MAINTAINERCLEANFILES = Makefile.in

