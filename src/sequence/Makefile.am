
TESTS = sequence.test submat.test translate.test splice.test \
        alphabet.test codonsubmat.test
noinst_PROGRAMS = $(TESTS)

INCLUDES = -I$(top_srcdir)/src/struct \
           -I$(top_srcdir)/src/general

noinst_HEADERS = sequence.h submat.h translate.h splice.h alphabet.h \
                 codonsubmat.h

splice_test_SOURCES = splice.test.c splice.c
SPLICE_OBJ = $(top_srcdir)/src/struct/matrix.o      \
             $(top_srcdir)/src/struct/sparsecache.o \
             $(top_srcdir)/src/general/argument.o   \
             $(top_srcdir)/src/general/lineparse.o  \
             -lm
splice_test_LDADD = $(SPLICE_OBJ)

sequence_test_SOURCES = sequence.test.c sequence.c alphabet.c translate.c \
                        splice.c
sequence_test_LDADD = $(SPLICE_OBJ)
submat_test_SOURCES = submat.test.c submat.c

translate_test_SOURCES = translate.test.c translate.c
translate_test_LDADD = $(top_srcdir)/src/general/argument.o

alphabet_test_SOURCES = alphabet.test.c alphabet.c
alphabet_test_LDADD = $(top_srcdir)/src/general/argument.o

codonsubmat_test_SOURCES = codonsubmat.test.c codonsubmat.c
codonsubmat_test_LDADD = $(top_srcdir)/src/general/argument.o -lm

# Files to clear away

MAINTAINERCLEANFILES = Makefile.in

