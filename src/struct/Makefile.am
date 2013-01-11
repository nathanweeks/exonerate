
TESTS = fsm.test matrix.test pqueue.test slist.test rangetree.test \
        vfsm.test recyclebin.test sparsecache.test dejavu.test     \
        bitarray.test splaytree.test noitree.test 
#       huffman.test fmindex.test

noinst_PROGRAMS = $(TESTS)

noinst_HEADERS = fsm.h matrix.h pqueue.h slist.h rangetree.h vfsm.h \
                 recyclebin.h sparsecache.h dejavu.h bitarray.h     \
                 splaytree.h noitree.h
#                huffman.h fmindex.h

AM_CPPFLAGS = -DCUSTOM_GUINT64_FORMAT="\"@custom_guint64_format@\""

fsm_test_SOURCES = fsm.test.c fsm.c
fsm_test_LDADD = recyclebin.o

matrix_test_SOURCES = matrix.test.c matrix.c
pqueue_test_SOURCES = pqueue.test.c pqueue.c recyclebin.c
slist_test_SOURCES = slist.test.c slist.c recyclebin.c
rangetree_test_SOURCES = rangetree.test.c rangetree.c
rangetree_test_LDADD = recyclebin.o

vfsm_test_SOURCES = vfsm.test.c vfsm.c
vfsm_test_LDADD = -lm

recyclebin_test_SOURCES = recyclebin.test.c recyclebin.c
sparsecache_test_SOURCES = sparsecache.test.c sparsecache.c
dejavu_test_SOURCES = dejavu.test.c dejavu.c
bitarray_test_SOURCES = bitarray.test.c bitarray.c

splaytree_test_SOURCES = splaytree.test.c splaytree.c
splaytree_test_LDADD = recyclebin.o

noitree_test_SOURCES = noitree.test.c noitree.c
noitree_test_LDADD = recyclebin.o splaytree.o

# huffman_test_SOURCES = huffman.test.c huffman.c
# huffman_test_LDADD = bitarray.o pqueue.o recyclebin.o

# fmindex_test_SOURCES = fmindex.test.c fmindex.c
# fmindex_test_LDADD = bitarray.o pqueue.o recyclebin.o huffman.o

# Files to clear away

MAINTAINERCLEANFILES = Makefile.in

