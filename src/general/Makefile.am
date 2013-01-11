
TESTS = argument.test lineparse.test compoundfile.test socket.test \
        jobqueue.test threadref.test
noinst_PROGRAMS = $(TESTS)

noinst_HEADERS = argument.h lineparse.h compoundfile.h socket.h \
                 jobqueue.h threadref.h

AM_CPPFLAGS = -DSOURCE_ROOT_DIR="\"@source_root_dir@\"" \
           -DCUSTOM_GUINT64_FORMAT="\"@custom_guint64_format@\"" \
           -I$(top_srcdir)/src/struct

argument_test_SOURCES = argument.test.c argument.c
lineparse_test_SOURCES = lineparse.test.c lineparse.c
compoundfile_test_SOURCES = compoundfile.test.c compoundfile.c
socket_test_SOURCES = socket.test.c socket.c
jobqueue_test_SOURCES = jobqueue.test.c jobqueue.c
jobqueue_test_LDADD = $(top_srcdir)/src/struct/pqueue.o     \
                      $(top_srcdir)/src/struct/recyclebin.o
threadref_test_SOURCES = threadref.test.c threadref.c
threadref_test_LDADD = -lpthread

# Files to clear away

MAINTAINERCLEANFILES = Makefile.in

