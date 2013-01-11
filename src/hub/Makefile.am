
TESTS = gam.test analysis.test bsam.test
noinst_PROGRAMS = $(TESTS)

AM_CPPFLAGS = -I$(top_srcdir)/src/sequence   \
           -I$(top_srcdir)/src/database   \
           -I$(top_srcdir)/src/comparison \
           -I$(top_srcdir)/src/struct     \
           -I$(top_srcdir)/src/c4         \
           -I$(top_srcdir)/src/bsdp       \
           -I$(top_srcdir)/src/sdp        \
           -I$(top_srcdir)/src/model      \
           -I$(top_srcdir)/src/general    \
           -DHOSTTYPE="\"@host@\""        \
           -DSOURCE_ROOT_DIR="\"@source_root_dir@\""

noinst_HEADERS = gam.h analysis.h bsam.h

C4_OBJECTS = $(top_srcdir)/src/c4/c4.o                \
             $(top_srcdir)/src/c4/codegen.o           \
             $(top_srcdir)/src/c4/cgutil.o            \
             $(top_srcdir)/src/c4/opair.o             \
             $(top_srcdir)/src/c4/alignment.o         \
             $(top_srcdir)/src/c4/optimal.o           \
             $(top_srcdir)/src/c4/viterbi.o           \
             $(top_srcdir)/src/c4/layout.o            \
             $(top_srcdir)/src/c4/region.o            \
             $(top_srcdir)/src/c4/subopt.o            \
             $(top_srcdir)/src/bsdp/bsdp.o            \
             $(top_srcdir)/src/bsdp/sar.o             \
             $(top_srcdir)/src/bsdp/heuristic.o       \
             $(top_srcdir)/src/bsdp/hpair.o           \
             $(top_srcdir)/src/sdp/boundary.o         \
             $(top_srcdir)/src/sdp/sdp.o              \
             $(top_srcdir)/src/sdp/scheduler.o        \
             $(top_srcdir)/src/sdp/straceback.o       \
             $(top_srcdir)/src/sdp/lookahead.o        \
             $(top_srcdir)/src/struct/matrix.o        \
             $(top_srcdir)/src/struct/slist.o         \
             $(top_srcdir)/src/struct/pqueue.o        \
             $(top_srcdir)/src/struct/recyclebin.o    \
             $(top_srcdir)/src/struct/rangetree.o     \
             $(top_srcdir)/src/sequence/sequence.o    \
             $(top_srcdir)/src/sequence/alphabet.o    \
             $(top_srcdir)/src/sequence/translate.o   \
             $(top_srcdir)/src/sequence/splice.o      \
             $(top_srcdir)/src/model/affine.o         \
             $(top_srcdir)/src/model/est2genome.o     \
             $(top_srcdir)/src/model/ner.o            \
             $(top_srcdir)/src/model/protein2dna.o    \
             $(top_srcdir)/src/model/protein2genome.o \
             $(top_srcdir)/src/model/coding2coding.o  \
             $(top_srcdir)/src/model/coding2genome.o  \
             $(top_srcdir)/src/model/cdna2genome.o    \
             $(top_srcdir)/src/model/genome2genome.o  \
             $(top_srcdir)/src/model/ungapped.o       \
             $(top_srcdir)/src/model/intron.o         \
             $(top_srcdir)/src/model/frameshift.o     \
             $(top_srcdir)/src/model/phase.o          \
             $(top_srcdir)/src/model/modeltype.o      \
             $(top_srcdir)/src/general/argument.o     \
             $(top_srcdir)/src/general/lineparse.o    \
             $(top_srcdir)/src/general/jobqueue.o     \
             $(top_srcdir)/src/comparison/match.o     \
             $(top_srcdir)/src/sequence/submat.o      \
             $(top_srcdir)/src/sequence/codonsubmat.o \
             -lm

gam_test_SOURCES = gam.test.c gam.c
gam_test_LDADD = $(top_srcdir)/src/comparison/hspset.o     \
                 $(top_srcdir)/src/comparison/comparison.o \
                 $(top_srcdir)/src/comparison/wordhood.o   \
                 $(top_srcdir)/src/struct/sparsecache.o    \
                 $(top_srcdir)/src/general/compoundfile.o  \
                 $(top_srcdir)/src/general/threadref.o     \
                 $(C4_OBJECTS)

bsam_test_SOURCES = bsam.test.c bsam.c
bsam_test_LDADD = $(top_srcdir)/src/general/argument.o       \
                  $(top_srcdir)/src/general/lineparse.o      \
                  $(top_srcdir)/src/general/threadref.o      \
                  $(top_srcdir)/src/comparison/hspset.o      \
                  $(top_srcdir)/src/comparison/comparison.o  \
                  $(top_srcdir)/src/comparison/match.o       \
                  $(top_srcdir)/src/comparison/wordhood.o    \
                  $(top_srcdir)/src/sequence/sequence.o      \
                  $(top_srcdir)/src/sequence/alphabet.o      \
                  $(top_srcdir)/src/sequence/translate.o     \
                  $(top_srcdir)/src/sequence/submat.o        \
                  $(top_srcdir)/src/sequence/codonsubmat.o   \
                  $(top_srcdir)/src/sequence/splice.o        \
                  $(top_srcdir)/src/struct/dejavu.o          \
                  $(top_srcdir)/src/struct/pqueue.o          \
                  $(top_srcdir)/src/struct/sparsecache.o     \
                  $(top_srcdir)/src/struct/recyclebin.o      \
                  $(top_srcdir)/src/struct/matrix.o          \
                  -lm
analysis_test_SOURCES = analysis.test.c analysis.c gam.c bsam.c
analysis_test_LDADD = $(top_srcdir)/src/comparison/hspset.o     \
                      $(top_srcdir)/src/comparison/wordhood.o   \
                      $(top_srcdir)/src/comparison/seeder.o     \
                      $(top_srcdir)/src/comparison/comparison.o \
                      $(top_srcdir)/src/database/fastapipe.o    \
                      $(top_srcdir)/src/database/fastadb.o      \
                      $(top_srcdir)/src/struct/dejavu.o         \
                      $(top_srcdir)/src/struct/fsm.o            \
                      $(top_srcdir)/src/struct/vfsm.o           \
                      $(top_srcdir)/src/struct/sparsecache.o    \
                      $(top_srcdir)/src/general/compoundfile.o  \
                      $(top_srcdir)/src/general/socket.o        \
                      $(top_srcdir)/src/general/threadref.o      \
                      $(C4_OBJECTS)

# Files to clear away

MAINTAINERCLEANFILES = Makefile.in

