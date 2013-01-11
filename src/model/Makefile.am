
TESTS = edit_distance.test affine.test protein2dna.test    \
        ner.test est2genome.test ungapped.test intron.test \
        protein2genome.test coding2coding.test             \
        frameshift.test coding2genome.test                 \
        cdna2genome.test genome2genome.test                \
        phase.test modeltype.test

noinst_PROGRAMS = $(TESTS) bootstrapper

AM_CPPFLAGS = -I$(top_srcdir)/src/c4                    \
           -I$(top_srcdir)/src/bsdp                  \
           -I$(top_srcdir)/src/sdp                   \
           -I$(top_srcdir)/src/sequence              \
           -I$(top_srcdir)/src/general               \
           -I$(top_srcdir)/src/comparison            \
           -I$(top_srcdir)/src/struct                \
           -DSOURCE_ROOT_DIR="\"@source_root_dir@\"" \
           -DGLIB_CFLAGS="\"@glib_cflags@\""

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

ALIGNMENT_OBJ = $(top_srcdir)/src/comparison/match.o     \
                $(top_srcdir)/src/sequence/codonsubmat.o \
                $(top_srcdir)/src/sequence/submat.o      \
                $(SEQUENCE_OBJ)

C4_OBJECTS = $(top_srcdir)/src/struct/slist.o       \
             $(top_srcdir)/src/struct/rangetree.o   \
             $(top_srcdir)/src/struct/recyclebin.o  \
             $(top_srcdir)/src/c4/optimal.o         \
             $(top_srcdir)/src/c4/viterbi.o         \
             $(top_srcdir)/src/c4/layout.o          \
             $(top_srcdir)/src/c4/codegen.o         \
             $(top_srcdir)/src/c4/cgutil.o          \
             $(top_srcdir)/src/c4/alignment.o       \
             $(top_srcdir)/src/c4/region.o          \
             $(top_srcdir)/src/c4/c4.o              \
             $(top_srcdir)/src/c4/subopt.o          \
             $(ALIGNMENT_OBJ)

# -lm should be at end of argument list for IRIX

noinst_HEADERS = edit_distance.h affine.h protein2dna.h        \
                 ner.h est2genome.h ungapped.h intron.h        \
                 protein2genome.h coding2coding.h frameshift.h \
                 coding2genome.h genome2genome.h cdna2genome.h \
                 phase.h modeltype.h

bootstrapper_SOURCES = bootstrapper.c modeltype.c                     \
                       ungapped.c affine.c est2genome.c ner.c         \
                       protein2dna.c protein2genome.c coding2coding.c \
                       coding2genome.c genome2genome.c cdna2genome.c  \
                       frameshift.c intron.c phase.c
bootstrapper_LDADD = $(top_srcdir)/src/struct/pqueue.o          \
                     $(top_srcdir)/src/comparison/comparison.o  \
                     $(top_srcdir)/src/comparison/hspset.o      \
                     $(top_srcdir)/src/comparison/wordhood.o    \
                     $(top_srcdir)/src/bsdp/heuristic.o         \
                     $(top_srcdir)/src/sdp/sdp.o                \
                     $(top_srcdir)/src/sdp/scheduler.o          \
                     $(top_srcdir)/src/sdp/straceback.o         \
                     $(top_srcdir)/src/sdp/lookahead.o          \
                     $(top_srcdir)/src/sdp/boundary.o           \
                     $(C4_OBJECTS)

edit_distance_test_SOURCES = edit_distance.test.c edit_distance.c
edit_distance_test_LDADD = $(C4_OBJECTS)

affine_test_SOURCES = affine.test.c affine.c ungapped.c
affine_test_LDADD = $(C4_OBJECTS)

protein2dna_test_SOURCES = protein2dna.test.c protein2dna.c \
                           affine.c ungapped.c frameshift.c
protein2dna_test_LDADD = $(C4_OBJECTS)

ner_test_SOURCES = ner.test.c ner.c affine.c ungapped.c
ner_test_LDADD = $(C4_OBJECTS)

est2genome_test_SOURCES = est2genome.test.c est2genome.c \
                          affine.c ungapped.c intron.c
est2genome_test_LDADD = $(C4_OBJECTS)

ungapped_test_SOURCES = ungapped.test.c ungapped.c
ungapped_test_LDADD = $(C4_OBJECTS)

intron_test_SOURCES = intron.test.c intron.c
intron_test_LDADD = $(C4_OBJECTS)

protein2genome_test_SOURCES = protein2genome.test.c protein2genome.c \
                             protein2dna.c ungapped.c intron.c       \
                             affine.c frameshift.c phase.c
protein2genome_test_LDADD = $(C4_OBJECTS)

coding2coding_test_SOURCES = coding2coding.test.c coding2coding.c \
                           affine.c ungapped.c frameshift.c
coding2coding_test_LDADD = $(C4_OBJECTS)

frameshift_test_SOURCES = frameshift.test.c frameshift.c
frameshift_test_LDADD = $(C4_OBJECTS)

coding2genome_test_SOURCES = coding2genome.test.c coding2genome.c \
                           coding2coding.c affine.c ungapped.c    \
                           frameshift.c intron.c phase.c
coding2genome_test_LDADD = $(C4_OBJECTS)

cdna2genome_test_SOURCES = cdna2genome.test.c cdna2genome.c     \
                           affine.c ungapped.c                  \
                           frameshift.c intron.c phase.c        \
                           coding2coding.c coding2genome.c      \
                           est2genome.c
cdna2genome_test_LDADD = $(C4_OBJECTS)

genome2genome_test_SOURCES = genome2genome.test.c genome2genome.c \
                           coding2coding.c affine.c ungapped.c    \
                           frameshift.c intron.c phase.c          \
                           coding2genome.c cdna2genome.c          \
                           est2genome.c
genome2genome_test_LDADD = $(C4_OBJECTS)

phase_test_SOURCES = phase.test.c phase.c intron.c
phase_test_LDADD = $(C4_OBJECTS)

modeltype_test_SOURCES = modeltype.test.c modeltype.c           \
                         ungapped.c affine.c est2genome.c ner.c \
                         protein2dna.c protein2genome.c         \
                         coding2coding.c coding2genome.c        \
                         genome2genome.c cdna2genome.c          \
                         frameshift.c phase.c intron.c
modeltype_test_LDADD = $(C4_OBJECTS)

# Files to clear away

MAINTAINERCLEANFILES = Makefile.in

