
bin_PROGRAMS = exonerate ipcress exonerate-server

AM_CPPFLAGS = -I$(top_srcdir)/src/sequence   \
           -I$(top_srcdir)/src/database   \
           -I$(top_srcdir)/src/comparison \
           -I$(top_srcdir)/src/struct     \
           -I$(top_srcdir)/src/c4         \
           -I$(top_srcdir)/src/bsdp       \
           -I$(top_srcdir)/src/sdp        \
           -I$(top_srcdir)/src/model      \
           -I$(top_srcdir)/src/hub        \
           -I$(top_srcdir)/src/general    \
           -DHOSTTYPE="\"@host@\""        \
           -DSOURCE_ROOT_DIR="\"@source_root_dir@\"" \
           -DCUSTOM_GUINT64_FORMAT="\"@custom_guint64_format@\""

# Missing viterb.o and scheduler.o so codegen versions may be substituted
C4_OBJECTS = $(top_srcdir)/src/c4/c4.o                \
             $(top_srcdir)/src/c4/codegen.o           \
             $(top_srcdir)/src/c4/cgutil.o            \
             $(top_srcdir)/src/c4/opair.o             \
             $(top_srcdir)/src/c4/alignment.o         \
             $(top_srcdir)/src/c4/optimal.o           \
             $(top_srcdir)/src/c4/layout.o            \
             $(top_srcdir)/src/c4/region.o            \
             $(top_srcdir)/src/c4/subopt.o            \
             $(top_srcdir)/src/bsdp/bsdp.o            \
             $(top_srcdir)/src/bsdp/sar.o             \
             $(top_srcdir)/src/bsdp/heuristic.o       \
             $(top_srcdir)/src/bsdp/hpair.o           \
             $(top_srcdir)/src/sdp/boundary.o         \
             $(top_srcdir)/src/sdp/sdp.o              \
             $(top_srcdir)/src/sdp/lookahead.o        \
             $(top_srcdir)/src/sdp/straceback.o       \
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
             $(top_srcdir)/src/general/socket.o       \
             $(top_srcdir)/src/general/jobqueue.o     \
             $(top_srcdir)/src/comparison/match.o     \
             $(top_srcdir)/src/sequence/submat.o      \
             $(top_srcdir)/src/sequence/codonsubmat.o

EXTRA_exonerate_SOURCES = c4_model_archive.a c4_model_archive.h

BUILT_SOURCES = c4_model_archive.a
BOOTSTRAPPER = $(top_srcdir)/src/model/bootstrapper
c4_model_archive.h: $(BOOTSTRAPPER)
	$(BOOTSTRAPPER) --compiled no
c4_model_archive.a: $(BOOTSTRAPPER)
	$(BOOTSTRAPPER) --compiled no

# Special viterbi.o and scheduler.o for compiled model linking
viterbi.o: $(top_srcdir)/src/c4/viterbi.c
	$(CC) $(CFLAGS) $(AM_CPPFLAGS) -I.   \
        -DUSE_COMPILED_MODELS             \
        -c $(top_srcdir)/src/c4/viterbi.c
scheduler.o: $(top_srcdir)/src/sdp/scheduler.c
	$(CC) $(CFLAGS) $(AM_CPPFLAGS) -I.   \
        -DUSE_COMPILED_MODELS             \
        -c $(top_srcdir)/src/sdp/scheduler.c

exonerate_COMPONENTS = @codegen_extra_ldadd@                     \
                       $(top_srcdir)/src/comparison/wordhood.o   \
                       $(top_srcdir)/src/comparison/hspset.o     \
                       $(top_srcdir)/src/comparison/seeder.o     \
                       $(top_srcdir)/src/comparison/comparison.o \
                       $(top_srcdir)/src/database/fastapipe.o    \
                       $(top_srcdir)/src/database/fastadb.o      \
                       $(top_srcdir)/src/struct/dejavu.o         \
                       $(top_srcdir)/src/struct/fsm.o            \
                       $(top_srcdir)/src/struct/vfsm.o           \
                       $(top_srcdir)/src/struct/sparsecache.o    \
                       $(top_srcdir)/src/hub/gam.o               \
                       $(top_srcdir)/src/hub/bsam.o              \
                       $(top_srcdir)/src/hub/analysis.o          \
                       $(top_srcdir)/src/general/compoundfile.o  \
                       $(top_srcdir)/src/general/threadref.o     \
                       $(C4_OBJECTS)

exonerate_SOURCES = exonerate.c
exonerate_DEPENDENCIES = @codegen_extra_sources@ Makefile \
                         $(exonerate_COMPONENTS)
exonerate_LDADD = $(exonerate_COMPONENTS) -lm

# Remove built sources from distribution
# (will be able to use nodist_exonerate_SOURCES in automake 1.5)
# nodist_exonerate_SOURCES = c4_model_archive.a
dist-hook:
	for file in $(EXTRA_exonerate_SOURCES) ; \
            do rm -f $(distdir)/$$file ; done ;

ipcress_SOURCES = ipcress.c
ipcress_LDADD = $(top_srcdir)/src/comparison/pcr.o       \
                $(top_srcdir)/src/comparison/wordhood.o  \
                $(top_srcdir)/src/database/fastadb.o     \
                $(top_srcdir)/src/sequence/alphabet.o    \
                $(top_srcdir)/src/sequence/sequence.o    \
                $(top_srcdir)/src/sequence/splice.o      \
                $(top_srcdir)/src/sequence/translate.o   \
                $(top_srcdir)/src/sequence/submat.o      \
                $(top_srcdir)/src/sequence/codonsubmat.o \
                $(top_srcdir)/src/struct/fsm.o           \
                $(top_srcdir)/src/struct/slist.o         \
                $(top_srcdir)/src/struct/matrix.o        \
                $(top_srcdir)/src/struct/recyclebin.o    \
                $(top_srcdir)/src/struct/sparsecache.o   \
                $(top_srcdir)/src/general/lineparse.o    \
                $(top_srcdir)/src/general/argument.o     \
                $(top_srcdir)/src/general/compoundfile.o \
                -lm
# Files to clear away

exonerate_server_SOURCES = exonerate-server.c
exonerate_server_LDADD   = $(top_srcdir)/src/general/socket.o       \
                           $(top_srcdir)/src/general/argument.o     \
                           $(top_srcdir)/src/comparison/hspset.o    \
                           $(top_srcdir)/src/comparison/match.o     \
                           $(top_srcdir)/src/comparison/wordhood.o  \
                           $(top_srcdir)/src/database/dataset.o     \
                           $(top_srcdir)/src/database/index.o       \
                           $(top_srcdir)/src/database/fastadb.o     \
                           $(top_srcdir)/src/struct/bitarray.o      \
                           $(top_srcdir)/src/struct/sparsecache.o   \
                           $(top_srcdir)/src/struct/matrix.o        \
                           $(top_srcdir)/src/struct/splaytree.o     \
                           $(top_srcdir)/src/struct/noitree.o       \
                           $(top_srcdir)/src/struct/vfsm.o          \
                           $(top_srcdir)/src/struct/pqueue.o        \
                           $(top_srcdir)/src/struct/rangetree.o     \
                           $(top_srcdir)/src/struct/recyclebin.o    \
                           $(top_srcdir)/src/sequence/sequence.o    \
                           $(top_srcdir)/src/sequence/alphabet.o    \
                           $(top_srcdir)/src/sequence/translate.o   \
                           $(top_srcdir)/src/sequence/splice.o      \
                           $(top_srcdir)/src/sequence/submat.o      \
                           $(top_srcdir)/src/sequence/codonsubmat.o \
                           $(top_srcdir)/src/general/compoundfile.o \
                           $(top_srcdir)/src/general/lineparse.o    \
                           $(top_srcdir)/src/general/threadref.o    \
                           -lm

CLEANFILES = $(EXTRA_exonerate_SOURCES) @codegen_extra_sources@ \
             @codegen_extra_ld_add@
MAINTAINERCLEANFILES = Makefile.in

