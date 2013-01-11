
bin_PROGRAMS = @installed_util_list@

EXTRA_PROGRAMS = fastaclean fastaclip fastachecksum fastacomposition \
               fastadiff fastaexplode fastafetch fastahardmask \
               fastaindex fastalength fastanrdb fastaoverlap \
               fastareformat fastaremove fastarevcomp fastasoftmask \
               fastasort fastasplit fastasubseq fastatranslate \
               fastavalidcds fasta2esd fastaannotatecdna \
               esd2esi

AM_CPPFLAGS = -I$(top_srcdir)/src/sequence   \
           -I$(top_srcdir)/src/database   \
           -I$(top_srcdir)/src/struct     \
           -I$(top_srcdir)/src/comparison \
           -I$(top_srcdir)/src/general

LDADD = $(top_srcdir)/src/sequence/sequence.o    \
        $(top_srcdir)/src/sequence/alphabet.o    \
        $(top_srcdir)/src/sequence/splice.o      \
        $(top_srcdir)/src/sequence/translate.o   \
        $(top_srcdir)/src/database/fastadb.o     \
        $(top_srcdir)/src/general/argument.o     \
        $(top_srcdir)/src/general/compoundfile.o \
        $(top_srcdir)/src/general/lineparse.o    \
        $(top_srcdir)/src/struct/sparsecache.o   \
        $(top_srcdir)/src/struct/matrix.o        \
        -lm

fastaclean_SOURCES = fastaclean.c
fastaclip_SOURCES = fastaclip.c
fastachecksum_SOURCES = fastachecksum.c
fastacomposition_SOURCES = fastacomposition.c
fastadiff_SOURCES = fastadiff.c
fastaexplode_SOURCES = fastaexplode.c
fastafetch_SOURCES = fastafetch.c
fastahardmask_SOURCES = fastahardmask.c
fastaindex_SOURCES = fastaindex.c
fastalength_SOURCES = fastalength.c
fastanrdb_SOURCES = fastanrdb.c
fastaoverlap_SOURCES = fastaoverlap.c
fastareformat_SOURCES = fastareformat.c
fastaremove_SOURCES = fastaremove.c
fastarevcomp_SOURCES = fastarevcomp.c
fastasoftmask_SOURCES = fastasoftmask.c
fastasort_SOURCES = fastasort.c
fastasplit_SOURCES = fastasplit.c
fastasubseq_SOURCES = fastasubseq.c
fastatranslate_SOURCES = fastatranslate.c
fasta2esd_SOURCES = fasta2esd.c
fastaannotatecdna_SOURCES = fastaannotatecdna.c
esd2esi_SOURCES = esd2esi.c

fasta2esd_LDADD  = $(top_srcdir)/src/database/dataset.o \
                   $(top_srcdir)/src/struct/bitarray.o  \
                   $(LDADD)

esd2esi_LDADD  = $(top_srcdir)/src/database/dataset.o     \
                 $(top_srcdir)/src/database/index.o       \
                 $(top_srcdir)/src/general/threadref.o    \
                 $(top_srcdir)/src/struct/bitarray.o      \
                 $(top_srcdir)/src/struct/vfsm.o          \
                 $(top_srcdir)/src/struct/pqueue.o        \
                 $(top_srcdir)/src/struct/recyclebin.o    \
                 $(top_srcdir)/src/struct/rangetree.o     \
                 $(top_srcdir)/src/struct/noitree.o       \
                 $(top_srcdir)/src/struct/splaytree.o     \
                 $(top_srcdir)/src/sequence/submat.o      \
                 $(top_srcdir)/src/sequence/codonsubmat.o \
                 $(top_srcdir)/src/comparison/wordhood.o  \
                 $(top_srcdir)/src/comparison/hspset.o    \
                 $(top_srcdir)/src/comparison/match.o     \
                 $(LDADD)

EXTRA_fastavalidcds_SOURCES = fastavalidcds.c


# Files to clear away

MAINTAINERCLEANFILES = Makefile.in

