/****************************************************************\
*                                                                *
*  Library for manipulation of FASTA format databases            *
*                                                                *
*  Guy St.C. Slater..   mailto:guy@ebi.ac.uk                     *
*  Copyright (C) 2000-2009.  All Rights Reserved.                *
*                                                                *
*  This source code is distributed under the terms of the        *
*  GNU General Public License, version 3. See the file COPYING   *
*  or http://www.gnu.org/licenses/gpl.txt for details            *
*                                                                *
*  If you use this code, please keep this notice intact.         *
*                                                                *
\****************************************************************/

#ifndef INCLUDED_FASTADB_H
#define INCLUDED_FASTADB_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <stdio.h>
#include <glib.h>

#include "compoundfile.h"
#include "sequence.h"
#include "argument.h"
#include "sparsecache.h"

typedef struct {
    gchar *suffix_filter;
} FastaDB_ArgumentSet;

FastaDB_ArgumentSet *FastaDB_ArgumentSet_create(Argument *arg);

typedef enum {
    FastaDB_Mask_ID  = (1<<1),
    FastaDB_Mask_DEF = (1<<2),
    FastaDB_Mask_SEQ = (1<<3),
    FastaDB_Mask_LEN = (1<<4),
    FastaDB_Mask_ALL = (~0)
} FastaDB_Mask;

typedef struct FastaDB {
            guint  ref_count;
         Alphabet *alphabet;
     CompoundFile *cf;
            gchar *out_buffer;
            guint  out_buffer_pos;
            guint  out_buffer_alloc;
             gint  line_length;
} FastaDB;
/* line_length is used for fasta file random-access
 * it is set to zero for irregular line lengths
 * is should not be used until the entire file has been parsed.
 */

typedef struct {
                    guint  ref_count;
                  FastaDB *source;
    CompoundFile_Location *location;
                 Sequence *seq;
} FastaDB_Seq;

typedef gboolean (*FastaDB_TraverseFunc)(FastaDB_Seq *fdbs,
                                         gpointer user_data);
/* Return TRUE to stop the traversal */

    FastaDB *FastaDB_open_list(GPtrArray *path_list,
                               Alphabet *alphabet);
    FastaDB *FastaDB_open_list_with_limit(GPtrArray *path_list,
             Alphabet *alphabet, gint chunk_id, gint chunk_total);
    FastaDB *FastaDB_open(gchar *path, Alphabet *alphabet);
    FastaDB *FastaDB_share(FastaDB *fdb);
    FastaDB *FastaDB_dup(FastaDB *fdb); /* For use in a separate thread */
       void  FastaDB_close(FastaDB *fdb);
       void  FastaDB_rewind(FastaDB *fdb);
   gboolean  FastaDB_is_finished(FastaDB *fdb);
       void  FastaDB_traverse(FastaDB *fdb, FastaDB_Mask mask,
                     FastaDB_TraverseFunc fdtf, gpointer user_data);
      gsize  FastaDB_memory_usage(FastaDB *fdb);

     FastaDB_Seq *FastaDB_next(FastaDB *fdb, FastaDB_Mask mask);
CompoundFile_Pos  FastaDB_find_next_start(FastaDB *fdb,
                                          CompoundFile_Pos pos);

    gboolean FastaDB_file_is_fasta(gchar *path);
    /* Returns true if first non-whitespace character in file is '>' */

typedef struct {
                FastaDB *source;
  CompoundFile_Location *location;
        Sequence_Strand  strand;
                   gint  seq_offset; /* for random access */
                   gint  length;     /* for random access */
} FastaDB_Key;

FastaDB_Seq *FastaDB_fetch(FastaDB *fdb, FastaDB_Mask mask,
                           CompoundFile_Pos pos);

FastaDB_Key *FastaDB_Key_create(FastaDB *source,
                                CompoundFile_Location *location,
                                Sequence_Strand strand,
                                gint seq_offset, gint length);
FastaDB_Key *FastaDB_Seq_get_key(FastaDB_Seq *fdbs);
FastaDB_Seq *FastaDB_Key_get_seq(FastaDB_Key *fdbk, FastaDB_Mask mask);
       void  FastaDB_Key_destroy(FastaDB_Key *fdbk);
      gchar *FastaDB_Key_get_def(FastaDB_Key *fdbk);
SparseCache *FastaDB_Key_get_SparseCache(FastaDB_Key *fdbk);
       void  FastaDB_SparseCache_compress(SparseCache_Page *page, gint len);


FastaDB_Seq **FastaDB_all(gchar *path, Alphabet *alphabet,
                          FastaDB_Mask mask, guint *total);

FastaDB_Seq *FastaDB_Seq_share(FastaDB_Seq *fdbs);
       void  FastaDB_Seq_destroy(FastaDB_Seq *fdbs);
FastaDB_Seq *FastaDB_Seq_revcomp(FastaDB_Seq *fdbs);
       void  FastaDB_Seq_all_destroy(FastaDB_Seq **fdbs);

       gint  FastaDB_Seq_print(FastaDB_Seq *fdbs, FILE *fp,
                               FastaDB_Mask mask);
       gint  FastaDB_Seq_all_print(FastaDB_Seq **fdbs, FILE *fp,
                                   FastaDB_Mask mask);

FastaDB_Seq *FastaDB_get_single(gchar *path, Alphabet *alphabet);
Alphabet_Type FastaDB_guess_type(gchar *path);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* INCLUDED_FASTADB_H */

