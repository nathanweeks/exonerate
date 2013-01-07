/****************************************************************\
*                                                                *
*  Library for manipulation of exonerate dataset files           *
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

#ifndef INCLUDED_DATASET_H
#define INCLUDED_DATASET_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <stdio.h>
#include <glib.h>

#ifdef USE_PTHREADS
#include <pthread.h>
#endif /* USE_PTHREADS */

#include "alphabet.h"
#include "sequence.h"
#include "fastadb.h"

/* File format:
   Header
   path_data: for each file
       path def\n
   seq_data: for each seq (in alphabetical order of id)
       id\n
   seq_info: for each seq
       database <ND>
       offset <MD>
       length <MS>
       gcg_checksum <14>
*/

typedef struct {
    guint64 magic;
    guint64 version;
    guint64 type;
    guint64 line_length;
    /**/
    guint64 number_of_dbs;
    guint64 max_db_len;
    guint64 total_db_len;
    /**/
    guint64 number_of_seqs;
    guint64 max_seq_len;
    guint64 total_seq_len;
    /**/
    guint64 path_data_offset;
    guint64 seq_data_offset;
    guint64 seq_info_offset;
    guint64 total_file_length;
} Dataset_Header;

typedef struct {
     gint num_db_width;
     gint max_db_len_width;
     gint max_seq_len_width;
    gsize seq_data_item_size;
} Dataset_Width;

typedef struct {
     FastaDB_Key  *key;
         guint64   gcg_checksum;
           gchar  *id;
           gchar  *def;
            gint   pos;
        Sequence  *cache_seq;
} Dataset_Sequence;

typedef struct {
             guint  ref_count;
          Alphabet *alphabet;
    Dataset_Header *header;
     Dataset_Width *width;
         GPtrArray *seq_list;  /* containing Dataset_Sequence objects */
           FastaDB *fdb;
#ifdef USE_PTHREADS
    pthread_mutex_t dataset_mutex;
#endif /* USE_PTHREADS */
} Dataset;

Dataset *Dataset_create(GPtrArray *path_list,
                        Alphabet_Type alphabet_type,
                        gboolean softmask_input);
 Dataset *Dataset_share(Dataset *dataset);
    void  Dataset_destroy(Dataset *dataset);
    void  Dataset_info(Dataset *dataset);
   gsize  Dataset_memory_usage(Dataset *dataset);
gboolean  Dataset_check_filetype(gchar *path);
/* Returns TRUE when magic number is correct for this filetype */
    void  Dataset_preload_sequences(Dataset *dataset);

   void  Dataset_write(Dataset *dataset, gchar *path);
Dataset *Dataset_read(gchar *path);
   void  Dataset_preload_seqs(Dataset *dataset);

    gint  Dataset_lookup_id(Dataset *dataset, gchar *id);
Sequence *Dataset_get_sequence(Dataset *dataset, gint dataset_pos);
/* Sequence returned will be Sequence_Type_EXTMEM */

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* INCLUDED_DATASET_H */

