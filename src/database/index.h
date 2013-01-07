/****************************************************************\
*                                                                *
*  Library for manipulation of exonerate index files             *
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

#ifndef INCLUDED_INDEX_H
#define INCLUDED_INDEX_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <stdio.h>
#include <glib.h>
#include <sys/types.h>
#include <unistd.h>

#ifdef USE_PTHREADS
#include <pthread.h>
#endif /* USE_PTHREADS */

#include "dataset.h"
#include "vfsm.h"
#include "hspset.h"
#include "bitarray.h"

/* File format:
   Header
   dataset path\n
   FW Strand
   RV Strand (if translated)

   Strand:
       Strand header
       WordList:
           For each word
               word_id <MW>
               freq_count <MI>
               index_offset <TI>
       Index:
           sequence <NS>
           pos <MS> (from dataset)
*/

typedef struct {
    guint64  magic;
    guint64  version;
    guint64  type;                /* plain | trans */
    guint64  dataset_path_len;
    /**/
    guint64  word_length;
    guint64  word_jump;
    guint64  word_ambiguity;
    guint64  saturate_threshold;
} Index_Header;

typedef struct {
    gint  max_word_width;       /* From vfsm->lrw : MW */
    gint  number_of_seqs_width; /* From dataset->header->number_of_seqs : NS */
} Index_Width;

typedef struct {
    gint sequence_id;
    gint position;
} Index_Address;

typedef struct {
      gint freq_count;
    gint64 index_offset;
} Index_Word;

typedef struct {
    guint64 max_index_length;    /* Filled by Index_survey_word_list() */
    guint64 word_list_length;    /* Filled by Index_survey_word_list() */
    guint64 total_index_length;  /* Filled by Index_find_offsets()     */
} Index_Strand_Header;

typedef struct {
    gint  max_index_len_width;    /* From   max_index_length : MI */
    gint  total_index_len_width;  /* From total_index_length : TI */
} Index_Strand_Width;

typedef struct {
     Index_Strand_Header  header;
      Index_Strand_Width  width;
                    /**/
                    gint *word_table; /* VFSM array */
              Index_Word *word_list;
                   off_t  strand_offset; /* Offset to strand header */
                BitArray *index_cache;
} Index_Strand;

typedef struct {
              guint  ref_count;
               FILE *fp;
              gchar *dataset_path;
            Dataset *dataset;
       Index_Header *header;
               VFSM *vfsm;
        Index_Width *width;
              /**/
       Index_Strand *forward;
       Index_Strand *revcomp; /* Only used when index is translated */
#ifdef USE_PTHREADS
     pthread_mutex_t index_mutex;
#endif /* USE_PTHREADS */
} Index;

/**/

   Index *Index_create(Dataset *dataset, gboolean is_translated, gint word_length,
                       gint word_jump, gint word_ambiguity, gint saturate_threshold,
                       gchar *index_path, gchar *dataset_path, gint memory_limit);
   Index *Index_share(Index *index);
    void  Index_destroy(Index *index);
    void  Index_info(Index *index);
   Index *Index_open(gchar *path);
 guint64  Index_memory_usage(Index *index);
    void  Index_preload_index(Index *index);
gboolean  Index_check_filetype(gchar *path);
/* Returns TRUE when magic number is correct for this filetype */

typedef struct {
    HSPset *hsp_set;
      gint  target_id;
} Index_HSPset;

void Index_HSPset_destroy(Index_HSPset *index_hsp_set);

GPtrArray *Index_get_HSPsets(Index *index, HSP_Param *hsp_param,
                             Sequence *query, gboolean revcomp_target);
/* Returns a GPtrArray containing Index_HSPset structs */

GPtrArray *Index_get_HSPsets_geneseed(Index *index, HSP_Param *hsp_param,
                                  Sequence *query, gboolean revcomp_target,
                                  gint geneseed_threshold, gint geneseed_repeat,
                                  gint max_query_span, gint max_target_span);

/**/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* INCLUDED_INDEX_H */

