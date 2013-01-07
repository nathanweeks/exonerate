/****************************************************************\
*                                                                *
*  Library for HSP sets (high-scoring segment pairs)             *
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

/* Version 0.4
 */

#ifndef INCLUDED_HSPSET_H
#define INCLUDED_HSPSET_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <glib.h>

#include "match.h"
#include "argument.h"
#include "recyclebin.h"
#include "pqueue.h"
#include "matrix.h"
#include "wordhood.h"
#include "recyclebin.h"
#include "threadref.h"

typedef struct {
        gint filter_threshold;
    gboolean use_wordhood_dropoff;
        gint seed_repeat;
    /**/
        gint dna_wordlen;
        gint protein_wordlen;
        gint codon_wordlen;
    /**/
        gint dna_hsp_dropoff;
        gint protein_hsp_dropoff;
        gint codon_hsp_dropoff;
    /**/
        gint dna_hsp_threshold;
        gint protein_hsp_threshold;
        gint codon_hsp_threshold;
    /**/
        gint dna_word_limit;
        gint protein_word_limit;
        gint codon_word_limit;
    /**/
        gint geneseed_threshold;
        gint geneseed_repeat;
    /**/
} HSPset_ArgumentSet;

HSPset_ArgumentSet *HSPset_ArgumentSet_create(Argument *arg);

/**/

typedef struct {
        guint  query_start;
        guint  target_start;
        guint  length;   /* Length is number of match state visits */
  Match_Score  score;
        guint  cobs;     /* cobs == Centre Offset By Score         */
struct HSPset *hsp_set;  /* Never included in hsp_set->ref_count   */
} HSP;
/* FIXME: remove hsp_set from HSP to save space ? */

void HSP_destroy(HSP *hsp);

#define HSP_query_advance(hsp) \
    ((hsp)->hsp_set->param->match->query->advance)

#define HSP_target_advance(hsp) \
    ((hsp)->hsp_set->param->match->target->advance)

#define HSP_query_end(hsp) \
    ((hsp)->query_start    \
  + ((hsp)->length*HSP_query_advance(hsp)))

#define HSP_target_end(hsp) \
    ((hsp)->target_start    \
  + ((hsp)->length*HSP_target_advance(hsp)))

#define HSP_query_cobs(hsp) \
    ((hsp)->query_start     \
   +((hsp)->cobs*HSP_query_advance(hsp)))

#define HSP_target_cobs(hsp) \
    ((hsp)->target_start     \
  + ((hsp)->cobs*HSP_target_advance(hsp)))

#define HSP_diagonal(hsp)                           \
    (((hsp)->target_start*HSP_query_advance(hsp))   \
    -((hsp)->query_start*HSP_target_advance(hsp)))
/* advance_{query,target} are swapped for position on diagonal */

#define HSP_get_score(hsp, query_pos, target_pos)    \
     ((hsp)->hsp_set->param->match->score_func(      \
      (hsp)->hsp_set->param->match,                  \
      (hsp)->hsp_set->query, (hsp)->hsp_set->target, \
      (query_pos), (target_pos)))

#define HSP_get_display(hsp, query_pos, target_pos, display_str)  \
     ((hsp)->hsp_set->param->match->display_func(                 \
      (hsp)->hsp_set->param->match,                               \
      (hsp)->hsp_set->query, (hsp)->hsp_set->target,              \
      (query_pos), (target_pos), display_str))

#define HSP_query_masked(hsp, query_pos)             \
    ((hsp)->hsp_set->param->match->query->mask_func( \
     (hsp)->hsp_set->param->match->query,            \
     (hsp)->hsp_set->query, (query_pos)))

#define HSP_target_masked(hsp, target_pos)            \
    ((hsp)->hsp_set->param->match->target->mask_func( \
     (hsp)->hsp_set->param->match->target,            \
     (hsp)->hsp_set->target, (target_pos)))

#define HSP_query_self(hsp, query_pos)               \
    ((hsp)->hsp_set->param->match->query->self_func( \
     (hsp)->hsp_set->param->match->query,            \
     (hsp)->hsp_set->query, (query_pos)))

#define HSP_target_self(hsp, target_pos)               \
    ((hsp)->hsp_set->param->match->target->self_func(  \
     (hsp)->hsp_set->param->match->target,             \
     (hsp)->hsp_set->target, (target_pos)))

#define HSPset_is_empty(hspset) ((hspset)->is_empty)

typedef struct HSP_Param {
         ThreadRef  *thread_ref;
HSPset_ArgumentSet  *has;
             Match  *match;
              gint   wordlen;
              gint   seedlen;
       Match_Score   dropoff;
       Match_Score   threshold;
       Match_Score   wordlimit;
          WordHood  *wordhood;
          gboolean   use_horizon;
              gint   seed_repeat;
        RecycleBin  *hsp_recycle;
#ifdef USE_PTHREADS
   pthread_mutex_t   hsp_recycle_lock;
#endif /* USE_PTHREADS */
} HSP_Param;

HSP_Param *HSP_Param_create(Match *match, gboolean use_horizon);
     void  HSP_Param_destroy(HSP_Param *hsp_param);
HSP_Param *HSP_Param_share(HSP_Param *hsp_param);
HSP_Param *HSP_Param_swap(HSP_Param *hsp_param);

     void  HSP_Param_set_wordlen(HSP_Param *hsp_param, gint wordlen);
     /**/
     void  HSP_Param_set_dna_hsp_threshold(HSP_Param *hsp_param,
                                           gint dna_hsp_threshold);
     void  HSP_Param_set_protein_hsp_threshold(HSP_Param *hsp_param,
                                               gint protein_hsp_threshold);
     void  HSP_Param_set_codon_hsp_threshold(HSP_Param *hsp_param,
                                             gint protein_hsp_threshold);
     /**/
     void  HSP_Param_set_dna_word_limit(HSP_Param *hsp_param,
                                        gint dna_word_limit);
     void  HSP_Param_set_protein_word_limit(HSP_Param *hsp_param,
                                            gint protein_word_limit);
     void  HSP_Param_set_codon_word_limit(HSP_Param *hsp_param,
                                          gint codon_word_limit);
     /**/
     void  HSP_Param_set_dna_hsp_dropoff(HSP_Param *hsp_param,
                                         gint dna_hsp_dropoff);
     void  HSP_Param_set_protein_hsp_dropoff(HSP_Param *hsp_param,
                                             gint dna_hsp_dropoff);
     void  HSP_Param_set_codon_hsp_dropoff(HSP_Param *hsp_param,
                                           gint dna_hsp_dropoff);
     /**/
     void  HSP_Param_set_hsp_threshold(HSP_Param *hsp_param,
                                       gint hsp_threshold);
     void  HSP_Param_set_seed_repeat(HSP_Param *hsp_param,
                                     gint seed_repeat);

typedef struct HSPset {
              gint    ref_count;
          Sequence    *query;
          Sequence    *target;
         HSP_Param    *param;
              gint ****horizon;
         GPtrArray    *hsp_list;
          /**/
          gboolean     is_finalised;
            PQueue   **filter;
          gboolean     is_empty;
         PQueueSet    *pqueue_set;
} HSPset;

HSPset *HSPset_create(Sequence *query, Sequence *target,
                      HSP_Param *hsp_param);
HSPset *HSPset_share(HSPset *hsp_set);
  void  HSPset_destroy(HSPset *hsp_set);
  void  HSPset_swap(HSPset *hsp_set, HSP_Param *hsp_param);
  void  HSPset_revcomp(HSPset *hsp_set);
  void  HSPset_seed_hsp(HSPset *hsp_set,
                        guint query_start, guint target_start);
  void  HSPset_add_known_hsp(HSPset *hsp_set,
                             guint query_start, guint target_start,
                             guint length);
  void  HSPset_seed_all_hsps(HSPset *hsp_set,
                        guint *seed_list, guint seed_list_len);
/* HSPset_seed_all_hsps() can only be called once on the HSPset.
 * It automatically finalises the HSPset
 * position_list should contain (qpos,tpos) pairs,
 * and may be sorted in place.
 */

HSPset *HSPset_finalise(HSPset *hsp_set);

void HSP_print(HSP *hsp, gchar *name);
void HSPset_print(HSPset *hsp_set);

void HSPset_filter_ungapped(HSPset *hsp_set);
/* Remove HSPs which overlap by more than 50% of the score.
 * This is used for 3:3 ungapped alignments.
 */

/**/

typedef struct HSPset_SList_Node {
    struct HSPset_SList_Node *next;
                        gint  query_pos;
                        gint  target_pos;
} HSPset_SList_Node;

RecycleBin *HSPset_SList_RecycleBin_create(void);

HSPset_SList_Node *HSPset_SList_append(RecycleBin *recycle_bin,
                                       HSPset_SList_Node *next,
                                       gint query_pos, gint target_pos);
void HSPset_seed_all_qy_sorted(HSPset *hsp_set, HSPset_SList_Node *seed_list);

/**/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* INCLUDED_HSPSET_H */

