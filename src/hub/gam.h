/****************************************************************\
*                                                                *
*  GAM: Gapped Alignment Manager                                 *
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

#ifndef INCLUDED_GAM_H
#define INCLUDED_GAM_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <glib.h>

#ifdef USE_PTHREADS
#include <pthread.h>
#endif /* USE_PTHREADS */

#include "sequence.h"
#include "c4.h"
#include "heuristic.h"
#include "hpair.h"
#include "sequence.h"
#include "hspset.h"
#include "submat.h"
#include "argument.h"
#include "translate.h"
#include "pqueue.h"
#include "modeltype.h"
#include "comparison.h"
#include "sdp.h"
#include "subopt.h"
#include "threadref.h"

typedef enum {
    GAM_Refinement_NONE,
    GAM_Refinement_REGION,
    GAM_Refinement_FULL,
    GAM_Refinement_TOTAL /* Just the total */
} GAM_Refinement;

typedef struct {
        Model_Type  type;
          C4_Score  threshold;
            gfloat  percent_threshold;
          gboolean  show_alignment;
          gboolean  show_sugar;
          gboolean  show_cigar;
          gboolean  show_vulgar;
          gboolean  show_query_gff;
          gboolean  show_target_gff;
             gchar *ryo;
              gint  best_n;
          gboolean  use_subopt;
          gboolean  use_gapped_extension;
      /**/
    GAM_Refinement  refinement;
              gint  refinement_boundary;
} GAM_ArgumentSet;

GAM_ArgumentSet *GAM_ArgumentSet_create(Argument *arg);

typedef struct {
       C4_Score score;
          glong pos;
          glong len;
} GAM_StoredResult;

typedef struct {
     gchar *query_id;
    PQueue *pq; /* Contains GAM_StoredResult */
      gint  tie_count; /* For best_n */
  C4_Score  tie_score; /* For best_n */
} GAM_QueryResult;

typedef struct {
       gchar *query_id;
    C4_Score threshold;
} GAM_QueryInfo;

typedef struct {
             ThreadRef *thread_ref;
         Alphabet_Type  query_type;
         Alphabet_Type  target_type;
              C4_Model *model;
             GPtrArray *match_list;
               Optimal *optimal;
             Heuristic *heuristic;
                   SDP *sdp;
                Submat *dna_submat;
                Submat *protein_submat;
             Translate *translate;
       GAM_ArgumentSet *gas;
                  void *bestn_tree; /* Contains GAM_QueryResult */
                  FILE *bestn_tmp_file;
                  gint  verbosity;
              gboolean  translate_both;
              gboolean  dual_match;
                  void *percent_threshold_tree; /* Contains GAM_QueryInfo */
             PQueueSet *pqueue_set;
                  gint  max_query_span;
                  gint  max_target_span;
#ifdef USE_PTHREADS
       pthread_mutex_t  gam_lock;
#endif /* USE_PTHREADS */
} GAM;

GAM *GAM_create(Alphabet_Type query_type, Alphabet_Type target_type,
                Submat *dna_submat, Submat *protein_submat,
                Translate *translate, gboolean use_exhaustive,
                gint verbosity);
GAM *GAM_share(GAM *gam);
void GAM_destroy(GAM *gam);
void GAM_report(GAM *gam);

typedef struct {
         gint  ref_count;
          GAM *gam;
     Sequence *query;
     Sequence *target;
    GPtrArray *alignment_list;
     gpointer  user_data;
     gpointer  self_data;
       SubOpt *subopt;
} GAM_Result;

GAM_Result *GAM_Result_ungapped_create(GAM *gam,
                                       Comparison *comparison);
/* Will return NULL when no alignments are produced */

GAM_Result *GAM_Result_heuristic_create(GAM *gam,
                                        Comparison *comparison);
/* Will return NULL when no alignments are produced */

GAM_Result *GAM_Result_exhaustive_create(GAM *gam,
                                         Sequence *query,
                                         Sequence *target);

GAM_Result *GAM_Result_share(GAM_Result *gam_result);
      void  GAM_Result_destroy(GAM_Result *gam_result);
      void  GAM_Result_submit(GAM_Result *gam_result);

void GAM_lock(GAM *gam);
void GAM_unlock(GAM *gam);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* INCLUDED_GAM_H */

