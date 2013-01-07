/****************************************************************\
*                                                                *
*  C4 dynamic programming library - Seeded Dynamic Programming   *
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

#ifndef INCLUDED_SDP_H
#define INCLUDED_SDP_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <glib.h>

#include "c4.h"
#include "comparison.h"
#include "slist.h"
#include "recyclebin.h"
#include "alignment.h"
#include "lookahead.h"
#include "subopt.h"
#include "region.h"
#include "sparsecache.h"
#include "boundary.h"
#include "scheduler.h"
#include "threadref.h"

/**/

typedef struct {
        gint dropoff;
    gboolean single_pass_subopt;
} SDP_ArgumentSet;

SDP_ArgumentSet *SDP_ArgumentSet_create(Argument *arg);

/**/

typedef struct {
                gint  query_pos;
                gint  target_pos;
            C4_Score  score;
     STraceback_Cell *cell;
} SDP_Terminal;

typedef struct SDP_Seed {
                gint  seed_id;
                 HSP *hsp;
        SDP_Terminal *max_start;     /* Best start from this HSP */
        SDP_Terminal *max_end;       /* Best end from this HSP   */
         PQueue_Node *pq_node;       /* To allow raising in PQ   */
} SDP_Seed;

/* first_pre->*->last_pre->[INDEX]->first_post->*->last_post
 */

/**/

typedef struct {
          ThreadRef *thread_ref;
    SDP_ArgumentSet *sas;
           C4_Model *model;
           gboolean  use_boundary;
          Scheduler *find_starts_scheduler;
          Scheduler *find_ends_scheduler;
} SDP;

SDP *SDP_create(C4_Model *model);
SDP *SDP_share(SDP *sdp);
void SDP_destroy(SDP *sdp);
GPtrArray *SDP_get_codegen_list(SDP *sdp);

/**/

typedef struct {
              SDP *sdp;
       Comparison *comparison;
             gint  alignment_count; /* Number reported so far */
         gpointer  user_data;
           SubOpt *subopt;
        GPtrArray *seed_list;
         gpointer *seed_list_by_score; /* Sorted in score order */
             gint  single_pass_pos;
         Boundary *boundary;
         C4_Score  last_score;
       STraceback *fwd_straceback;
       STraceback *rev_straceback;
} SDP_Pair;

SDP_Pair *SDP_Pair_create(SDP *sdp, SubOpt *subopt,
                          Comparison *comparison, gpointer user_data);
void SDP_Pair_destroy(SDP_Pair *sdp_pair);

/**/

Alignment *SDP_Pair_next_path(SDP_Pair *sdp_pair, C4_Score threshold);

/**/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* INCLUDED_SDP_H */

