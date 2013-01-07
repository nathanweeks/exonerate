/****************************************************************\
*                                                                *
*  C4 dynamic programming library - code for heuristic pairs     *
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

#ifndef INCLUDED_HPAIR_H
#define INCLUDED_HPAIR_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <glib.h>
#include "c4.h"
#include "alignment.h"
#include "bsdp.h"
#include "heuristic.h"
#include "hspset.h"
#include "subopt.h"

typedef struct {
         gint   query_length;
         gint   target_length;
     gpointer   user_data;
    Heuristic  *heuristic;
     gboolean   is_finalised;
    GPtrArray  *portal_data_list; /* HSPset indexed by portal id */
         gint  *bsdp_node_offset; /* indexed by match */
         BSDP  *bsdp;
       SubOpt  *subopt;
         gint   verbosity;
} HPair;

HPair *HPair_create(Heuristic *heuristic, SubOpt *subopt,
                    gint query_length, gint target_length,
                    gint verbosity, gpointer user_data);
void HPair_destroy(HPair *hpair);
/* user_data is provided for the Optimal functions */

void HPair_add_hspset(HPair *hpair, C4_Portal *portal,
                      HSPset *hsp_set);
/* HSPsets must be finalised when provided */

void HPair_finalise(HPair *hpair, C4_Score threshold);

Alignment *HPair_next_path(HPair *hpair, C4_Score threshold);

/* FIXME: need C4_Score HPair_next_score(HPair *hpair);
 */

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* INCLUDED_HPAIR_H */

