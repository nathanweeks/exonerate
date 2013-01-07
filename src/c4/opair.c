/****************************************************************\
*                                                                *
*  C4 dynamic programming library - code for optimal pairs       *
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

#include "opair.h"

/**/

OPair *OPair_create(Optimal *optimal, SubOpt *subopt,
                    gint query_length, gint target_length, gpointer user_data){
    register OPair *opair = g_new(OPair, 1);
    g_assert(optimal);
    g_assert(query_length >= 0);
    g_assert(target_length >= 0);
    opair->optimal = Optimal_share(optimal);
    opair->user_data = user_data;
    opair->subopt = SubOpt_share(subopt);
    opair->region = Region_create(0, 0, query_length, target_length);
    opair->prev_score = C4_IMPOSSIBLY_HIGH_SCORE;
    return opair;
    }

void OPair_destroy(OPair *opair){
    SubOpt_destroy(opair->subopt);
    Optimal_destroy(opair->optimal);
    Region_destroy(opair->region);
    g_free(opair);
    return;
    }

Alignment *OPair_next_path(OPair *opair, C4_Score threshold){
    register Alignment *alignment;
    alignment = Optimal_find_path(opair->optimal, opair->region,
                                  opair->user_data, threshold,
                                  opair->subopt);
    if(alignment){
        if(alignment->score > opair->prev_score)
            g_warning("Score [%d] greater than previous [%d]",
                      alignment->score, opair->prev_score);
        g_assert(alignment->score <= opair->prev_score);
        opair->prev_score = alignment->score;
        }
    return alignment;
    }

/**/

