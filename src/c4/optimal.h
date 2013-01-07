/****************************************************************\
*                                                                *
*  C4 dynamic programming library - optimal alignment code       *
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

#ifndef INCLUDED_OPTIMAL_H
#define INCLUDED_OPTIMAL_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <glib.h>

#include "viterbi.h"
#include "subopt.h"

/**/

typedef enum {
    Optimal_Type_SCORE         = (1 << 0),
    Optimal_Type_PATH          = (1 << 1),
    Optimal_Type_REDUCED_SPACE = (1 << 2)
} Optimal_Type;

typedef struct {
          gchar *name;
           gint  ref_count;
   Optimal_Type   type;
        Viterbi  *find_score;
        Viterbi  *find_path;
        Viterbi  *find_region;
        Viterbi  *find_checkpoint_continuation;
        Viterbi  *find_path_continuation;
} Optimal;

/**/

Optimal *Optimal_create(C4_Model *model, gchar *name,
                        Optimal_Type type, gboolean use_codegen);

/* If name is NULL, "optimal:model->name" is used */
   void  Optimal_destroy(Optimal *optimal);
Optimal *Optimal_share(Optimal *optimal);

GPtrArray *Optimal_make_Codegen_list(Optimal *optimal);
/* Returns a list of the Codegen objects created */

C4_Score Optimal_find_score(Optimal *optimal, Region *region,
                            gpointer user_data, SubOpt *subopt);
/* Subopt is required for Optimal_find_score() in BSDP */

Alignment *Optimal_find_path(Optimal *optimal, Region *region,
                             gpointer user_data, C4_Score threshold,
                             SubOpt *subopt);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* INCLUDED_OPTIMAL_H */

