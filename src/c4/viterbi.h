/****************************************************************\
*                                                                *
*  C4 dynamic programming library - the Viterbi implementation   *
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

#ifndef INCLUDED_VITERBI_H
#define INCLUDED_VITERBI_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <glib.h>

#include "c4.h"
#include "region.h"
#include "alignment.h"
#include "argument.h"
#include "layout.h"
#include "codegen.h"
#include "subopt.h"

typedef struct {
        gint traceback_memory_limit;
} Viterbi_ArgumentSet;

Viterbi_ArgumentSet *Viterbi_ArgumentSet_create(Argument *arg);

/**/

typedef struct {
    C4_Score ****score_matrix;              /* 1st element is score */
        gint     cell_size;        /* score, model,built-in shadows */
        gint     region_start_query_id;       /* Cell for qy shadow */
        gint     region_start_target_id;      /* Cell for tg shadow */
        gint     checkpoint_id;              /* Cell for checkpoint */
} Viterbi_Row;
/* Score matrix[max_target_advance]
 *             [query_length]
 *             [state_list->len]
 *             [1 + shadow_total]
 *
 * Cell ids are set to -1 when not in use.
 */

typedef struct {
   GPtrArray *checkpoint_list;
        gint  cell_size;
    C4_Score  last_srp;
        gint  section_length;
        gint  counter;
} Viterbi_Checkpoint;
/* Each checkpoint is a (C4_Score****) as in Viterbi_Row.
 * cell_size is the same as in Viterbi_Row
 */

typedef struct {
    C4_State *first_state; /* Must be set before DP */
    C4_Score *first_cell;  /* Must be set before DP */
    C4_State *final_state; /* Will be set after  DP */
    C4_Score *final_cell;  /* Will be set after  DP */
} Viterbi_Continuation;

typedef struct {
               Viterbi_Row *vr;
/* For Viterbi_Mode_FIND_{REGION,PATH} */
                   gint     curr_query_end;
                   gint     curr_target_end;
/* For Viterbi_Mode_FIND_REGION */
                   gint     curr_query_start;
                   gint     curr_target_start;
                 Region    *alignment_region;
/* For Viterbi_Mode_FIND_PATH */
          C4_Transition ****traceback;
/* For Viterbi_Mode_FIND_CHECKPOINTS */
     Viterbi_Checkpoint    *checkpoint;
   Viterbi_Continuation    *continuation;
} Viterbi_Data;

#define Viterbi_DP_Func_ARGS_STR               \
        "(C4_Model *model, Region *region,     \
          Viterbi_Data *vd, SubOpt_Index *soi, \
          gpointer user_data)"

typedef C4_Score (*Viterbi_DP_Func)( \
        C4_Model *model, Region *region,     \
        Viterbi_Data *vd, SubOpt_Index *soi, \
        gpointer user_data);

/* Avoiding stringify this with #
 * because of gcc 3.1 preprocessor weirdness.
 */

typedef enum {
    Viterbi_Mode_FIND_SCORE,
    Viterbi_Mode_FIND_PATH,
    Viterbi_Mode_FIND_REGION,
    Viterbi_Mode_FIND_CHECKPOINTS
} Viterbi_Mode;

typedef struct {
    Viterbi_ArgumentSet *vas;
                  gchar *name;
               C4_Model *model;
        Viterbi_DP_Func  func;
           Viterbi_Mode  mode;
               gboolean  use_continuation;
                  gsize  cell_size;
                 Layout *layout;
} Viterbi;

Viterbi *Viterbi_create(C4_Model *model, gchar *name,
                        Viterbi_Mode mode, gboolean use_continuation,
                        gboolean use_codegen);
   void  Viterbi_destroy(Viterbi *viterbi);

/* Will issue a warning if a compiled version of a requested
 * function is not available unless codegen is disabled.
 */

gboolean Viterbi_use_reduced_space(Viterbi *viterbi, Region *region);

/**/

Viterbi_Data *Viterbi_Data_create(Viterbi *Viterbi, Region *region);
        void  Viterbi_Data_destroy(Viterbi_Data *vd);
   Alignment *Viterbi_Data_create_Alignment(Viterbi_Data *vd,
                  C4_Model *model, C4_Score score, Region *region);

        void  Viterbi_Data_set_continuation(Viterbi_Data *vd,
                  C4_State *first_state, C4_Score *first_cell,
                  C4_State *final_state, C4_Score *final_cell);
        void  Viterbi_Data_clear_continuation(Viterbi_Data *vd);

typedef struct {
    Region *region;
    C4_State *first_state;
    C4_Score *final_cell;
} Viterbi_SubAlignment;

GPtrArray *Viterbi_Checkpoint_traceback(Viterbi *viterbi,
        Viterbi_Data *vd, Region *region,
        C4_State *first_state, C4_Score *final_cell);
void Viterbi_SubAlignment_destroy(Viterbi_SubAlignment *vsa);

C4_Score Viterbi_calculate(Viterbi *viterbi, Region *region,
                           Viterbi_Data *vd, gpointer user_data,
                           SubOpt *subopt);
Codegen *Viterbi_make_Codegen(Viterbi *viterbi);
/* Creates and compiles codegen */

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* INCLUDED_VITERBI_H */

