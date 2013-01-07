/****************************************************************\
*                                                                *
*  C4 dynamic programming library - DP Scheduler                 *
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

#ifndef INCLUDED_SCHEDULER_H
#define INCLUDED_SCHEDULER_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <glib.h>
#include "c4.h"
#include "region.h"
#include "lookahead.h"
#include "boundary.h"
#include "slist.h"
#include "recyclebin.h"
#include "subopt.h"
#include "sparsecache.h"
#include "codegen.h"
#include "straceback.h"

/**/

typedef struct {
         gint query_pos;
         gint target_pos;
         gint seed_id;
     C4_Score start_score;
} Scheduler_Seed;

typedef void (*Scheduler_Seed_init_Func)(gpointer seed_data);
typedef void (*Scheduler_Seed_next_Func)(gpointer seed_data);
typedef gboolean (*Scheduler_Seed_get_Func)(gpointer seed_data,
                                            Scheduler_Seed *seed);
/*
typedef void(*Scheduler_obstruct_Func)(gpointer seed_data,
                                       gint src_seed_id,
                                       gint dst_seed_id);
*/
typedef void(*Scheduler_start_Func)(gpointer seed_data,
                                    gint seed_id, C4_Score score,
                                    gint start_query_pos,
                                    gint start_target_pos,
                                    STraceback_Cell *stcell);
typedef void(*Scheduler_end_Func)(gpointer seed_data,
                                  gint seed_id, C4_Score score,
                                  gint end_query_pos,
                                  gint end_target_pos,
                                  STraceback_Cell *stcell);

struct Scheduler_Pair; /* To avoid compiler warning */
struct Scheduler_Row;  /* To avoid compiler warning */
struct Scheduler_Cell; /* To avoid compiler warning */
typedef void (*Scheduler_DP_Func)(struct Scheduler_Cell *cell,
                                  struct Scheduler_Row *row,
                                  struct Scheduler_Pair *spair);
#define Scheduler_DP_Func_ARGS_STR      \
        "(struct Scheduler_Cell *cell,  \
          struct Scheduler_Row  *row,   \
          struct Scheduler_Pair *pair)" \

typedef struct {
                    gboolean   is_forward;
                    gboolean   has_traceback;
                    gboolean   use_boundary;
                        gint   shadow_start;
                       gsize   cell_score_offset;
                       gsize   cell_traceback_offset;
                       gsize   cell_size;
    Scheduler_Seed_init_Func   init_func;
    Scheduler_Seed_next_Func   next_func;
     Scheduler_Seed_get_Func   get_func;
     /* Scheduler_obstruct_Func   obstruct_func; */
        Scheduler_start_Func   start_func;
          Scheduler_end_Func   end_func;
                    C4_Model  *model;
                     C4_Span **span_map;
                    C4_Score   dropoff;
                       gchar  *name;
           Scheduler_DP_Func   cell_func;
} Scheduler;

Scheduler *Scheduler_create(C4_Model *model,
                      gboolean is_forward, gboolean has_traceback,
                      gboolean use_boundary,
                      Scheduler_Seed_init_Func init_func,
                      Scheduler_Seed_next_Func next_func,
                      Scheduler_Seed_get_Func get_func,
                      /* Scheduler_obstruct_Func obstruct_func, */
                      Scheduler_start_Func start_func,
                      Scheduler_end_Func end_func,
                      C4_Score dropoff);
void Scheduler_destroy(Scheduler *scheduler);
Codegen *Scheduler_make_Codegen(Scheduler *scheduler);

/**/

typedef struct {
           C4_Score  score;
           C4_Score  max;
               gint  seed_id;
               gint  query_entry;
               gint  target_entry;
    STraceback_Cell *cell;
           C4_Score *shadow_data;
} Scheduler_SpanSeed;

typedef struct {
               C4_Span *span;
    Scheduler_SpanSeed *curr_span_seed;
} Scheduler_SpanData;

typedef struct {
    Scheduler_SpanData **span_data;
              C4_Model  *model;      /* FIXME: move to Scheduler ? */
} Scheduler_SpanHistory;

/**/

typedef struct Scheduler_Pair {
                Scheduler  *scheduler;
                 SListSet  *slist_set;
                   Region  *region;
                     gint   best_seed_id;
                 gpointer   user_data;
                 gpointer   seed_data;
                Lookahead  *row_index;
               RecycleBin  *cell_recycle;
               STraceback  *straceback;
/**/
                 Boundary  *boundary;
              SparseCache  *span_seed_cache;
    Scheduler_SpanHistory  *span_history;
/**/
             SubOpt_Index  *subopt_index;
} Scheduler_Pair;

Scheduler_Pair *Scheduler_Pair_create(
                           Scheduler *scheduler,
                           STraceback *straceback,
                           gint query_len, gint target_len,
                           SubOpt *subopt,
                           Boundary *boundary, gint best_seed_id,
                           gpointer seed_data, gpointer user_data);

void Scheduler_Pair_destroy(Scheduler_Pair *spair);
void Scheduler_Pair_calculate(Scheduler_Pair *spair);

/**/

typedef struct Scheduler_Cell {
               gint   query_pos;
           gboolean   permit_span_thaw; /* FIXME: remove */
           C4_Score **score;
    STraceback_Cell **traceback; /* Only when performing traceback */
} Scheduler_Cell;
/* score matrix[state->list->len  ]
 *             [ 4 + shadow_total ]
 *      all: 0:[score]
 *      all: 1:[max]
 *      all: 2:[seed_id]
 * {other shadows ...}
 *
 * traceback[state->list->len]
 *
 */

/**/

typedef struct Scheduler_Row {
        gint   target_pos;
   Lookahead  *cell_index;
       SList  *used_list;
       SList  *unused_list;
} Scheduler_Row;

/* FIXME: BEGIN temp for codegen */
void Scheduler_Cell_process(Scheduler_Cell *cell,
                            Scheduler_Row *row,
                            Scheduler_Pair *spair);
Scheduler_Row *Scheduler_Row_create(gint target_pos,
                               Scheduler_Pair *spair);
Scheduler_Cell *Scheduler_Cell_create(gint query_pos,
                                      gboolean permit_span_thaw,
                                      Scheduler_Pair *spair);
void Scheduler_Cell_assign(Scheduler_Pair *spair,
                           Scheduler_Cell *src_cell, gint input_pos,
                           Scheduler_Cell *dst_cell, gint output_pos,
                           C4_Score dst_score, C4_Score max_score,
                           C4_Transition *transition, gint seed_id,
                           gint dst_query_pos, gint dst_target_pos);
void Scheduler_SpanData_get_curr(Scheduler_SpanData *span_data,
                                 SparseCache *cache,
                                 gint query_pos, gint target_pos,
                                 STraceback *straceback);
void Scheduler_SpanData_submit(Scheduler_SpanData *span_data,
                                Scheduler_SpanSeed *span_seed,
                                SparseCache *cache,
                                Scheduler_Pair *spair);
STraceback_Cell *Scheduler_Cell_add_span(STraceback_Cell *prev,
                                         STraceback *straceback,
                                         C4_Span *span,
                                         gint query_len, gint target_len);
/* FIXME: END temp for codegen */

/**/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* INCLUDED_SCHEDULER_H */

