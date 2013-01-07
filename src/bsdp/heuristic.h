/****************************************************************\
*                                                                *
*  C4 dynamic programming library - code for heuristics          *
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

#ifndef INCLUDED_HEURISTIC_H
#define INCLUDED_HEURISTIC_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <glib.h>
#include "c4.h"
#include "optimal.h"

#if 0
typedef enum {
    Heuristic_Refinement_NONE,
    Heuristic_Refinement_REGION,
    Heuristic_Refinement_FULL,
    Heuristic_Refinement_TOTAL /* Just the total */
} Heuristic_Refinement;

gchar *Heuristic_Refinement_to_string(Heuristic_Refinement refinement);
Heuristic_Refinement Heuristic_Refinement_from_string(gchar *str);
#endif /* 0 */

typedef struct {
                    gint terminal_range_internal;
                    gint terminal_range_external;
                    gint join_range_internal;
                    gint join_range_external;
                    gint span_range_internal;
                    gint span_range_external;
#if 0
    Heuristic_Refinement refinement;
                    gint refinement_boundary;
#endif /* 0 */
} Heuristic_ArgumentSet;

Heuristic_ArgumentSet *Heuristic_ArgumentSet_create(Argument *arg);

typedef struct {
    C4_Score **matrix;
         gint  query_range;
         gint  target_range;
} Heuristic_Bound;

typedef struct {
    gint internal_query;
    gint internal_target;
    gint external_query;
    gint external_target;
} Heuristic_Range;

typedef struct {
    Heuristic_Range *src_range;
    Heuristic_Range *dst_range;
    C4_DerivedModel *join_model;
            Optimal *optimal;
    Heuristic_Bound *bound;
} Heuristic_Join;

typedef struct {
    gint query_pos;
    gint target_pos;
} Heuristic_Span_Cell;
/* Correspond to positions in the src_region */
/* FIXME: should correspond to positions
 * in the src_integration_matrix
 */

typedef struct {
        Heuristic_Range   *src_range;
        Heuristic_Range   *dst_range;
        C4_DerivedModel   *src_model;
        C4_DerivedModel   *src_traceback_model;
        C4_DerivedModel   *dst_model;
        Heuristic_Bound   *src_bound;
        Heuristic_Bound   *dst_bound;
                Optimal   *src_optimal;
                Optimal   *src_traceback_optimal;
                Optimal   *dst_optimal;
                C4_Span   *span;
/**/
               C4_Score ***src_integration_matrix;
    Heuristic_Span_Cell  **dst_integration_matrix;
                 Region   *curr_src_region;
                 Region   *curr_dst_region;
               C4_Score   *dummy_cell;
} Heuristic_Span;
/* The sizes of the integration matrices are as in {src,dst}_bound
 */

typedef struct {
    Heuristic_Range *range;
    C4_DerivedModel *terminal_model;
    Heuristic_Bound *bound;
            Optimal *optimal;
} Heuristic_Terminal;

typedef struct {
                  gint  id;
             C4_Portal *portal;
         C4_Transition *transition;
    Heuristic_Terminal *start_terminal;
    Heuristic_Terminal *end_terminal;
} Heuristic_Match;
/* A Heuristic_Match for each portal/transition combination */

typedef struct {
    Heuristic_Match *src;
    Heuristic_Match *dst;
     Heuristic_Join *join;
          GPtrArray *span_list; /* Contains (*Heuristic_Span) */
} Heuristic_Pair;
/* A Heuristic_Pair for each valid Heuristic_Match pair */

void Heuristic_Pair_get_max_range(Heuristic_Pair *pair,
                                  gint *max_query_join_range,
                                  gint *max_target_join_range);

typedef struct {
                     gchar   *name;
                      gint    ref_count;
     Heuristic_ArgumentSet   *has;
                  C4_Model   *model;
#if 0
                   Optimal   *optimal;
#endif /* 0 */
                      gint    match_total;
           Heuristic_Match  **match_list;
            Heuristic_Pair ***pair_matrix;
} Heuristic;

/**/

Heuristic *Heuristic_create(C4_Model *model);
     void  Heuristic_destroy(Heuristic *heuristic);
Heuristic *Heuristic_share(Heuristic *heuristic);

/**/

/* FIXME: move to hpair ? */

typedef struct {
    Heuristic_Span *heuristic_span;
} Heuristic_Data;

void Heuristic_Span_register(Heuristic_Span *heuristic_span,
                             Region *src_region, Region *dst_region);
void Heuristic_Span_integrate(Heuristic_Span *heuristic_span,
                              Region *src_region, Region *dst_region);
void Heuristic_Span_add_traceback(Heuristic_Span *heuristic_span,
         Alignment *alignment, gint query_length, gint target_length);

GPtrArray *Heuristic_get_Optimal_list(Heuristic *heuristic);
/* Required for bootstrapper */

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#endif /* INCLUDED_HEURISTIC_H */

