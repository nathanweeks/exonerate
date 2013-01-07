/****************************************************************\
*                                                                *
*  C4 dynamic programming library - suboptimal alignments        *
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

#ifndef INCLUDED_SUBOPT_H
#define INCLUDED_SUBOPT_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <glib.h>

#include "alignment.h"
#include "region.h"
#include "rangetree.h"

/* Points are stored in a RangeTree,
 * then put in row/col lists for dp lookup.
 */

typedef struct {
          gint  ref_count;
          gint  query_length;
          gint  target_length;
     RangeTree *range_tree;
          gint  path_count;
} SubOpt;

SubOpt *SubOpt_create(gint query_length, gint target_length);
SubOpt *SubOpt_share(SubOpt *subopt);
  void  SubOpt_destroy(SubOpt *subopt);
  void  SubOpt_add_alignment(SubOpt *subopt, Alignment *alignment);

typedef gboolean (*SubOpt_FindFunc)(gint query_pos, gint target_pos,
                                    gint path_id, gpointer user_data);

gboolean SubOpt_find(SubOpt *subopt, Region *region,
                     SubOpt_FindFunc find_func, gpointer user_data);
gboolean SubOpt_overlaps_alignment(SubOpt *subopt, Alignment *alignment);

/**/

typedef struct {
    gint *query_pos;
    gint  target_pos;
    gint  total;       /* Number of query positions in this row */
} SubOpt_Index_Row;

typedef struct {
              Region *region;
           GPtrArray *row_list; /* Contains SubOpt_Index_Row */
    SubOpt_Index_Row *curr_row;
                gint  curr_row_index;
                gint  curr_query_index;
    SubOpt_Index_Row *blank_row;
} SubOpt_Index;

SubOpt_Index *SubOpt_Index_create(SubOpt *subopt, Region *region);
        void  SubOpt_Index_destroy(SubOpt_Index *soi);
        void  SubOpt_Index_set_row(SubOpt_Index *soi, gint target_pos);
/* FIXME: could make SubOpt_Index_set_row a macro for speed */
    gboolean  SubOpt_Index_is_blocked(SubOpt_Index *soi,
                                      gint query_pos);

#define SubOpt_Index_is_blocked_fast(soi, q_pos)                 \
   ((soi->curr_row->query_pos[  soi->curr_query_index] <  q_pos) \
   ?(soi->curr_row->query_pos[++soi->curr_query_index] == q_pos) \
   :(soi->curr_row->query_pos[  soi->curr_query_index] == q_pos))
/* For SubOpt_Index_set_row() and SubOpt_Index_is_blocked_fast(),
 * query_pos and target_pos should use the coordinates
 * of the indexed region, not of the whole sequence
 *
 * It is assumed that target_pos increments between
 * each call to SubOpt_Index_set_row().
 *
 * Multiple calls to SubOpt_Index_is_blocked() must be allowed,
 * for models with multiple match states.
 */

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* INCLUDED_SUBOPT_H */

