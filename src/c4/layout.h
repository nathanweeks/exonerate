/****************************************************************\
*                                                                *
*  C4 dynamic programming library - DP layout code               *
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

#ifndef INCLUDED_LAYOUT_H
#define INCLUDED_LAYOUT_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <glib.h>

#include "c4.h"

/**/

typedef struct {
    GArray *transition_mask;
} Layout_Mask;

typedef struct {
    Layout_Mask *normal;
    Layout_Mask *end_query;
    Layout_Mask *end_target;
    Layout_Mask *corner;
} Layout_Cell;

typedef struct {
    GPtrArray *cell_list;
} Layout_Row;

typedef struct {
    GPtrArray *row_list;
} Layout;

/**/

Layout *Layout_create(C4_Model *model);
  void  Layout_destroy(Layout *layout);
  void  Layout_info(Layout *layout, C4_Model *model);

gboolean Layout_is_transition_valid(Layout *layout, C4_Model *model,
                              C4_Transition *transition,
                              gint dst_query_pos, gint dst_target_pos,
                              gint query_length, gint target_length);

/**/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* INCLUDED_LAYOUT_H */
