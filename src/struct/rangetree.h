/****************************************************************\
*                                                                *
*  Trees for 2D range searching.                                 *
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

#ifndef INCLUDED_RANGETREE_H
#define INCLUDED_RANGETREE_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <glib.h>
#include "recyclebin.h"

typedef struct RangeTree_Node {
    struct RangeTree_Node *left;
    struct RangeTree_Node *right;
                     gint  x;
                     gint  y;
                 gpointer  info;
} RangeTree_Node;

typedef struct {
    RangeTree_Node *root;
              void *recent_data; /* binary search tree */
        RecycleBin *node_recycle;
} RangeTree;

RangeTree *RangeTree_create(void);

typedef void (*RangeTree_FreeFunc)(gpointer info, gpointer user_data);

void RangeTree_destroy(RangeTree *rt, RangeTree_FreeFunc rtff,
                       gpointer user_data);

void RangeTree_add(RangeTree *rt, gint x, gint y, gpointer info);

typedef gboolean (*RangeTree_ReportFunc)(gint x, gint y, gpointer info,
                                         gpointer user_data);

gboolean RangeTree_find(RangeTree *rt, gint x_start, gint x_length,
                                   gint y_start, gint y_length,
                    RangeTree_ReportFunc rtrf, gpointer user_data);
/* Will call rtrf for all points within supplied range.
 * If [xy]_start is zero, will only work on [xy]
 * If rtrf() returns TRUE, will exit early, and return TRUE.
 */

gboolean RangeTree_check_pos(RangeTree *rt, gint x, gint y);
/* Check a single position without refilling RangeTree */

gboolean RangeTree_is_empty(RangeTree *rt);

gboolean RangeTree_traverse(RangeTree *rt,
                            RangeTree_ReportFunc rtrf, gpointer user_data);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* INCLUDED_RANGETREE_H */

