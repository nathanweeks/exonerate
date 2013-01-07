/****************************************************************\
*                                                                *
*  Splay Tree Library                                            *
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

/* Splay Trees
 *
 * For a description, see:
 * "Self-adjusting Binary Search Trees"
 * D.D.Sleator and R.E.Tarjan
 * JACM 32:3 pp 652-686, July 1985
 */

#ifndef INCLUDED_SPLAYTREE_H
#define INCLUDED_SPLAYTREE_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <glib.h>
#include "recyclebin.h"

typedef struct SplayTree {
    struct SplayTree *left;
    struct SplayTree *right;
            gpointer  data;
} SplayTree;

typedef gint (*SplayTree_CompareFunc)(gpointer data_a, gpointer data_b,
                                      gpointer user_data);
/* Returns 0 if equal or overlapping, <0 if before, >0 if after */
typedef void (*SplayTree_FreeFunc)(gpointer data, gpointer user_data);

typedef struct {
               RecycleBin *recycle;
    SplayTree_CompareFunc  compare_func;
       SplayTree_FreeFunc  free_func;
                 gpointer  user_data;
} SplayTree_Set;

typedef gboolean (*SplayTree_Traverse_Func)(gpointer data,
                                            gpointer user_data);

SplayTree_Set *SplayTree_Set_create(SplayTree_CompareFunc compare_func,
                                    SplayTree_FreeFunc free_func,
                                    gpointer user_data);
         void  SplayTree_Set_destroy(SplayTree_Set *sts);
        gsize  SplayTree_Set_memory_usage(SplayTree_Set *sts);

SplayTree *SplayTree_splay(SplayTree *st, SplayTree_Set *sts,
                           gpointer data);
SplayTree *SplayTree_insert(SplayTree *st, SplayTree_Set *sts,
                            gpointer data);
SplayTree *SplayTree_remove(SplayTree *st, SplayTree_Set *sts,
                            gpointer data, gpointer *result);
SplayTree *SplayTree_lookup(SplayTree *st, SplayTree_Set *sts,
                            gpointer data);
      void SplayTree_destroy(SplayTree *st, SplayTree_Set *sts);
      void SplayTree_traverse(SplayTree *st, SplayTree_Set *sts,
                              SplayTree_Traverse_Func traverse_func);

/**/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* INCLUDED_SPLAYTREE_H */

