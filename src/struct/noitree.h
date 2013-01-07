/****************************************************************\
*                                                                *
*  NOI_Tree : non-overlapping interval tree                      *
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

#ifndef INCLUDED_NOITREE_H
#define INCLUDED_NOITREE_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "splaytree.h"

typedef struct {
    SplayTree_Set *sts;
} NOI_Tree_Set;

typedef struct NOI_Tree {
           SplayTree *st;
     struct NOI_Tree *delta;
} NOI_Tree;

/**/

NOI_Tree_Set *NOI_Tree_Set_create(gpointer user_data);
        void  NOI_Tree_Set_destroy(NOI_Tree_Set *nts);

/**/

typedef void (*NOI_Tree_Traverse_Func)(gint start, gint length,
                                       gpointer user_data);

NOI_Tree *NOI_Tree_create(NOI_Tree_Set *nts);
    void  NOI_Tree_destroy(NOI_Tree *nt, NOI_Tree_Set *nts);

    void  NOI_Tree_insert(NOI_Tree *nt, NOI_Tree_Set *nts,
                          gint start, gint length);
    void  NOI_Tree_traverse(NOI_Tree *nt, NOI_Tree_Set *nts,
                            NOI_Tree_Traverse_Func ntf);

    void  NOI_Tree_delta_init(NOI_Tree *nt, NOI_Tree_Set *nts);

    void  NOI_Tree_delta_traverse(NOI_Tree *nt, NOI_Tree_Set *nts,
                                  NOI_Tree_Traverse_Func ntf);

/**/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* INCLUDED_NOITREE_H */

