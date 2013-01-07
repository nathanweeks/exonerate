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

#include "noitree.h"

static void test_traverse_func(gint start, gint length,
                                gpointer user_data){
    g_message("traverse with [%d,%d]", start, length);
    return;
    }

int main(void){
    register NOI_Tree_Set *nts = NOI_Tree_Set_create(NULL);
    register NOI_Tree *nt = NOI_Tree_create(nts);
    /**/
#if 0
    NOI_Tree_insert(nt, nts,  1, 3);
    NOI_Tree_insert(nt, nts,  7, 3);
    NOI_Tree_insert(nt, nts, 13, 3);
    NOI_Tree_insert(nt, nts, 19, 3);
    NOI_Tree_insert(nt, nts, 25, 3);
    g_message("traverse 1:");
    NOI_Tree_traverse(nt, nts, test_traverse_func);
    /**/
    NOI_Tree_delta_init(nt, nts);
    NOI_Tree_insert(nt, nts,  6, 19);
    g_message("traverse 2:");
    NOI_Tree_traverse(nt, nts, test_traverse_func);
    g_message("traverse delta:");
    NOI_Tree_delta_traverse(nt, nts, test_traverse_func);
    /**/
    NOI_Tree_destroy(nt, nts);
#endif /* 0 */

    NOI_Tree_delta_init(nt, nts);
    NOI_Tree_insert(nt, nts, 100, 10);
    NOI_Tree_insert(nt, nts, 300, 10);
    NOI_Tree_insert(nt, nts, 200, 10);
    g_message("traverse:");
    NOI_Tree_traverse(nt, nts, test_traverse_func);
    g_message("delta traverse:");
    NOI_Tree_delta_traverse(nt, nts, test_traverse_func);

    /**/
    NOI_Tree_Set_destroy(nts);
    g_message("done");
    return 0;
    }

/**/

