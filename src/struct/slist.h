/****************************************************************\
*                                                                *
*  Efficient single-linked list routines.                        *
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

/* Note: SList is used instead of the glib GSList object,
 * as it provides:
 *   o All constant time operations (including adding and joining)
 *     ( except SList_traverse() and SList_destroy_with_data() )
 *   o Ease of destroying a set of lists at once
 *   o Ability to add to either end of list.
 *   o SList does not change on append (has a container)
 */

#ifndef INCLUDED_SLIST_H
#define INCLUDED_SLIST_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <glib.h>
#include "recyclebin.h"

typedef struct SListNode {
    struct SListNode *next;
            gpointer  data;
} SListNode;

typedef struct {
    RecycleBin *list_recycle;
    RecycleBin *node_recycle;
} SListSet;

typedef struct {
    SListNode *head;
    SListNode *tail;
     SListSet *set;
} SList;
/* FIXME: optimisation : remove slist->set to save space */

#define SList_isempty(slist) ((slist)->head == (slist)->tail)
#define SList_first(slist) \
        (SList_isempty(slist)?NULL:(slist)->head->next->data)
#define SList_last(slist) \
        (SList_isempty(slist)?NULL:(slist)->tail->data)

typedef void (*SList_DestroyFunc)(gpointer data, gpointer user_data);
typedef gboolean (*SList_TraverseFunc)(gpointer data,
                                       gpointer user_data);
typedef gboolean (*SList_CompareFunc)(gpointer a, gpointer b,
                                      gpointer user_data);
/* SList_TraverseFunc should return TRUE to stop the traversal
 */

SListSet *SListSet_create(void);
    void  SListSet_destroy(SListSet *slist_set);
   gsize  SListSet_memory_usage(SListSet *slist_set);
/* SListSet also destroy will all SLists in this set */

SList *SList_create(SListSet *slist_set);
 void  SList_empty(SList *slist);
 void  SList_destroy(SList *slist);
 void  SList_destroy_with_data(SList *slist,
                               SList_DestroyFunc destroy_func,
                               gpointer user_data);
 void  SList_empty_with_data(SList *slist,
                             SList_DestroyFunc destroy_func,
                             gpointer user_data);

    void SList_queue(SList *slist, gpointer data);
    void SList_stack(SList *slist, gpointer data);
gpointer SList_pop(SList *slist);
  SList *SList_join(SList *left, SList *right);
/* SList_stack() adds to start of list
 * SList_queue() adds to end of list
 * SList_pop() removes from the start of the list
 * SList_join() appends right onto end of left.
 *              right is emptied but not destroyed.
 */

SList *SList_merge(SList *left, SList *right,
                   SList_CompareFunc compare_func,
                   gpointer user_data);
/* SList_merge() assumes left and right are already sorted,
 * and merge and sort the results into left.
 * right is emptied but not destroyed.
 */

void SList_remove(SList *slist, SListNode *prev_node);
void SList_insert(SList *slist, SListNode *prev_node, SList *new_slist);
/* SList_remove() removes prev_node->next from slist
 *                (or the first node if prev_node is NULL)
 * SList_insert() inserts new_slist into slist after prev_node
 *                (or at the beggining if prev_node is NULL)
 *                new_slist is destroyed in the process.
 */

void SList_traverse(SList *slist, SList_TraverseFunc traverse_func,
                    gpointer user_data);

#define SList_for_each(slist, sln) \
    for((sln) = (slist)->head->next; (sln); (sln) = (sln)->next)

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* INCLUDED_SLIST_H */

