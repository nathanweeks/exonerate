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

#include "slist.h"

/**/

SListSet *SListSet_create(void){
    register SListSet *slist_set = g_new(SListSet, 1);
    slist_set->list_recycle = RecycleBin_create("SList",
                                                sizeof(SList), 256);
    slist_set->node_recycle = RecycleBin_create("SListNode",
                                           sizeof(SListNode), 1024);
    return slist_set;
    }

void SListSet_destroy(SListSet *slist_set){
    g_assert(slist_set);
    RecycleBin_destroy(slist_set->list_recycle);
    RecycleBin_destroy(slist_set->node_recycle);
    g_free(slist_set);
    return;
    }

gsize SListSet_memory_usage(SListSet *slist_set){
    return sizeof(SListSet)
         + RecycleBin_memory_usage(slist_set->list_recycle)
         + RecycleBin_memory_usage(slist_set->node_recycle);
    }

/**/

static SListNode *SListNode_alloc(SListSet *slist_set){
    register SListNode *sln = RecycleBin_alloc(slist_set->node_recycle);
    sln->next = NULL;
    sln->data = NULL;
    return sln;
    }

static void SListNode_free(SListSet *slist_set, SListNode *sln){
    RecycleBin_recycle(slist_set->node_recycle, sln);
    return;
    }

/**/

SList *SList_create(SListSet *slist_set){
    register SList *slist;
    g_assert(slist_set);
    slist = RecycleBin_alloc(slist_set->list_recycle);
    slist->head = slist->tail = SListNode_alloc(slist_set);
    slist->set = slist_set;
    return slist;
    }

void SList_empty(SList *slist){
    g_assert(slist);
    while(!SList_isempty(slist))
        SList_pop(slist);
    return;
    }

void SList_destroy(SList *slist){
    g_assert(slist);
    SList_empty(slist);
    RecycleBin_recycle(slist->set->list_recycle, slist);
    return;
    }

void SList_empty_with_data(SList *slist,
                           SList_DestroyFunc destroy_func,
                           gpointer user_data){
    register gpointer *data;
    g_assert(slist);
    g_assert(destroy_func);
    while(!SList_isempty(slist)){
        data = SList_pop(slist);
        destroy_func(data, user_data);
        }
    return;
    }

void SList_destroy_with_data(SList *slist,
                             SList_DestroyFunc destroy_func,
                             gpointer user_data){
    g_assert(slist);
    g_assert(destroy_func);
    SList_empty_with_data(slist, destroy_func, user_data);
    SList_destroy(slist);
    return;
    }

void SList_queue(SList *slist, gpointer data){
    g_assert(slist);
    slist->tail->next = SListNode_alloc(slist->set);
    slist->tail = slist->tail->next;
    slist->tail->data = data;
    return;
    }

void SList_stack(SList *slist, gpointer data){
    register SListNode *sln;
    g_assert(slist);
    sln = SListNode_alloc(slist->set);
    slist->head->data = data;
    sln->next = slist->head;
    slist->head = sln;
    return;
    }

gpointer SList_pop(SList *slist){
    register gpointer data;
    register SListNode *sln;
    g_assert(slist);
    data = slist->head->next->data;
    sln = slist->head;
    slist->head = sln->next;
    SListNode_free(slist->set, sln);
    return data;
    }

SList *SList_join(SList *left, SList *right){
    register SListNode *head, *tail;
    g_assert(left);
    g_assert(right);
    g_assert(left->set == right->set);
    if(SList_isempty(right))
        return left;
    if(SList_isempty(left)){
        head = left->head;
        tail = left->tail;
        left->head = right->head;
        left->tail = right->tail;
        right->head = head;
        right->tail = tail;
        return left;
        }
    g_assert(!SList_isempty(right));
    left->tail->next = right->head->next;
    left->tail = right->tail;
    right->tail = right->head;
    right->head->next = NULL;
    return left;
    }

SList *SList_merge(SList *left, SList *right,
                   SList_CompareFunc compare_func, gpointer user_data){
    register SList *merged = SList_create(left->set);
    register SList swap;
    g_assert(left->set == right->set);
    while(!(SList_isempty(left) || SList_isempty(right))){
        if(compare_func(left->head->next->data,
                        right->head->next->data, user_data)){
            SList_queue(merged, SList_pop(left));
        } else {
            SList_queue(merged, SList_pop(right));
            }
        }
    SList_join(merged, left);
    SList_join(merged, right);
    SList_destroy(left);
    /* Swap left and merged */
    swap.head = left->head;
    swap.tail = left->tail;
    left->head = merged->head;
    left->tail = merged->tail;
    merged->head = swap.head;
    merged->tail = swap.tail;
    return merged;
    }
/* FIXME: optimisation : more efficient without the pop/queue */

void SList_remove(SList *slist, SListNode *prev_node){
    register SListNode *sln;
    g_assert(slist);
    if(!prev_node){
        SList_pop(slist);
        return;
        }
    if(prev_node->next){
        sln = prev_node->next;
        prev_node->next = prev_node->next->next;
        SListNode_free(slist->set, sln);
        }
    return;
    }

void SList_insert(SList *slist, SListNode *prev_node, SList *new_slist){
    register SListNode *sln;
    g_assert(slist);
    g_assert(new_slist);
    g_assert(slist->set == new_slist->set);
    if(!prev_node)
        prev_node = slist->head;
    sln = prev_node->next;
    if(!sln)
        slist->tail = new_slist->tail;
    prev_node->next = new_slist->head->next;
    new_slist->tail->next = sln;
    SListNode_free(new_slist->set, new_slist->head);
    RecycleBin_recycle(new_slist->set->list_recycle, new_slist);
    return;
    }

void SList_traverse(SList *slist, SList_TraverseFunc traverse_func,
                    gpointer user_data){
    register SListNode *sln;
    g_assert(slist);
    g_assert(traverse_func);
    for(sln = slist->head->next; sln; sln = sln->next)
        if(traverse_func(sln->data, user_data))
            break;
    return;
    }

