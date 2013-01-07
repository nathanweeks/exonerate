/****************************************************************\
*                                                                *
*  Priority queue library using pairing heaps.                   *
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

#include "pqueue.h"

PQueueSet *PQueueSet_create(void){
    register PQueueSet *pq_set = g_new(PQueueSet, 1);
    pq_set->pq_recycle = RecycleBin_create("PQueue",
                                           sizeof(PQueue), 64);
    pq_set->node_recycle = RecycleBin_create("PQueue_Node",
                                           sizeof(PQueue), 64);
    return pq_set;
    }

void PQueueSet_destroy(PQueueSet *pq_set){
    g_assert(pq_set);
    RecycleBin_destroy(pq_set->pq_recycle);
    RecycleBin_destroy(pq_set->node_recycle);
    g_free(pq_set);
    return;
    }

PQueue *PQueue_create(PQueueSet *pq_set,
                      PQueue_Compare_Func comp_func,
                      gpointer user_data){
    register PQueue *pq;
    g_assert(pq_set);
    g_assert(comp_func);
    pq = RecycleBin_alloc(pq_set->pq_recycle);
    pq->root = NULL;
    pq->total = 0;
    pq->comp_func = comp_func;
    pq->set = pq_set;
    pq->user_data = user_data;
    return pq;
    }

static void PQueue_Node_destroy(PQueueSet *pq_set, PQueue_Node *pqn){
    if(pqn){
        PQueue_Node_destroy(pq_set, pqn->left);
        PQueue_Node_destroy(pq_set, pqn->next);
        RecycleBin_recycle(pq_set->node_recycle, pqn);
        }
    return;
    }

static void PQueue_Node_destroy_with_func(PQueueSet *pq_set,
                   PQueue_Node *pqn,
                   PQueue_Free_Func free_func, gpointer *user_data){
    if(pqn){
        PQueue_Node_destroy_with_func(pq_set, pqn->left,
                                      free_func, user_data);
        PQueue_Node_destroy_with_func(pq_set, pqn->next,
                                      free_func, user_data);
        free_func(pqn->data, user_data);
        RecycleBin_recycle(pq_set->node_recycle, pqn);
        }
    return;
    }

void PQueue_destroy(PQueue *pq, PQueue_Free_Func free_func,
                    gpointer user_data){
    if(free_func)
        PQueue_Node_destroy_with_func(pq->set, pq->root,
                                      free_func, user_data);
    else
        PQueue_Node_destroy(pq->set, pq->root);
    RecycleBin_recycle(pq->set->pq_recycle, pq);
    return;
    }

#ifndef Swap
#define Swap(x,y,temp) ((temp)=(x),(x)=(y),(y)=(temp))
#endif /* Swap */

static PQueue_Node *PQueue_Node_order(PQueue *pq,
                    PQueue_Node *a, PQueue_Node *b){
    register PQueue_Node *tmp;
    g_assert(a);
    g_assert(b);
    if(pq->comp_func(a->data, b->data, pq->user_data)){
        /* Add b before a */
        a->next = b->next;
        if(a->next)
            a->next->prev = a;
        Swap(a, b, tmp);
    } else { /* Add a before b */
        b->prev = a->prev;
        }
    a->prev = b;
    a->next = b->left;
    if(a->next)
        a->next->prev = a;
    b->left = a;
    return b;
    }

PQueue_Node *PQueue_push(PQueue *pq, gpointer data){
    register PQueue_Node *pqn = RecycleBin_alloc(pq->set->node_recycle);
    /* Create empty heap and meld with existing heap */
    pqn->data = data;
    pqn->left = pqn->next = pqn->prev = NULL;
    if(pq->root)
        pq->root = PQueue_Node_order(pq, pq->root, pqn);
    else
        pq->root = pqn; /* Is the first node */
    pq->total++;
    return pqn;
    }

static void PQueue_Node_extract(PQueue_Node *pqn){
    if(pqn->next)
        pqn->next->prev = pqn->prev;
    if(pqn->prev->left == pqn)
        pqn->prev->left = pqn->next;
    else
        pqn->prev->next = pqn->next;
    pqn->next = NULL;
    return;
    }

void PQueue_raise(PQueue *pq, PQueue_Node *pqn){
    if(pq->root == pqn) /* If the root node */
        return;
    PQueue_Node_extract(pqn);
    pq->root = PQueue_Node_order(pq, pq->root, pqn);
    return;
    }

void PQueue_change(PQueue *pq, PQueue_Node *pqn){
    register PQueue_Node *npqn;
    npqn = PQueue_push(pq, PQueue_remove(pq, pqn));
    g_assert(npqn == pqn);
    return;
    }

static PQueue_Node *PQueue_Node_combine(PQueue *pq,
                                        PQueue_Node *pqn){
    register gint i, count;
    register GPtrArray *combine;
    register PQueue_Node *npqn;
    if(!pqn->next){ /* Single member */
        return pqn;
        }
    combine = g_ptr_array_new();
    for(count = 0; pqn; count++){
        g_ptr_array_add(combine, pqn);
        pqn->prev->next = NULL; /* Break links */
        pqn = pqn->next;
        }
    count--;
    for(i = 0; i < count; i+=2)
        combine->pdata[i] = PQueue_Node_order(pq, combine->pdata[i],
                                                  combine->pdata[i+1]);
    if(!(count&1)) /* Pick up spare if odd count */
        combine->pdata[i-2] = PQueue_Node_order(pq, combine->pdata[i-2],
                                                    combine->pdata[i]);
    for(i-=2; i>=2; i-=2)
        combine->pdata[i-2] = PQueue_Node_order(pq, combine->pdata[i-2],
                                                    combine->pdata[i]);
    npqn = combine->pdata[0];
    g_ptr_array_free(combine, TRUE);
    return npqn;
    }

gpointer PQueue_pop(PQueue *pq){
    register PQueue_Node *pqn = pq->root;
    register gpointer data;
    if(!pq->root)
        return NULL; /* Queue is empty */
    data = pq->root->data;
    g_assert((pq->total <= 1) || pq->root->left);
    pq->root = pq->root->left
             ? PQueue_Node_combine(pq, pq->root->left)
             : NULL;
    g_assert(pqn);
    RecycleBin_recycle(pq->set->node_recycle, pqn);
    if(!--pq->total) /* Set root as NULL if PQueue is now empty */
        pq->root = NULL;
    return data;
    }

PQueue *PQueue_join(PQueue *a, PQueue *b){
    g_assert(a->set == b->set);
    g_assert(a->comp_func == b->comp_func);
    g_assert(a->user_data == b->user_data);
    a->root = PQueue_Node_combine(a, b->root);
    a->total += b->total;
    RecycleBin_recycle(a->set->pq_recycle, b);
    return a;
    }

gpointer PQueue_remove(PQueue *pq, PQueue_Node *pqn){
    register gpointer data = pqn->data;
    register PQueue_Node *nb;
    g_assert(pqn);
    g_assert(pq->root);
    if(pqn == pq->root)
        return PQueue_pop(pq);
    g_assert(pq->total > 1);
    PQueue_Node_extract(pqn);
    if(--pq->total){
        if(pqn->left){
            nb = PQueue_Node_combine(pq, pqn->left);
            pq->root = PQueue_Node_order(pq, pq->root, nb);
            }
    } else {
        pq->root = NULL;
        }
    RecycleBin_recycle(pq->set->node_recycle, pqn);
    return data;
    }

static void PQueue_traverse_recur(PQueue *pq, PQueue_Node *pqn,
                       PQueue_Traverse_Func tf, gpointer user_data,
                       gboolean *stop){
    if((!pqn || (*stop)))
        return;
    (*stop) = tf(pqn->data, user_data);
    PQueue_traverse_recur(pq, pqn->left, tf, user_data, stop);
    PQueue_traverse_recur(pq, pqn->next, tf, user_data, stop);
    return;
    }

void PQueue_traverse(PQueue *pq, PQueue_Traverse_Func tf,
                     gpointer user_data){
    gboolean stop = FALSE;
    PQueue_traverse_recur(pq, pq->root, tf, user_data, &stop);
    return;
    }

