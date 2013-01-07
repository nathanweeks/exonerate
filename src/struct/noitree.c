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

/**/

typedef struct {
    gint start;
    gint length;
} NOI_Tree_Interval;

#define NOI_Tree_Interval_end(interval) ((interval)->start+(interval)->length)

static gint NOI_Tree_compare_func(gpointer data_a, gpointer data_b,
                                  gpointer user_data){
    register NOI_Tree_Interval *int_a = data_a, *int_b = data_b;
    if(NOI_Tree_Interval_end(int_a) < int_b->start)
        return -1; /* a before b */
    if(NOI_Tree_Interval_end(int_b) < int_a->start)
        return  1; /* b before a */
    return 0; /* overlapping */
    }

static void NOI_Tree_free_func(gpointer data, gpointer user_data){
    register NOI_Tree_Interval *interval = data;
    g_free(interval);
    return;
    }

typedef struct {
                  gpointer user_data;
    NOI_Tree_Traverse_Func ntf;
} NOI_Tree_Traverse_Data;

NOI_Tree_Set *NOI_Tree_Set_create(gpointer user_data){
    register NOI_Tree_Set *nts = g_new(NOI_Tree_Set, 1);
    register NOI_Tree_Traverse_Data *nttd = g_new(NOI_Tree_Traverse_Data, 1);
    nttd->user_data = user_data;
    nttd->ntf = NULL;
    nts->sts = SplayTree_Set_create(NOI_Tree_compare_func,
                                    NOI_Tree_free_func, nttd);
    return nts;
    }

void NOI_Tree_Set_destroy(NOI_Tree_Set *nts){
    register NOI_Tree_Traverse_Data *nttd = nts->sts->user_data;
    g_free(nttd);
    SplayTree_Set_destroy(nts->sts);
    g_free(nts);
    return;
    }

/**/

NOI_Tree *NOI_Tree_create(NOI_Tree_Set *nts){
    register NOI_Tree *nt = g_new(NOI_Tree, 1);
    nt->st = NULL;
    nt->delta = NULL;
    return nt;
    }

void NOI_Tree_destroy(NOI_Tree *nt, NOI_Tree_Set *nts){
    if(nt->delta)
        NOI_Tree_destroy(nt->delta, nts);
    SplayTree_destroy(nt->st, nts->sts);
    g_free(nt);
    return;
    }

static NOI_Tree_Interval *NOI_Tree_Interval_merge(NOI_Tree_Interval *int_a,
                                                  NOI_Tree_Interval *int_b){
    register gint start = MIN(int_a->start, int_b->start),
                    end = MAX(NOI_Tree_Interval_end(int_a),
                              NOI_Tree_Interval_end(int_b));
    /*
    g_message("MERGING [%d,%d] [%d,%d] ==> [%d,%d]",
              int_a->start, int_a->length,
              int_b->start, int_b->length,
              start, end-start);
    */
    int_a->start = start;
    int_a->length = end-start;
    g_free(int_b);
    return int_a;
    }

void NOI_Tree_insert(NOI_Tree *nt, NOI_Tree_Set *nts,
                     gint start, gint length){
    NOI_Tree_Interval *result = NULL, *prev = NULL;
    register NOI_Tree_Interval *interval = g_new(NOI_Tree_Interval, 1);
    register gint cursor = start;
    /*
    g_message("[%s], (%d,%d) [%p] [%p]", __FUNCTION__, start, length,
                                         nt, interval);
    */
    interval->start = start;
    interval->length = 0;
    /* Splay with zero length interval to bring 1st overlapping node to root
     * This ensures removal of overlapping nodes will be in order
     */
    if(nt->st){
        nt->st = SplayTree_splay(nt->st, nts->sts, interval);
        if(nts->sts->compare_func(interval, nt->st->data,
                                            nts->sts->user_data) > 0){
            if(nt->st->right){ /* remove 1st if before interval */
                nt->st = SplayTree_remove(nt->st, nts->sts, nt->st->data,
                                          (gpointer)&prev);
                cursor = NOI_Tree_Interval_end(prev);
                }
            }
        interval->length = length;
        while(nt->st && (!nts->sts->compare_func(interval, nt->st->data,
                                                  nts->sts->user_data))){
            nt->st = SplayTree_remove(nt->st, nts->sts, nt->st->data,
                                     (gpointer)&result);
            g_assert(result);
            if(nt->delta && (cursor < result->start))
                NOI_Tree_insert(nt->delta, nts, cursor, result->start-cursor);
            cursor = NOI_Tree_Interval_end(result);
            /**/
            interval = NOI_Tree_Interval_merge(interval, result);
            result = NULL;
            }
        if(nt->delta && (cursor < NOI_Tree_Interval_end(interval))){
            NOI_Tree_insert(nt->delta, nts, cursor,
                            NOI_Tree_Interval_end(interval) - cursor);
            }
    } else {
        interval->length = length;
        if(nt->delta)
            NOI_Tree_insert(nt->delta, nts, start, length);
        }
    nt->st = SplayTree_insert(nt->st, nts->sts, interval);
    if(prev)
        nt->st = SplayTree_insert(nt->st, nts->sts, prev);
    return;
    }
/* Delta Algorithm:
 *     start with cursor at start of insert_interval
 *     for each overlapping existing interval
 *          if(cursor < curr_interval->start)
 *             add new_interval from cursort to curr_interval->start to delta_tree
 *          move cursor to end of curr_interval
 *     if(cursor < insert_interval_end)
 *         add new_interval from cursor to insert_interval_end to delta_tree
 */

static gboolean NOI_Tree_SplayTree_traverse_func(gpointer data,
                                                 gpointer user_data){
    register NOI_Tree_Interval *interval = data;
    register NOI_Tree_Traverse_Data *nttd = user_data;
    g_assert(nttd->ntf);
    nttd->ntf(interval->start, interval->length, nttd->user_data);
    return FALSE;
    }

void NOI_Tree_traverse(NOI_Tree *nt, NOI_Tree_Set *nts,
                       NOI_Tree_Traverse_Func ntf){
    register NOI_Tree_Traverse_Data *nttd = nts->sts->user_data;
    g_assert(!nttd->ntf);
    nttd->ntf = ntf;
    SplayTree_traverse(nt->st, nts->sts, NOI_Tree_SplayTree_traverse_func);
    nttd->ntf = NULL;
    return;
    }

void NOI_Tree_delta_init(NOI_Tree *nt, NOI_Tree_Set *nts){
    if(nt->delta)
        NOI_Tree_destroy(nt->delta, nts);
    nt->delta = NOI_Tree_create(nts);
    return;
    }

void NOI_Tree_delta_traverse(NOI_Tree *nt, NOI_Tree_Set *nts,
                             NOI_Tree_Traverse_Func ntf){
    if(nt->delta)
        NOI_Tree_traverse(nt->delta, nts, ntf);
    return;
    }

/**/

