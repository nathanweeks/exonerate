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

#include <stdlib.h> /* For lrand48() */

#include "rangetree.h"

static gint RangeTree_recent_data_compare(gconstpointer a,
                                          gconstpointer b){
    register RangeTree_Node *rtn_a = (RangeTree_Node*)a,
                            *rtn_b = (RangeTree_Node*)b;
    if(rtn_a->x == rtn_b->x)
        return rtn_b->y - rtn_a->y;
    return rtn_b->x - rtn_a->x;
    }

RangeTree *RangeTree_create(void){
    register RangeTree *rt = g_new(RangeTree, 1);
    rt->root = NULL;
    rt->recent_data = g_tree_new(RangeTree_recent_data_compare);
    rt->node_recycle = RecycleBin_create("RangeTree_Node",
                                         sizeof(RangeTree_Node), 256);
    return rt;
    }

void RangeTree_destroy(RangeTree *rt, RangeTree_FreeFunc rtff,
                              gpointer user_data){
    g_tree_destroy(rt->recent_data);
    RecycleBin_destroy(rt->node_recycle);
    g_free(rt);
    return;
    }

void RangeTree_add(RangeTree *rt, gint x, gint y, gpointer info){
    register RangeTree_Node *rtn = RecycleBin_alloc(rt->node_recycle);
    g_assert(!RangeTree_check_pos(rt, x, y));
    rtn->x = x;
    rtn->y = y;
    rtn->info = info;
    rtn->left = NULL;
    rtn->right = NULL;
    g_assert(!g_tree_lookup(rt->recent_data, rtn));
    g_tree_insert(rt->recent_data, rtn, rtn);
    return;
    }

typedef struct {
    gint tl_x;
    gint tl_y;
    gint br_x;
    gint br_y;
} RangeTree_Rectangle;

static gboolean RangeTree_inside_rectangle(RangeTree_Rectangle *r,
                                           RangeTree_Node *n){
    if((n->x < r->tl_x)
    || (n->y < r->tl_y)
    || (n->x >= r->br_x)
    || (n->y >= r->br_y))
        return FALSE;
    return TRUE;
    }

static void RangeTree_find_recur(RangeTree_Node *n, gboolean dir,
                                 RangeTree_Rectangle *r,
                                 RangeTree_ReportFunc rtrf,
                                 gpointer user_data, gboolean *found){
    if(!n)
        return;
    if(dir ? (r->tl_x < n->x) : (r->tl_y < n->y))
        RangeTree_find_recur(n->left, dir^1, r, rtrf, user_data, found);
    if(*found)
        return;
    if(RangeTree_inside_rectangle(r, n)){
        if(rtrf(n->x, n->y, n->info, user_data)){
            (*found) = TRUE;
            return;
            }
        }
    if(dir ? (n->x <= r->br_x) : (n->y <= r->br_y))
        RangeTree_find_recur(n->right, dir^1, r, rtrf,
                             user_data, found);
    return;
    }

static void RangeTree_insert(RangeTree *rt, RangeTree_Node *rtn){
    register RangeTree_Node *n, *parent;
    register gboolean dir, dim = FALSE;
    if(!rt->root){
        rt->root = rtn;
        return;
        }
    for(n = parent = rt->root; n; dim ^= 1){
        dir = dim ? (rtn->x < n->x)
                  : (rtn->y < n->y);
        parent = n;
        n = dir ? n->left : n->right;
        }
    if(dir)
        parent->left = rtn;
    else
        parent->right = rtn;
    return;
    }

static void RangeTree_insert_list(RangeTree *rt, GPtrArray *list){
    register gint i, j;
    register RangeTree_Node *rtn, *swap_rtn;
    srand(list->len);
    for(i = 0; i < list->len; i++){ /* Shuffle list */
        rtn = list->pdata[i];
        j = ((gint)lrand48()) % list->len;
        swap_rtn = list->pdata[j];
        list->pdata[j] = rtn;
        list->pdata[i] = swap_rtn;
        }
    for(i = 0; i < list->len; i++){ /* Insert list */
        rtn = list->pdata[i];
        RangeTree_insert(rt, rtn);
        }
    return;
    }
/* The list contains RangeTree_Data objects to be inserted.
 * As the range tree is not balanced, this allows
 * a set of points to be inserted in random order
 * to avoid the worst-case performance of the RangeTree.
 */

static gint RangeTree_insert_recent_collect(gpointer key,
                                            gpointer value,
                                            gpointer data){
    register GPtrArray *list = data;
    g_ptr_array_add(list, value);
    return FALSE;
    }

static void RangeTree_insert_recent(RangeTree *rt){
    register GPtrArray *node_list = g_ptr_array_new();
    g_tree_traverse(rt->recent_data, RangeTree_insert_recent_collect,
                    G_IN_ORDER, node_list);
    if(node_list->len){
        RangeTree_insert_list(rt, node_list);
        g_tree_destroy(rt->recent_data);
        rt->recent_data = g_tree_new(RangeTree_recent_data_compare);
        }
    g_ptr_array_free(node_list, TRUE);
    return;
    }

static gboolean RangeTree_find_internal(RangeTree *rt,
                    gint x_start, gint x_length,
                    gint y_start, gint y_length,
                    RangeTree_ReportFunc rtrf, gpointer user_data){
    gboolean found = FALSE;
    RangeTree_Rectangle r;
    g_assert(x_length >= 0);
    g_assert(y_length >= 0);
    r.tl_x = x_start;
    r.tl_y = y_start;
    r.br_x = x_start + x_length;
    r.br_y = y_start + y_length;
    RangeTree_find_recur(rt->root, FALSE, &r, rtrf, user_data, &found);
    return found;
    }

gboolean RangeTree_find(RangeTree *rt, gint x_start, gint x_length,
                                   gint y_start, gint y_length,
                    RangeTree_ReportFunc rtrf, gpointer user_data){
    RangeTree_insert_recent(rt);
    return RangeTree_find_internal(rt, x_start, x_length,
                                       y_start, y_length,
                                       rtrf, user_data);
    }

static gboolean RangeTree_find_point(gint x, gint y, gpointer info,
                              gpointer user_data){
    return TRUE;
    }

gboolean RangeTree_check_pos(RangeTree *rt,
                             gint x, gint y){
    RangeTree_Node rtn;
    rtn.x = x;
    rtn.y = y;
    if(g_tree_lookup(rt->recent_data, &rtn))
        return TRUE;
    return RangeTree_find_internal(rt, x, 1, y, 1,
                                   RangeTree_find_point, NULL);
    }

gboolean RangeTree_is_empty(RangeTree *rt){
    RangeTree_insert_recent(rt);
    return rt->root?FALSE:TRUE;
    }

static gboolean RangeTree_traverse_recur(RangeTree_Node *rtn,
                                     RangeTree_ReportFunc rtrf,
                                     gpointer user_data){
    if(!rtn)
        return FALSE;
    if(RangeTree_traverse_recur(rtn->left, rtrf, user_data)
    || rtrf(rtn->x, rtn->y, rtn->info, user_data)
    || RangeTree_traverse_recur(rtn->right, rtrf, user_data))
        return TRUE;
    return FALSE;
    }

gboolean RangeTree_traverse(RangeTree *rt,
                        RangeTree_ReportFunc rtrf, gpointer user_data){
    return RangeTree_traverse_recur(rt->root, rtrf, user_data);
    }

