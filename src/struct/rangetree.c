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

#include <search.h> /* For tdelete(), tfind(), tsearch(), twalk() */
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
    rt->recent_data = NULL;
    rt->node_recycle = RecycleBin_create("RangeTree_Node",
                                         sizeof(RangeTree_Node), 256);
    return rt;
    }

void RangeTree_destroy(RangeTree *rt, RangeTree_FreeFunc rtff,
                              gpointer user_data){
    while (rt->recent_data)
        tdelete((void *)*(RangeTree_Node **)rt->recent_data, &rt->recent_data,
                RangeTree_recent_data_compare);
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
    g_assert(!tfind((void *)rtn, &rt->recent_data,
           RangeTree_recent_data_compare));
    tsearch((void *)rtn, &rt->recent_data, RangeTree_recent_data_compare);
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


static void RangeTree_insert_recent(RangeTree *rt){
    /* Remove root node with each iteration */
    while (rt->recent_data) {
        RangeTree_insert(rt, *((RangeTree_Node **)(rt->recent_data)));
        tdelete(*(RangeTree_Node **)rt->recent_data, &rt->recent_data, 
                RangeTree_recent_data_compare);
    }
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
    if(tfind((void *)&rtn, &rt->recent_data, RangeTree_recent_data_compare))
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

