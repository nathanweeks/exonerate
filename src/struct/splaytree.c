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

#include "splaytree.h"

SplayTree_Set *SplayTree_Set_create(SplayTree_CompareFunc compare_func,
                                    SplayTree_FreeFunc free_func,
                                    gpointer user_data){
    register SplayTree_Set *sts = g_new(SplayTree_Set, 1);
    sts->recycle = RecycleBin_create("SplayTree", sizeof(SplayTree), 1024);
    g_assert(compare_func);
    sts->compare_func = compare_func;
    sts->free_func = free_func;
    sts->user_data = user_data;
    return sts;
    }

void SplayTree_Set_destroy(SplayTree_Set *sts){
    RecycleBin_destroy(sts->recycle);
    g_free(sts);
    return;
    }

gsize SplayTree_Set_memory_usage(SplayTree_Set *sts){
    return RecycleBin_memory_usage(sts->recycle);
    }

SplayTree *SplayTree_splay(SplayTree *st, SplayTree_Set *sts, gpointer data){
    register SplayTree *l, *r, *t;
    register gint comp;
    SplayTree n;
    if(!st)
        return st;
    n.left = n.right = NULL;
    l = r = &n;
    while((comp = sts->compare_func(data, st->data, sts->user_data))){
        if(comp < 0){
            if(!st->left)
                break;
            if(sts->compare_func(data, st->left->data, sts->user_data) < 0){
                /* Rotate right */
                t = st->left;
                st->left = t->right;
                t->right = st;
                st = t;
                if(!st->left)
                    break;
                }
            r->left = st; /* Link right */
            r = st;
            st = st->left;
        } else {
            if(!st->right)
                break;
            if(sts->compare_func(data, st->right->data, sts->user_data) > 0){
                /* Rotate left */
                t = st->right;
                st->right = t->left;
                t->left = st;
                st = t;
                if(!st->right)
                    break;
                }
            l->right = st; /* Link left */
            l = st;
            st = st->right;
            }
        }
    l->right = st->left;  /* Assemble */
    r->left = st->right;
    st->left = n.right;
    st->right = n.left;
    return st;
    }

SplayTree *SplayTree_insert(SplayTree *st, SplayTree_Set *sts, gpointer data){
    register SplayTree *n;
    register gint comp;
    if(!st){
        n = RecycleBin_alloc(sts->recycle);
        n->data = data;
        n->left = n->right = NULL;
        return n;
        }
    st = SplayTree_splay(st, sts, data);
    comp = sts->compare_func(data, st->data, sts->user_data);
    if(!comp) /* Already in tree */
        return st;
    n = RecycleBin_alloc(sts->recycle);
    n->data = data;
    if(comp < 0){
        n->left = st->left;
        n->right = st;
        st->left = NULL;
    } else {
        n->right = st->right;
        n->left = st;
        st->right = NULL;
        }
    return n;
    }
/* If the data is already present, the new data is not used.
 */

SplayTree *SplayTree_remove(SplayTree *st, SplayTree_Set *sts,
                            gpointer data, gpointer *result){
    register SplayTree *n;
    register gint comp;
    if(!st)
        return NULL;
    st = SplayTree_splay(st, sts, data);
    comp = sts->compare_func(data, st->data, sts->user_data);
    if(!comp){ /* found */
        if(st->right){
            n = SplayTree_splay(st->right, sts, data);
            g_assert(!n->left);
            n->left = st->left;
        } else {
            n = st->left;
            }
        if(result)
            (*result) = st->data;
        RecycleBin_recycle(sts->recycle, st);
        return n;
        }
    return st; /* not found */
    }
/* If the data is found, (*result) is set to the removed data */

SplayTree *SplayTree_lookup(SplayTree *st, SplayTree_Set *sts, gpointer data){
    return st?SplayTree_splay(st, sts, data):NULL;
    }

void SplayTree_destroy(SplayTree *st, SplayTree_Set *sts){
    if(!st)
        return;
    SplayTree_destroy(st->left, sts);
    SplayTree_destroy(st->right, sts);
    if(sts->free_func)
        sts->free_func(st->data, sts->user_data);
    RecycleBin_recycle(sts->recycle, st);
    return;
    }

static gboolean SplayTree_traverse_recur(SplayTree *st,
                                         SplayTree_Traverse_Func traverse_func,
                                         gpointer user_data){
    if(!st)
        return FALSE;
    if(SplayTree_traverse_recur(st->left, traverse_func, user_data)
    || traverse_func(st->data, user_data)
    || SplayTree_traverse_recur(st->right, traverse_func, user_data))
        return TRUE;
    return FALSE;
    }

void SplayTree_traverse(SplayTree *st, SplayTree_Set *sts,
                        SplayTree_Traverse_Func traverse_func){
    SplayTree_traverse_recur(st, traverse_func, sts->user_data);
    return;
    }

/**/

