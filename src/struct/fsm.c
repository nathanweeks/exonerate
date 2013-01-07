/****************************************************************\
*                                                                *
*  Library for FSM-based word matching.                          *
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

/* Simplified FSM code derived from
 * my earlier implementation in the ESTate package:
 * http://www.hgmp.mrc.ac.uk/~gslater/ESTate/
 */

#include "fsm.h"

#include <string.h> /* For memcpy() */

FSM *FSM_create(gchar *alphabet, FSM_Join_Func merge_func,
                FSM_Join_Func combine_func, gpointer user_data){
    register FSM *f = g_new0(FSM, 1);
    register guchar *p;
    g_assert(alphabet);
    g_assert(merge_func);
    g_assert(combine_func);
    for(p = ((guchar*)alphabet); *p; p++)
        f->index[*p] = 1;
    for(p = f->index+ALPHABETSIZE; p > f->index; p--)
        if(*p)
            *p = ++f->width;
    f->width++; /* Add one for NULL */
    memcpy(f->insertion_filter, f->index, sizeof(gchar)*ALPHABETSIZE);
    memcpy(f->traversal_filter, f->index, sizeof(gchar)*ALPHABETSIZE);
    f->merge_func = merge_func;
    f->combine_func = combine_func;
    f->user_data = user_data;
    f->recycle = RecycleBin_create("FSM", sizeof(FSM_Node)*f->width, 1024);
    f->chunk_count = 1;
    f->root = RecycleBin_alloc_blank(f->recycle);
    f->is_compiled = FALSE;
    return f;
    }

gsize FSM_memory_usage(FSM *f){
    return sizeof(FSM)
         + (sizeof(FSM_Node)
            * f->width
            * f->chunk_count);
    }

void FSM_info(FSM *f){
    g_print("FSM information:\n"
            " Alphabet size   = %d\n"
            " Size of node    = %d\n"
            " Size of chunk   = %d\n"
            " Chunks used     = %lu\n"
            " Bytes used      = %lu\n",
            f->width,
            (gint)sizeof(FSM_Node),
            (gint)sizeof(FSM_Node)*f->width,
            f->chunk_count,
            (gulong)FSM_memory_usage(f));
    return;
    }

void FSM_destroy(FSM *f){
    RecycleBin_destroy(f->recycle);
    g_free(f);
    return;
    }

static void FSM_destroy_with_data_recur(FSM *f, FSM_Node *n,
                                        FSM_Destroy_Func fdf,
                                        gpointer user_data){
    register gint i;
    register FSM_Node *t;
    for(i = 1; i < f->width; i++){
        if(n[i].data){
            fdf(n[i].data, user_data);
            n[i].data = NULL;
            }
        if(n[i].next){
            t = n[i].next;
            n[i].next = NULL;
            FSM_destroy_with_data_recur(f, t, fdf, user_data);
            }
        }
    return;
    }

void FSM_destroy_with_data(FSM *f, FSM_Destroy_Func fdf,
                           gpointer user_data){
    FSM_destroy_with_data_recur(f, f->root, fdf, user_data);
    FSM_destroy(f);
    return;
    }

gpointer FSM_add(FSM *f, gchar *seq, guint len, gpointer node_data){
    register FSM_Node **ptp, *n = f->root;
    register guchar *useq, *end;
    g_assert(!f->is_compiled);
    for(useq = (guchar*)seq, end = useq+len-1; useq < end; useq++){
        if(!f->insertion_filter[*useq])
            g_error("Illegal char (%d)[%c] in FSM word [%s]",
                    *useq, *useq, seq);
        if(!*(ptp = &n[f->insertion_filter[*useq]].next)){
            *ptp = RecycleBin_alloc_blank(f->recycle);
            f->chunk_count++;
            }
        n = *ptp;
        }
    if(n[f->insertion_filter[*useq]].data){
        if(node_data)
            n[f->insertion_filter[*useq]].data =
                f->merge_func(n[f->insertion_filter[*useq]].data,
                              node_data, f->user_data);
    } else {
        n[f->insertion_filter[*useq]].data = node_data;
        }
    return n[f->insertion_filter[*useq]].data;
    }

/* FSM_compile : Converts the pre-FSM trie to the FSM.
   1: All position zero nodes must have no score and point to root.
   2: All unused nodes must me made to point to the node which
      would be reached by insertion of their longest suffix
      (this will always be above the current node, and it's location
      is independent of any visits to root on this path).
   3: The scores from subsequences must also be inherited.
   4: This is achieved with a level-order trie traversal.
   5: A queue is used to remove traversal recursion and backtracking.
   6: No space is required for the queue, as it is makes temporary
      use of the position zero nodes during compilation.
   7: This algorithm is linear with the number of states.
*/
void FSM_compile(FSM *f){
    register gint i;
    FSM_Node a;
    register FSM_Node *suffix, *prev;
    register FSM_Node *out = &a, *in = out;
    g_assert(!f->is_compiled);
    out->next = in;              /* Initialise queue */
    f->root->next = f->root;
    for(i = 1; i < f->width; i++)
        if(f->root[i].next){
            in->data = f->root;
            in->next = f->root[i].next;
            in       = f->root[i].next;
        } else
            f->root[i].next   = f->root;
    while(out != in){
        prev   = out;
        suffix = out->data;
        out    = out->next;
        for(i = 1; i < f->width; i++){
            if(suffix[i].data){
                if(out[i].data)
                    out[i].data = f->combine_func(out[i].data,
                                                  suffix[i].data,
                                                  f->user_data);
                else
                    out[i].data = suffix[i].data;
                }
            if(out[i].next){
                in->data = suffix[i].next;
                in->next = out[i].next;
                in       = out[i].next;
            } else
                out[i].next = suffix[i].next;
            }
        prev->next = f->root;
        prev->data = NULL;
        }
    out->next = f->root;
    out->data = NULL;
    f->is_compiled = TRUE;
    return;
    }

void FSM_traverse(FSM *f, gchar *seq, FSM_Traverse_Func ftf,
                  gpointer user_data){
    register FSM_Node *n = f->root;
    register guchar *p = (guchar*)seq;
    register gint c;
    g_assert(f->is_compiled);
    do {
        if(n[c = f->traversal_filter[*p]].data)
           ftf(p-(guchar*)seq, n[c].data, user_data);
        n = n[c].next;
    } while(*++p);
    return;
    }

static void FSM_add_filter(FSM *f, guchar *src, guchar *dst){
    register gint i;
    if(dst){ /* Merge filter */
        for(i = 0; i < ALPHABETSIZE; i++)
            dst[i] = dst[src[i]];
    } else { /* Reset to original index */
        memcpy(dst, f->index, sizeof(gchar)*ALPHABETSIZE);
        }
    return;
    }

void FSM_add_insertion_filter(FSM *f, guchar *filter){
    FSM_add_filter(f, filter, f->insertion_filter);
    return;
    }

void FSM_add_traversal_filter(FSM *f, guchar *filter){
    FSM_add_filter(f, filter, f->traversal_filter);
    return;
    }

/**/

