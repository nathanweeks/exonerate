/****************************************************************\
*                                                                *
*  Efficient Memory Allocation Routines                          *
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

#include <search.h> /* For tdelete(), tfind(), tsearch(), and twalk() */
#include <string.h>
#include "recyclebin.h"

static void *global_recycle_bin_tree = NULL;

static int RecycleBin_compare(const void *a,
                              const void *b){
    if ((char*)a - (char*)b == 0)
        return 0;
    else if ((char*)a - (char*)b < 0)
        return -1;
    else
        return 1;
    }

RecycleBin *RecycleBin_create(gchar *name, gsize node_size,
                              gint nodes_per_chunk){
    register RecycleBin *recycle_bin = g_new(RecycleBin, 1);
    g_assert(name);
    g_assert(node_size >= sizeof(gpointer));
    recycle_bin->ref_count = 1;
    recycle_bin->name = g_strdup(name);
    recycle_bin->chunk_list = g_ptr_array_new();
    recycle_bin->chunk_pos = nodes_per_chunk;
    recycle_bin->nodes_per_chunk = nodes_per_chunk;
    recycle_bin->node_size = node_size;
    recycle_bin->count = 0;
    recycle_bin->recycle = NULL;
    tsearch((void *)recycle_bin, &global_recycle_bin_tree, RecycleBin_compare);
    return recycle_bin;
    }

void RecycleBin_destroy(RecycleBin *recycle_bin){
    register gint i;
    if(--recycle_bin->ref_count)
        return;
    g_assert(global_recycle_bin_tree);
    g_assert(tfind((void *)recycle_bin, &global_recycle_bin_tree, 
           RecycleBin_compare));
    tdelete((void *)recycle_bin, &global_recycle_bin_tree, RecycleBin_compare);
    for(i = 0; i < recycle_bin->chunk_list->len; i++)
        g_free(recycle_bin->chunk_list->pdata[i]);
    g_ptr_array_free(recycle_bin->chunk_list, TRUE);
    g_free(recycle_bin->name);
    g_free(recycle_bin);
    return;
    }

RecycleBin *RecycleBin_share(RecycleBin *recycle_bin){
    g_assert(recycle_bin);
    recycle_bin->ref_count++;
    return recycle_bin;
    }

gsize RecycleBin_memory_usage(RecycleBin *recycle_bin){
    g_assert(recycle_bin);
    return recycle_bin->node_size
         * recycle_bin->nodes_per_chunk
         * recycle_bin->chunk_list->len;
    }

gpointer RecycleBin_alloc(RecycleBin *recycle_bin){
    register RecycleBin_Node *node;
    register gchar *chunk;
    g_assert(recycle_bin);
    if(recycle_bin->recycle){
        node = (gpointer)recycle_bin->recycle;
        recycle_bin->recycle = node->next;
    } else {
        if(recycle_bin->chunk_pos == recycle_bin->nodes_per_chunk){
            chunk = g_malloc(recycle_bin->nodes_per_chunk
                           * recycle_bin->node_size);
            g_ptr_array_add(recycle_bin->chunk_list, chunk);
            recycle_bin->chunk_pos = 1;
            node = (RecycleBin_Node*)chunk;
        } else {
            chunk = recycle_bin->chunk_list->pdata
                   [recycle_bin->chunk_list->len-1];
            node = (RecycleBin_Node*) (chunk
                    + (recycle_bin->node_size
                     * recycle_bin->chunk_pos++));
            }
        }
    recycle_bin->count++;
    return (gpointer)node;
    }

gpointer RecycleBin_alloc_blank(RecycleBin *recycle_bin){
    register gpointer data = RecycleBin_alloc(recycle_bin);
    memset(data, 0, recycle_bin->node_size);
    return data;
    }

void RecycleBin_recycle(RecycleBin *recycle_bin, gpointer data){
    register RecycleBin_Node *node = data;
    g_assert(node != recycle_bin->recycle);
    node->next = recycle_bin->recycle;
    recycle_bin->recycle = node;
    recycle_bin->count--;
    return;
    }

static void RecycleBin_profile_traverse(const void *ptr,
                                             VISIT order,
                                               int level){
    if (order == postorder || order == leaf) {
        RecycleBin *recycle_bin = *(RecycleBin **)ptr;
        g_assert(recycle_bin);
        g_message("  RecycleBin [%s] %d Mb (%d items)",
                  recycle_bin->name,
                  (gint)(RecycleBin_memory_usage(recycle_bin)>>20),
                  recycle_bin->count);
        }
    }

void RecycleBin_profile(void){
    g_message("BEGIN RecycleBin profile");
    if(global_recycle_bin_tree)
        twalk(global_recycle_bin_tree, RecycleBin_profile_traverse);
    else
        g_message("no active RecycleBins");
    g_message("END RecycleBin profile");
    return;
    }

/* FIXME: could use Trash Stacks when migrating to glib-2 */

