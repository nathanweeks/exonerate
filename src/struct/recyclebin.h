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

#ifndef INCLUDED_RECYCLEBIN_H
#define INCLUDED_RECYCLEBIN_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <glib.h>

typedef struct RecycleBin_Node {
     struct RecycleBin_Node *next;
} RecycleBin_Node;

typedef struct {
                gint  ref_count;
               gchar *name;
           GPtrArray *chunk_list;
                gint  chunk_pos;
                gint  nodes_per_chunk;
               gsize  node_size;
                gint  count;
     RecycleBin_Node *recycle;
} RecycleBin;

RecycleBin *RecycleBin_create(gchar *name, gsize node_size,
                              gint nodes_per_chunk);
      void  RecycleBin_destroy(RecycleBin *recycle_bin);
RecycleBin *RecycleBin_share(RecycleBin *recycle_bin);

   gsize RecycleBin_memory_usage(RecycleBin *recycle_bin);

gpointer RecycleBin_alloc(RecycleBin *recycle_bin);
gpointer RecycleBin_alloc_blank(RecycleBin *recycle_bin);
    void RecycleBin_recycle(RecycleBin *recycle_bin, gpointer data);
    void RecycleBin_profile(void);
#define  RecycleBin_total(rb) ((rb)->count)

/* FIXME: optimisation: use macros */

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* INCLUDED_RECYCLEBIN_H */

