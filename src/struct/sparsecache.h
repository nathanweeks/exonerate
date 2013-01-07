/****************************************************************\
*                                                                *
*  Sparse Cache Object                                           *
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

#ifndef INCLUDED_SPARSECACHE_H
#define INCLUDED_SPARSECACHE_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <glib.h>

/* Currently using simple demand paging
 *
 * FIXME: optimisation: add option to flush of LRU pages ?
 */

typedef gpointer (*SparseCache_GetFunc)(gint pos, gpointer page_data,
                                        gpointer user_data);
typedef gint (*SparseCache_CopyFunc)(gint start, gint length, gchar *dst,
                                     gpointer page_data, gpointer user_data);

typedef struct {
                   gpointer data;
        SparseCache_GetFunc get_func;
       SparseCache_CopyFunc copy_func;
                      gsize data_size;
} SparseCache_Page;

typedef SparseCache_Page *(*SparseCache_FillFunc)(gint start,
                                                  gpointer user_data);

typedef void (*SparseCache_EmptyFunc)(SparseCache_Page *page, gpointer user_data);
typedef void (*SparseCache_FreeFunc)(gpointer user_data);

#define SparseCache_PAGE_SIZE_BIT_WIDTH 12
#define SparseCache_PAGE_SIZE (1 << SparseCache_PAGE_SIZE_BIT_WIDTH)
#define SparseCache_pos2page(pos) ((pos) >> SparseCache_PAGE_SIZE_BIT_WIDTH)

typedef struct {
                    gint   ref_count;
        SparseCache_Page **page_list;
                    gint   page_total;
                    gint   page_used;
                    gint   length;
                   gsize   page_memory_usage;
    SparseCache_FillFunc   fill_func;
   SparseCache_EmptyFunc   empty_func;
    SparseCache_FreeFunc   free_func;
                gpointer   user_data;
} SparseCache;

SparseCache *SparseCache_create(gint length,
                                SparseCache_FillFunc fill_func,
                                SparseCache_EmptyFunc empty_func,
                                SparseCache_FreeFunc free_func,
                                gpointer user_data);
        void SparseCache_destroy(SparseCache *sc);
SparseCache *SparseCache_share(SparseCache *sc);
    gpointer SparseCache_get(SparseCache *sc, gint pos);
       gsize SparseCache_memory_usage(SparseCache *sc);
        void SparseCache_copy(SparseCache *sc, gint start, gint length,
                              gchar *dst);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* INCLUDED_SPARSECACHE_H */

