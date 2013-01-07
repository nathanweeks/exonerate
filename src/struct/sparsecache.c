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

#include "sparsecache.h"

/**/

static void SparseCache_Page_destroy(SparseCache_Page *page){
    if(page->data)
        g_free(page->data);
    g_free(page);
    return;
    }

SparseCache *SparseCache_create(gint length,
                                SparseCache_FillFunc fill_func,
                                SparseCache_EmptyFunc empty_func,
                                SparseCache_FreeFunc free_func,
                                gpointer user_data){
    register SparseCache *sc = g_new(SparseCache, 1);
    g_assert(fill_func);
    sc->ref_count = 1;
    sc->page_total = (length >> SparseCache_PAGE_SIZE_BIT_WIDTH) + 1;
    sc->page_used = 0;
    sc->page_list = g_new0(SparseCache_Page*, sc->page_total);
    sc->page_memory_usage = 0;
    sc->length = length;
    sc->fill_func = fill_func;
    sc->empty_func = empty_func;
    sc->free_func = free_func;
    sc->user_data = user_data;
    return sc;
    }

void SparseCache_destroy(SparseCache *sc){
    register SparseCache_Page *page;
    register gint i;
    if(--sc->ref_count)
        return;
    for(i = 0; i < sc->page_total; i++){
        page = sc->page_list[i];
        if(!page)
            continue;
        sc->page_memory_usage -= (sizeof(SparseCache_Page)+page->data_size);
        if(sc->empty_func)
           sc->empty_func(page, sc->user_data);
        else
            SparseCache_Page_destroy(page);
        }
    g_free(sc->page_list);
    if(sc->free_func && sc->user_data)
        sc->free_func(sc->user_data);
    g_free(sc);
    return;
    }

SparseCache *SparseCache_share(SparseCache *sc){
    sc->ref_count++;
    return sc;
    }

static SparseCache_Page *SparseCache_fill(SparseCache *sc, gint page_id){
    register SparseCache_Page *page;
    g_assert(!sc->page_list[page_id]);
    page = sc->fill_func((page_id << SparseCache_PAGE_SIZE_BIT_WIDTH),
                         sc->user_data);
    sc->page_list[page_id] = page;
    sc->page_used++;
    sc->page_memory_usage += (sizeof(SparseCache_Page)+page->data_size);
    return page;
    }

gpointer SparseCache_get(SparseCache *sc, gint pos){
    register gint page_id = SparseCache_pos2page(pos);
    register SparseCache_Page *page = sc->page_list[page_id];
    g_assert(pos >= 0);
    g_assert(pos < sc->length);
    if(!page)
        page = SparseCache_fill(sc, page_id);
    return page->get_func((pos & (SparseCache_PAGE_SIZE-1)),
                          page->data, sc->user_data);
    }

gsize SparseCache_memory_usage(SparseCache *sc){
    return sizeof(SparseCache)
         + (sizeof(SparseCache_Page*)*sc->page_total)
         + sc->page_memory_usage;
    }

void SparseCache_copy(SparseCache *sc, gint start, gint length,
                      gchar *dst){
    register gint pos = start, dst_pos = 0, page_id, copy_len,
                  read_remain = length, page_remain, page_start;
    register SparseCache_Page *page;
    do {
        page_id = SparseCache_pos2page(pos);
        page = sc->page_list[page_id];
        if(!page)
            page = SparseCache_fill(sc, page_id);
        page_start = page_id*SparseCache_PAGE_SIZE;
        page_remain = page_start+SparseCache_PAGE_SIZE-pos;
        copy_len = MIN(page_remain, read_remain);
        g_assert(page->copy_func);
        dst_pos += page->copy_func(pos-page_start, copy_len, dst+dst_pos,
                                   page->data, sc->user_data);
        read_remain -= copy_len;
        pos += copy_len;
    } while(read_remain > 0);
    return;
    }

/**/

