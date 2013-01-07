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

static gpointer test_get_func(gint pos, gpointer page_data,
                          gpointer user_data){
    register gint *data = page_data;
    return GINT_TO_POINTER(data[pos]);
    }

static SparseCache_Page *test_fill_func(gint start, gpointer user_data){
    register SparseCache_Page *page = g_new(SparseCache_Page, 1);
    register gint i, *data = g_new(gint, SparseCache_PAGE_SIZE);
    g_message("Filling from [%d] with [%d]", start, SparseCache_PAGE_SIZE);
    for(i = 0; i < SparseCache_PAGE_SIZE; i++){
        data[i] = start + i;
        g_message(" ... filling [%d]", data[i]);
        }
    page->data = (gpointer)data;
    page->data_size = sizeof(gint)*SparseCache_PAGE_SIZE;
    page->get_func = test_get_func;
    g_message("done filling");
    return page;
    }

int main(void){
    register SparseCache *sc = SparseCache_create(1000, test_fill_func,
                                                  NULL, NULL, NULL);
    register gint i;
    for(i = 0; i < 10; i++){
        g_message("checking [%d]", i);
        g_assert(GPOINTER_TO_INT(SparseCache_get(sc, i)) == i);
        }
    for(i = 20; i >= 0; i--){
        g_message("checking [%d]", i);
        g_assert(GPOINTER_TO_INT(SparseCache_get(sc, i)) == i);
        }
    for(i = 30; i >= 10; i--){
        g_message("checking [%d]", i);
        g_assert(GPOINTER_TO_INT(SparseCache_get(sc, i)) == i);
        }
    for(i = 10; i < 50; i++){
        g_message("checking [%d]", i);
        g_assert(GPOINTER_TO_INT(SparseCache_get(sc, i)) == i);
        }
    for(i = 60; i >= 0; i--){
        g_message("checking [%d]", i);
        g_assert(GPOINTER_TO_INT(SparseCache_get(sc, i)) == i);
        }
    for(i = 0; i < 60; i++){
        g_message("checking [%d]", i);
        g_assert(GPOINTER_TO_INT(SparseCache_get(sc, i)) == i);
        }
    for(i = 0; i < 100; i++){
        g_message("checking [%d]", i);
        g_assert(GPOINTER_TO_INT(SparseCache_get(sc, i)) == i);
        }
    SparseCache_destroy(sc);
    g_message("SparseCache OK");
    return 0;
    }

