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

static gint compare_func(gpointer data_a, gpointer data_b, gpointer user_data){
    /* g_message("comparing [%d][%d]", (gint)data_a, (gint)data_b); */
    return data_a - data_b;
    }

static gboolean traverse_func(gpointer data, gpointer user_data){
    g_message("traverse [%d]", GPOINTER_TO_INT(data));
    return FALSE;
    }

int main(void){
    register SplayTree_Set *sts = SplayTree_Set_create(compare_func,
                                                       NULL, NULL);
    register SplayTree *st = NULL;
    int i;
    int num[100] = { 87,  18,  95,  86,  32,  96,  25,  21,  36,  65,
                     42,  78,  72,  29,  85,  66,  22,  13,  61,  53,
                     88,  14,  28,  84,  47,  67,  40,  26,  24,   1,
                     70,  60,  15,  97,   9,   6,  63,  45,  37,  64,
                     79,   3,  41,  46,  16,  38,  35,  10,  59,  71,
                     17,  93,   8,  81,  51,  31,  82,  98,  80,  12,
                     27,  62,  30,   5,  90,  52,   2,  68,  89,  34,
                     94,  20,  11,   4,  92,  58,  99,  77,   7,  23,
                     43,  33,   0,  73,  57,  39,  91,  56,  83,  75,
                     74,  49,  50,  44,  55,  19,  54,  76,  69,  48 };
    for(i = 0; i < 100; i++){
        g_message("insert [%d]", num[i]);
        st = SplayTree_insert(st, sts, GINT_TO_POINTER(num[i]));
        }
    for(i = 0; i < 100; i++){
        g_message("lookup [%d]", num[i]);
        st = SplayTree_lookup(st, sts, GINT_TO_POINTER(num[i]));
        g_assert(GPOINTER_TO_INT(st->data) == num[i]);
        }
#if 0
    for(i = 0; i < 100; i++){
        g_message("remove [%d]", num[i]);
        st = SplayTree_remove(st, sts, GINT_TO_POINTER(num[i]), NULL);
        }
#endif /* 0 */
    SplayTree_traverse(st, sts, traverse_func);
    st = SplayTree_splay(st, sts, GINT_TO_POINTER(0));
    for(i = 0; i < 100; i++){
        g_message("removing [%d]", GPOINTER_TO_INT(st->data));
        st = SplayTree_remove(st, sts, st->data, NULL);
        }
    /**/
    g_assert(!st);
    SplayTree_destroy(st, sts);
    SplayTree_Set_destroy(sts);
    g_message("done");
    return 0;
    }

/**/

