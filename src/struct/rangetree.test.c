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

#include "rangetree.h"

/*
   9     1
   8         1
   7     +-------+
   6     |     1 |
   5     1   2   |   1
   4     |       |
   3     +-------1
   2
   1   1
   0 1 2 3 4 5 6 7 8 9
*/

static gboolean test_rtrf(gint x, gint y, gpointer info,
                          gpointer user_data){
    register gint *total = user_data;
    (*total)++;
    g_message("Found (%d,%d) (count:%d)", x, y, *total);
    return FALSE;
    }

int main(void){
    register RangeTree *rt = RangeTree_create();
    gint total = 0;
    RangeTree_add(rt, 3, 5, NULL); /* within */
    RangeTree_add(rt, 5, 8, NULL);
    RangeTree_add(rt, 5, 5, NULL); /* within */
    RangeTree_add(rt, 2, 1, NULL);
    RangeTree_add(rt, 9, 5, NULL);
    RangeTree_add(rt, 7, 3, NULL); /* within */
    RangeTree_add(rt, 3, 9, NULL);
    RangeTree_add(rt, 6, 6, NULL); /* within */
    /**/
    RangeTree_find(rt, 3, 5, 3, 5, test_rtrf, &total);
    g_message("total is [%d] (expect 4)", total);
    g_assert(total == 4);
    /**/
    total = 0;
    RangeTree_find(rt, 5, 1, 5, 1, test_rtrf, &total);
    g_message("total is [%d] (expect 1)", total);
    g_assert(total == 1);
    /**/
    RangeTree_destroy(rt, NULL, NULL);
    return 0;
    }

