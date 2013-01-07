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

#include "recyclebin.h"

int main(void){
    register RecycleBin *rb = RecycleBin_create("test",
                                     sizeof(gpointer), 3);
    register gint i;
    register GPtrArray *list = g_ptr_array_new();
    for(i = 0; i < 12; i++)
        g_ptr_array_add(list, RecycleBin_alloc(rb));
    for(i = 0; i < 12; i++)
        RecycleBin_recycle(rb, list->pdata[i]);
    RecycleBin_destroy(rb);
    g_ptr_array_free(list, TRUE);
    return 0;
    }

