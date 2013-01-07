/****************************************************************\
*                                                                *
*  Deja-vu library for fast linear space repeat finding          *
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

#include <string.h> /* For strlen() */

#include "dejavu.h"

static void test_display_repeat(gint first_pos, gint curr_pos,
                                gint length, gchar *seq, gint len,
                                gpointer user_data){
    register gint i;
    register gint *total = user_data;
    g_print(">repeat_%d (%d, %d, %d)\n[",
            (*total)++, first_pos, curr_pos, length);
    for(i = 0; i < length; i++)
        g_print("%c", seq[first_pos+i]);
    g_print("]\n");
    return;
    }

int main(void){
    register gchar *seq = "TAGACTGTAAGTCCTCTGTGTAAGCTCGTTA";
    register gint min_wordlen = 1, max_wordlen = 7;
    register DejaVu *dv = DejaVu_create(seq, strlen(seq));
    gint total = 0;
    DejaVu_traverse(dv, min_wordlen, max_wordlen,
                    test_display_repeat, &total, NULL, 1);
    g_message("total = %d", total);
    g_assert(total == 82);
    DejaVu_destroy(dv);
    return 0;
    }

