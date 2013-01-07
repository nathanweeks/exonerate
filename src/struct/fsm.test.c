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

#include <string.h> /* For strlen() */
#include "fsm.h"

static gpointer test_merge_Func(gpointer a, gpointer b,
                                gpointer user_data){
    register gint ia = GPOINTER_TO_INT(a),
                  ib = GPOINTER_TO_INT(b);
    g_assert(a);
    g_assert(b);
    return GINT_TO_POINTER(ia + ib);
    }

static void test_traverse_Func(guint seq_pos,
                               gpointer node_data,
                               gpointer user_data){
    register gint node_count = GPOINTER_TO_INT(node_data),
                 *total      = (gint*)user_data;
    *total += node_count;
    return;
    }

int main(void){
    gint count = 0;
    register FSM *f = FSM_create("abcdefghijklmnopqrstuvwxyz",
                      test_merge_Func, test_merge_Func, NULL);
#define TESTSETSIZE 32
    char *testset[TESTSETSIZE] = { /* C KEYWORDS */
        "auto",     "break",    "case",     "char",
        "const",    "continue", "default",  "do",
        "double",   "else",     "enum",     "extern",
        "float",    "for",      "goto",     "if",
        "int",      "long",     "register", "return",
        "short",    "signed",   "sizeof",   "static",
        "struct",   "switch",   "typedef",  "union",
        "unsigned", "void",     "volatile", "while"  };
    char *testsen = "-unSIGNed--switCHAR---DOUBLElse--sizeoFOR";
/*        Total=9            2       1 1    1   1  1       1 1 */
    register gint i;
    register guchar *tolower_filter = (guchar*)
    "----------------------------------------------------------------"
    "-abcdefghijklmnopqrstuvwxyz------abcdefghijklmnopqrstuvwxyz-----"
    "----------------------------------------------------------------"
    "----------------------------------------------------------------";
    for(i = 0; i < TESTSETSIZE; i++)
        FSM_add(f, testset[i], strlen(testset[i]), GINT_TO_POINTER(1));
    FSM_compile(f);
    FSM_add_traversal_filter(f, tolower_filter);
    FSM_traverse(f, testsen, test_traverse_Func, &count);
    g_print("FSM count is [%d]\n", count);
    g_assert(count == 9);
    FSM_destroy(f);
    return 0;
    }

