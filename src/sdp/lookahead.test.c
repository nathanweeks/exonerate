/****************************************************************\
*                                                                *
*  C4 dynamic programming library - SDP Lookahead Object         *
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

#include "lookahead.h"

static void test_free_func(gpointer data, gpointer user_data){
    register gchar *word = data;
    g_assert(word);
    g_message("freeing [%s]", word);
    g_free(word);
    return;
    }

int main(void){
    register Lookahead *lookahead = Lookahead_create(0, 32,
                                  test_free_func, NULL);
    register gchar *word;
    word = Lookahead_get(lookahead, 0);
    g_assert(!word);
    /**/
    Lookahead_set(lookahead, 0, g_strdup("zero"));
    word = Lookahead_get(lookahead, 0);
    g_message("Have zero [%s]", word);
    /**/
    Lookahead_set(lookahead, 7, g_strdup("seven"));
    word = Lookahead_get(lookahead, 7);
    g_message("Have seven [%s]", word);
    /**/
    Lookahead_set(lookahead, 12, g_strdup("twelve"));
    word = Lookahead_get(lookahead, 12);
    g_message("Have twelve [%s]", word);
    /**/
    Lookahead_set(lookahead, 24, g_strdup("twentyfour"));
    word = Lookahead_get(lookahead, 24);
    g_message("Have twentyfour [%s]", word);
    /**/
    Lookahead_next(lookahead);
    word = Lookahead_get(lookahead, 0);
    g_message("Moved to [%s]", word);
    /**/
    Lookahead_next(lookahead);
    word = Lookahead_get(lookahead, 0);
    g_message("Moved to [%s]", word);
    /**/
    Lookahead_get(lookahead, 12);
    word = Lookahead_get(lookahead, 12);
    g_message("Should have twentyfour [%s]", word);
    /**/
    Lookahead_move(lookahead, 20);
    word = Lookahead_get(lookahead, 4);
    g_message("Should still have twentyfour [%s]", word);
    /**/
    Lookahead_set(lookahead, 20, g_strdup("forty"));
    Lookahead_move(lookahead, 24);
    Lookahead_next(lookahead);
    word = Lookahead_get(lookahead, 0);
    g_message("Should now have forty [%s]", word);
    Lookahead_destroy(lookahead);
    return 0;
    }

