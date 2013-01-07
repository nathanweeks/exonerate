/****************************************************************\
*                                                                *
*  Efficient single-linked list routines.                        *
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

#include <string.h> /* For strcmp() */
#include "slist.h"

typedef struct {
    gchar *crib[9];
     gint  count;
} SList_TestInfo;

static gboolean slist_test_traverse(gpointer data, gpointer user_data){
    register gchar *word = data;
    SList_TestInfo *sti = user_data;
    g_message("Word [%5s] crib[%5s] count[%d]",
              word, sti->crib[sti->count], sti->count);
    g_assert(!strcmp(word, sti->crib[sti->count]));
    sti->count++;
    return FALSE;
    }

int main(void){
    register SListSet *slist_set = SListSet_create();
    register SList *a = SList_create(slist_set),
                   *b = SList_create(slist_set), *c;
    SList_TestInfo sti = { {"The", "quick", "brown", "fox",
                            "jumps", "over", "a", "lazy", "dog"},
                          0};
    g_message("check [%d][%d]", SList_isempty(a), SList_isempty(b));
    SList_queue(a, "jumps");
    SList_stack(a, "fox");
    SList_queue(a, "over");
    SList_stack(a, "brown");
    SList_queue(a, "a");
    SList_stack(b, "quick");
    SList_queue(a, "lazy");
    SList_stack(b, "The");
    SList_queue(a, "dog");
    c = SList_join(b, a);
    SList_destroy(a);
    SList_traverse(c, slist_test_traverse, &sti);
    SList_destroy(c);
    SListSet_destroy(slist_set);
    return 0;
    }

