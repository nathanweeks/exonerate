/****************************************************************\
*                                                                *
*  edit distance model                                           *
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

#include "edit_distance.h"

EditDistance_Data *EditDistance_Data_create(gchar *query,
                                            gchar *target){
    register EditDistance_Data *edd = g_new(EditDistance_Data, 1);
    edd->query = g_strdup(query);
    edd->target = g_strdup(target);
    edd->query_len = strlen(query);
    edd->target_len = strlen(target);
    return edd;
    }

void EditDistance_Data_destroy(EditDistance_Data *edd){
    g_free(edd->query);
    g_free(edd->target);
    g_free(edd);
    return;
    }

static C4_Score edit_match_calc_func(gint query_pos, gint target_pos,
                                     gpointer user_data){
    register EditDistance_Data *edd = (EditDistance_Data*)user_data;
    g_assert(query_pos >= 0);
    g_assert(target_pos >= 0);
    g_assert(query_pos < edd->query_len);
    g_assert(target_pos < edd->target_len);
    if(edd->query[query_pos] == edd->target[target_pos])
         return 0;
    else
         return -1;
    }

C4_Model *EditDistance_create(void){
    register C4_Model *edit_distance = C4_Model_create("edit distance");
    register C4_State
        *main_state   = C4_Model_add_state(edit_distance, "main");
    register C4_Calc
     *indel_calc = C4_Model_add_calc(edit_distance, "indel",  -1,
         NULL, NULL, NULL, NULL, NULL, NULL, C4_Protect_NONE),
     *match_calc = C4_Model_add_calc(edit_distance, "match", 0,
         edit_match_calc_func, NULL, NULL, NULL, NULL, NULL,
         C4_Protect_NONE);
    /**/
    C4_Model_configure_start_state(edit_distance,
                                   C4_Scope_CORNER, NULL, NULL);
    C4_Model_configure_end_state(edit_distance, C4_Scope_CORNER,
                                 NULL, NULL);
    /**/
    C4_Model_add_transition(edit_distance, "start to main",
                     NULL, main_state, 0, 0, NULL, C4_Label_NONE, NULL);
    C4_Model_add_transition(edit_distance, "main to end",
                     main_state, NULL, 0, 0, NULL, C4_Label_NONE, NULL);
    C4_Model_add_transition(edit_distance, "match",
                     main_state, main_state, 1, 1, match_calc,
                     C4_Label_MATCH, NULL);
    C4_Model_add_transition(edit_distance, "query insert",
                     main_state, main_state, 1, 0, indel_calc,
                     C4_Label_GAP, NULL);
    C4_Model_add_transition(edit_distance, "target insert",
                     main_state, main_state, 0, 1, indel_calc,
                     C4_Label_GAP, NULL);
    /**/
    C4_Model_add_portal(edit_distance, "match portal", match_calc,
                        1, 1);
    /**/
    C4_Model_close(edit_distance);
    return edit_distance;
    }

