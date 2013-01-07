/****************************************************************\
*                                                                *
*  Module for frameshift modelling                               *
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

#include "frameshift.h"
#include "ungapped.h"

Frameshift_ArgumentSet *Frameshift_ArgumentSet_create(Argument *arg){
    register ArgumentSet *as;
    static Frameshift_ArgumentSet fas;
    if(arg){
        as = ArgumentSet_create("Frameshift Options");
        ArgumentSet_add_option(as, 'f', "frameshift", "penalty",
            "Frameshift creation penalty", "-28",
            Argument_parse_int, &fas.frameshift_penalty);
        Argument_absorb_ArgumentSet(arg, as);
        if(fas.frameshift_penalty > 0)
            g_warning("Frameshift penalty should be negative [%d]",
                      fas.frameshift_penalty);
        }
    return &fas;
    }

/**/

Frameshift_Data *Frameshift_Data_create(void){
    register Frameshift_Data *fd = g_new0(Frameshift_Data, 1);
    fd->fas = Frameshift_ArgumentSet_create(NULL);
    return fd;
    }

void Frameshift_Data_destroy(Frameshift_Data *fd){
    g_free(fd);
    return;
    }

/**/

static C4_Score frameshift_penalty_calc_func(gint query_pos,
                                             gint target_pos,
                                             gpointer user_data){
    register Ungapped_Data *ud = user_data;
    register Frameshift_Data *fd
            = Ungapped_Data_get_Frameshift_Data(ud);
    return fd->fas->frameshift_penalty;
    }

static gchar *frameshift_penalty_calc_macro
    = "(fd->fas->frameshift_penalty)";

/**/

static C4_Calc *Frameshift_find_existing_calc(C4_Model *model){
    register gint i;
    register C4_Calc *calc;
    for(i = 0; i < model->calc_list->len; i++){
        calc = model->calc_list->pdata[i];
        if(calc->calc_func == frameshift_penalty_calc_func)
            return calc;
        }
    return NULL;
    }

void Frameshift_add(C4_Model *model, C4_State *match_state,
                    gchar *suffix, gboolean apply_to_query){
    register gchar *frameshift_name = g_strconcat("frameshift ",
                                                  suffix, NULL);
    register C4_State *frameshift_state = C4_Model_add_state(model,
                                              frameshift_name);
    register C4_Calc *frameshift_calc
             = Frameshift_find_existing_calc( model);
    register Frameshift_ArgumentSet *fas
           = Frameshift_ArgumentSet_create(NULL);
    register gchar *name;
    if(!frameshift_calc)
        frameshift_calc = C4_Model_add_calc(model, "frameshift",
                fas->frameshift_penalty,
                frameshift_penalty_calc_func,
                frameshift_penalty_calc_macro,
                NULL, NULL, NULL, NULL, C4_Protect_NONE);
    /* Transitions from match */
    name = g_strconcat("frameshift open 1 ", suffix, NULL);
    C4_Model_add_transition(model, name,
            match_state, frameshift_state,
            apply_to_query?1:0,
            apply_to_query?0:1,
            frameshift_calc, C4_Label_FRAMESHIFT, NULL);
    g_free(name);
    name = g_strconcat("frameshift open 2 ", suffix, NULL);
    C4_Model_add_transition(model, name,
            match_state, frameshift_state,
            apply_to_query?2:0,
            apply_to_query?0:2,
            frameshift_calc, C4_Label_FRAMESHIFT, NULL);
    g_free(name);
    /* Transition from frameshift */
    name = g_strconcat("frameshift close 0 ", suffix, NULL);
    C4_Model_add_transition(model, name,
            frameshift_state, match_state,
            0, 0, NULL, C4_Label_NONE, NULL);
    g_free(name);
    name = g_strconcat("frameshift close 3 ", suffix, NULL);
    C4_Model_add_transition(model, name,
            frameshift_state, match_state,
            apply_to_query?3:0,
            apply_to_query?0:3,
            NULL, C4_Label_FRAMESHIFT, NULL);
    g_free(name);
    C4_Model_append_codegen(model, NULL,
            "register Frameshift_Data *fd = ud->frameshift_data;\n",
            NULL);
    g_free(frameshift_name);
    return;
    }
/* A frameshift state is used, rather than using 0/4 and 0/5
 * transitions to keep max_advance for most models down to 3
 * as this is more efficient the current viterbi implementation.
 */

