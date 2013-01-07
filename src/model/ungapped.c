/****************************************************************\
*                                                                *
*  Module for various ungapped alignment models                  *
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

#include "ungapped.h"

void Ungapped_Data_init(Ungapped_Data *ud,
                        Sequence *query, Sequence *target,
                        Match_Type match_type){
    g_assert(query);
    g_assert(target);
    if(!ud->query)
        ud->query = Sequence_share(query);
    if(!ud->target)
        ud->target = Sequence_share(target);
    if(!ud->mas)
        ud->mas = Match_ArgumentSet_create(NULL);
    if(!ud->match_list[match_type])
        ud->match_list[match_type] = Match_find(match_type);
    return;
    }

Ungapped_Data *Ungapped_Data_create(Sequence *query, Sequence *target,
                                    Match_Type match_type){
    register Ungapped_Data *ud = g_new0(Ungapped_Data, 1);
    Ungapped_Data_init(ud, query, target, match_type);
    return ud;
    }

void Ungapped_Data_clear(Ungapped_Data *ud){
    register gint i;
    if(ud->query){
        Sequence_destroy(ud->query);
        ud->query = NULL;
        }
    if(ud->target){
        Sequence_destroy(ud->target);
        ud->target = NULL;
        }
    for(i = 0; i < Match_Type_TOTAL; i++)
        if(ud->match_list[i])
            ud->match_list[i] = NULL;
    return;
    }

void Ungapped_Data_destroy(Ungapped_Data *ud){
    Ungapped_Data_clear(ud);
    g_free(ud);
    return;
    }

/**/

#define Ungapped_CalcFunc(type)                                    \
    static C4_Score Ungapped_calc_func_##type(gint query_pos,      \
                                              gint target_pos,     \
                                              gpointer user_data){ \
        register Ungapped_Data *ud = (Ungapped_Data*)user_data;    \
        register Match *match = ud->match_list[Match_Type_##type]; \
        g_assert(ud);                                              \
        g_assert(query_pos >= 0);                                  \
        g_assert(target_pos >= 0);                                 \
        g_assert(query_pos < ud->query->len);                      \
        g_assert(target_pos < ud->target->len);                    \
        g_assert(match);                                           \
        return match->score_func(match, ud->query, ud->target,     \
                                 query_pos, target_pos);           \
        }

#if 0
#define Ungapped_CalcMacro(type)                              \
    static gchar *Ungapped_calc_macro_##type =                \
        "ud->match_list[Match_Type_"#type"]->score_func(\n"   \
        "              ud->match_list[Match_Type_"#type"],\n" \
        "              ud->query, ud->target, %QP, %TP)";
/* FIXME: implement Ungapped_get_calc_macro() */
/* FIXME: add proper macro support for Match
 *        ie. without calling match->score_func()
 *        by calling Match_macro_create() etc
 */
/* Ungapped_CalcMacro(type) */
#endif /* 0 */

#define Ungapped_CalcElements(type) \
    Ungapped_CalcFunc(type)

Ungapped_CalcElements(DNA2DNA)
Ungapped_CalcElements(PROTEIN2PROTEIN)
Ungapped_CalcElements(DNA2PROTEIN)
Ungapped_CalcElements(PROTEIN2DNA)
Ungapped_CalcElements(CODON2CODON)

/**/

C4_Model *Ungapped_create(Match_Type match_type){
    register gchar *model_name;
    register C4_Model *ungapped;
    register C4_State *match_state;
    register C4_Calc *match_calc;
    register C4_CalcFunc calc_func = NULL;
    register gchar *calc_macro = Match_Type_get_score_macro(match_type);
    register Match *match = Match_find(match_type);
    model_name = g_strdup_printf("ungapped:%s",
                                 Match_Type_get_name(match->type));
    ungapped = C4_Model_create(model_name);
    g_free(model_name);
    match_state  = C4_Model_add_state(ungapped, "match");
    switch(match_type){
        case Match_Type_DNA2DNA:
            calc_func = Ungapped_calc_func_DNA2DNA;
            break;
        case Match_Type_PROTEIN2PROTEIN:
            calc_func = Ungapped_calc_func_PROTEIN2PROTEIN;
            break;
        case Match_Type_DNA2PROTEIN:
            calc_func = Ungapped_calc_func_DNA2PROTEIN;
            break;
        case Match_Type_PROTEIN2DNA:
            calc_func = Ungapped_calc_func_PROTEIN2DNA;
            break;
        case Match_Type_CODON2CODON:
            calc_func = Ungapped_calc_func_CODON2CODON;
            break;
        default:
            g_error("Unknown Match_Type [%d]", match_type);
            break;
        }
    match_calc = C4_Model_add_calc(ungapped, "match",
        Match_max_score(match), calc_func, calc_macro,
        NULL, NULL, NULL, NULL, C4_Protect_NONE);
/**/
    C4_Model_add_transition(ungapped, "start to match",
                 NULL, match_state, 0, 0, NULL, C4_Label_NONE, NULL);
    C4_Model_add_transition(ungapped, "match to end",
                 match_state, NULL, 0, 0, NULL, C4_Label_NONE, NULL);
    /* Transitions from match */
    C4_Model_add_transition(ungapped, "match",
                         match_state, match_state,
                         match->query->advance,
                         match->target->advance,
                         match_calc, C4_Label_MATCH, match);
    /* FIXME: need to keep match allocated */
/**/
    C4_Model_append_codegen(ungapped,
      "#include \"ungapped.h\"\n"
      "#include \"submat.h\"\n",
      "register Ungapped_Data *ud = (Ungapped_Data*)user_data;\n",
      " -I" SOURCE_ROOT_DIR "/src/model"
      " -I" SOURCE_ROOT_DIR "/src/comparison"
      " -I" SOURCE_ROOT_DIR "/src/c4"
      " -I" SOURCE_ROOT_DIR "/src/bsdp"
      " -I" SOURCE_ROOT_DIR "/src/sdp");
    C4_Model_close(ungapped);
    return ungapped;
    }

Alignment *Ungapped_Alignment_create(C4_Model *model,
                                     Ungapped_Data *ud, HSP *hsp){
    register gint i;
    register C4_Transition *transition, *start2match = NULL,
                           *match2match = NULL, *match2end = NULL;
    register Region *region = Region_create(hsp->query_start,
                                            hsp->target_start,
                       HSP_query_end(hsp)-hsp->query_start,
                       HSP_target_end(hsp)-hsp->target_start);
    register Alignment *alignment = Alignment_create(model,
                                           region, hsp->score);
    g_assert(model->transition_list->len == 3);
    for(i = 0; i < model->transition_list->len; i++){
        transition = model->transition_list->pdata[i];
        if(transition->input == model->start_state->state)
            start2match = transition;
        else if(transition->output == model->end_state->state)
            match2end = transition;
        else
            match2match = transition;
        }
    g_assert(start2match);
    g_assert(match2match);
    g_assert(match2end);
    Alignment_add(alignment, start2match, 1);
    Alignment_add(alignment, match2match, hsp->length);
    Alignment_add(alignment, match2end, 1);
    g_assert(Alignment_is_valid(alignment, region, ud));
    Region_destroy(region);
    return alignment;
    }

