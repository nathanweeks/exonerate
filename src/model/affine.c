/****************************************************************\
*                                                                *
*  Module for various affine gapped models                       *
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
#include "affine.h"

Affine_ArgumentSet *Affine_ArgumentSet_create(Argument *arg){
    register ArgumentSet *as;
    static Affine_ArgumentSet aas;
    if(arg){
        as = ArgumentSet_create("Affine Model Options");
        ArgumentSet_add_option(as, 'o', "gapopen", "penalty",
            "Affine gap open penalty", "-12",
            Argument_parse_int, &aas.gap_open);
        ArgumentSet_add_option(as, 'e', "gapextend", "penalty",
            "Affine gap extend penalty", "-4",
            Argument_parse_int, &aas.gap_extend);
        ArgumentSet_add_option(as, '\0', "codongapopen", "penalty",
            "Codon affine gap open penalty", "-18",
            Argument_parse_int, &aas.codon_gap_open);
        ArgumentSet_add_option(as, '\0', "codongapextend", "penalty",
            "Codon affine gap extend penalty", "-8",
            Argument_parse_int, &aas.codon_gap_extend);
        Argument_absorb_ArgumentSet(arg, as);
        if(aas.gap_open > 0)
            g_warning("Gap open penalty [%d] should be negative",
                    aas.gap_open);
        if(aas.gap_extend > 0)
            g_warning("Gap extend penalty [%d] should be negative",
                    aas.gap_extend);
        if(aas.codon_gap_open > 0)
            g_warning("Codon gap open penalty [%d] should be negative",
                    aas.codon_gap_open);
        if(aas.codon_gap_extend > 0)
            g_warning("Codon gap extend penalty [%d] should be negative",
                    aas.codon_gap_extend);
        }
    return &aas;
    }

/**/

void Affine_Data_init(Affine_Data *ad,
                      Sequence *query, Sequence *target,
                      gboolean translate_both){
    register Match_Type match_type
        = Match_Type_find(query->alphabet->type,
                          target->alphabet->type,
                          translate_both);
    Ungapped_Data_init(&ad->ud, query, target, match_type);
    if(!ad->aas)
        ad->aas = Affine_ArgumentSet_create(NULL);
    return;
    }

Affine_Data *Affine_Data_create(Sequence *query, Sequence *target,
                                gboolean translate_both){
    register Affine_Data *ad = g_new0(Affine_Data, 1);
    Affine_Data_init(ad, query, target, translate_both);
    return ad;
    }

void Affine_Data_clear(Affine_Data *ad){
    Ungapped_Data_clear(&ad->ud);
    return;
    }

void Affine_Data_destroy(Affine_Data *ad){
    Affine_Data_clear(ad);
    g_free(ad);
    return;
    }

/**/

static C4_Score affine_gap_open_calc_func(gint query_pos,
                                          gint target_pos,
                                          gpointer user_data){
    register Affine_Data *ad = (Affine_Data*)user_data;
    return ad->aas->gap_open;
    }

static gchar *affine_gap_open_calc_macro = "(ad->aas->gap_open)";

static C4_Score affine_gap_extend_calc_func(gint query_pos,
                                            gint target_pos,
                                            gpointer user_data){
    register Affine_Data *ad = (Affine_Data*)user_data;
    return ad->aas->gap_extend;
    }

static gchar *affine_gap_extend_calc_macro = "(ad->aas->gap_extend)";

/**/

static C4_Score affine_codon_gap_open_calc_func(gint query_pos,
                                                gint target_pos,
                                                gpointer user_data){
    register Affine_Data *ad = (Affine_Data*)user_data;
    return ad->aas->codon_gap_open;
    }

static gchar *affine_codon_gap_open_calc_macro = "(ad->aas->codon_gap_open)";

static C4_Score affine_codon_gap_extend_calc_func(gint query_pos,
                                                  gint target_pos,
                                                  gpointer user_data){
    register Affine_Data *ad = (Affine_Data*)user_data;
    return ad->aas->codon_gap_extend;
    }

static gchar *affine_codon_gap_extend_calc_macro = "(ad->aas->codon_gap_extend)";

/**/

gchar *Affine_Model_Type_get_name(Affine_Model_Type type){
    register gchar *model_name = NULL;
    switch(type){
        case Affine_Model_Type_GLOBAL:
            model_name = "global";
            break;
        case Affine_Model_Type_BESTFIT:
            model_name = "bestfit";
            break;
        case Affine_Model_Type_LOCAL:
            model_name = "local";
            break;
        case Affine_Model_Type_OVERLAP:
            model_name = "overlap";
            break;
        default:
            g_error("Unknown affine model type [%d]", type);
            break;
        }
    return model_name;
    }

C4_Model *Affine_create(Affine_Model_Type type,
                        Alphabet_Type query_type,
                        Alphabet_Type target_type,
                        gboolean translate_both){
    register C4_Scope scope = C4_Scope_ANYWHERE;
    register C4_Model *affine;
    register C4_State *insert_state, *delete_state;
    register C4_Transition *match_transition;
    register C4_Calc *gap_open_calc, *gap_extend_calc;
    register Affine_ArgumentSet *aas = Affine_ArgumentSet_create(NULL);
    register Match_Type match_type = Match_Type_find(query_type,
                                       target_type, translate_both);
    register C4_CalcFunc gap_open_calc_func = NULL,
                         gap_extend_calc_func = NULL;
    register gchar *gap_open_calc_macro = NULL,
                   *gap_extend_calc_macro = NULL;
    register gchar *model_name = Affine_Model_Type_get_name(type);
    affine = Ungapped_create(match_type);
    switch(type){
        case Affine_Model_Type_GLOBAL:
            scope = C4_Scope_CORNER;
            break;
        case Affine_Model_Type_BESTFIT:
            scope = C4_Scope_QUERY;
            break;
        case Affine_Model_Type_LOCAL:
            scope = C4_Scope_ANYWHERE;
            break;
        case Affine_Model_Type_OVERLAP:
            scope = C4_Scope_EDGE;
            break;
        default:
            g_error("Unknown affine model type [%d]", type);
            break;
        }
    g_assert(model_name);
    model_name = g_strdup_printf("affine:%s:%s", model_name,
                                 Match_Type_get_name(match_type));
    C4_Model_rename(affine, model_name);
    g_free(model_name);
    C4_Model_configure_start_state(affine, scope, NULL, NULL);
    C4_Model_configure_end_state(affine, scope, NULL, NULL);
    C4_Model_open(affine);
    insert_state = C4_Model_add_state(affine, "insert");
    delete_state = C4_Model_add_state(affine, "delete");
    match_transition = C4_Model_select_single_transition(affine,
                                                 C4_Label_MATCH);
    g_assert(match_transition);
    if(MAX(match_transition->advance_query,
           match_transition->advance_target) == 3){
        gap_open_calc_func = affine_codon_gap_open_calc_func;
        gap_extend_calc_func = affine_codon_gap_extend_calc_func;
        gap_open_calc_macro = affine_codon_gap_open_calc_macro;
        gap_extend_calc_macro = affine_codon_gap_extend_calc_macro;
    } else {
        gap_open_calc_func = affine_gap_open_calc_func;
        gap_extend_calc_func = affine_gap_extend_calc_func;
        gap_open_calc_macro = affine_gap_open_calc_macro;
        gap_extend_calc_macro = affine_gap_extend_calc_macro;
        }
    gap_open_calc = C4_Model_add_calc(affine, "gap open",
        aas->gap_open,
        gap_open_calc_func, gap_open_calc_macro,
        NULL, NULL, NULL, NULL, C4_Protect_NONE);
    gap_extend_calc = C4_Model_add_calc(affine, "gap extend",
        aas->gap_extend,
        gap_extend_calc_func, gap_extend_calc_macro,
        NULL, NULL, NULL, NULL, C4_Protect_NONE);
/**/
    /* Transitions from match */
    C4_Model_add_transition(affine, "match to insert",
                       match_transition->input, insert_state,
                       match_transition->advance_query, 0,
                       gap_open_calc, C4_Label_GAP, NULL);
    C4_Model_add_transition(affine, "match to delete",
                       match_transition->input, delete_state,
                       0, match_transition->advance_target,
                       gap_open_calc, C4_Label_GAP, NULL);
    /* Transitions from insert */
    C4_Model_add_transition(affine, "insert",
                    insert_state, insert_state,
                    match_transition->advance_query, 0,
                    gap_extend_calc, C4_Label_GAP, NULL);
    C4_Model_add_transition(affine, "insert to match",
                    insert_state, match_transition->output, 0, 0,
                    NULL, C4_Label_NONE, NULL);
    /* Transitions from delete */
    C4_Model_add_transition(affine, "delete",
                    delete_state, delete_state,
                    0, match_transition->advance_target,
                    gap_extend_calc, C4_Label_GAP, NULL);
    C4_Model_add_transition(affine, "delete to match",
                    delete_state, match_transition->output, 0, 0,
                    NULL, C4_Label_NONE, NULL);
/**/
    C4_Model_add_portal(affine, "match portal",
                        match_transition->calc,
                        match_transition->advance_query,
                        match_transition->advance_target);
/**/
    C4_Model_append_codegen(affine,
      "#include \"affine.h\"\n",
      "register Affine_Data *ad = (Affine_Data*)user_data;\n", NULL);
    C4_Model_close(affine);
    return affine;
    }

