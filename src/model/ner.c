/****************************************************************\
*                                                                *
*  Model for alignments with non-equivalenced regions            *
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

#include "ner.h"

NER_ArgumentSet *NER_ArgumentSet_create(Argument *arg){
    register ArgumentSet *as;
    static NER_ArgumentSet nas;
    if(arg){
        as = ArgumentSet_create("NER Model Options");
        ArgumentSet_add_option(as, 0, "minner", "length",
           "Minimum NER length", "10",
           Argument_parse_int, &nas.min_ner);
        ArgumentSet_add_option(as, 0, "maxner", "length",
           "Maximum NER length", "50000",
           Argument_parse_int, &nas.max_ner);
        ArgumentSet_add_option(as, 0, "neropen", "score",
           "NER open penalty", "-20",
           Argument_parse_int, &nas.ner_open_penalty);
        Argument_absorb_ArgumentSet(arg, as);
        }
    return &nas;
    }

NER_Data *NER_Data_create(Sequence *query, Sequence *target){
    register NER_Data *nd = g_new0(NER_Data, 1);
    Affine_Data_init(&nd->ad, query, target, FALSE);
    nd->nas = NER_ArgumentSet_create(NULL);
    return nd;
    }

void NER_Data_destroy(NER_Data *nd){
    Affine_Data_clear(&nd->ad);
    g_free(nd);
    return;
    }

/* Calc functions */

static C4_Score ner_ner_open_calc_func(gint query_pos,
                                               gint target_pos,
                                               gpointer user_data){
    register NER_Data *nd = (NER_Data*)user_data;
    return nd->nas->ner_open_penalty;
    }

static gchar *ner_ner_open_calc_macro
    = "(nd->nas->ner_open_penalty)\n";

/**/

C4_Model *NER_create(Alphabet_Type query_type,
                     Alphabet_Type target_type){
    register C4_Model *ner = Affine_create(Affine_Model_Type_LOCAL,
                                           query_type, target_type,
                                           FALSE);
    register C4_State *ner_state;
    register C4_Calc *ner_open_calc;
    register NER_ArgumentSet *nas
           = NER_ArgumentSet_create(NULL);
    register C4_Transition *match_transition;
    register gchar *model_name = g_strdup_printf("NER:%s", ner->name);
    C4_Model_rename(ner, model_name);
    g_free(model_name);
    C4_Model_open(ner);
    match_transition = C4_Model_select_single_transition(ner,
                                               C4_Label_MATCH);
    g_assert(match_transition);
    ner_state = C4_Model_add_state(ner, "ner");
    ner_open_calc = C4_Model_add_calc(ner, "ner open",
        nas->ner_open_penalty,
        ner_ner_open_calc_func,
        ner_ner_open_calc_macro,
        NULL, NULL, NULL, NULL, C4_Protect_NONE);
/**/
    /* Transitions from match */
    C4_Model_add_transition(ner, "match to ner",
                       match_transition->input, ner_state, 1, 1,
                       ner_open_calc, C4_Label_NER, NULL);
    /* Transitions from ner */
    C4_Model_add_transition(ner, "ner to match",
                       ner_state, match_transition->input, 0, 0,
                       NULL, C4_Label_NONE, NULL);
    C4_Model_add_transition(ner, "ner loop insert",
                       ner_state, ner_state, 1, 0,
                       NULL, C4_Label_NER, NULL);
    C4_Model_add_transition(ner, "ner loop delete",
                       ner_state, ner_state, 0, 1,
                       NULL, C4_Label_NER, NULL);
/**/
    C4_Model_add_span(ner, "ner span", ner_state,
                      nas->min_ner, nas->max_ner,
                      nas->min_ner, nas->max_ner);
/**/
    C4_Model_append_codegen(ner,
        "#include \"ner.h\"\n",
        "register NER_Data *nd = (NER_Data*)user_data;\n", NULL);
    C4_Model_close(ner);
    return ner;
    }

