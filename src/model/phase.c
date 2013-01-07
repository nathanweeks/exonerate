/****************************************************************\
*                                                                *
*  Model for phased introns.                                     *
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

#include "protein2genome.h"

/**/

/* 16:
 * Phase_{dna,protein}_phase{12}_{query,target}_calc_{func,macro}
 */

/*...*/

#if 0

static gchar Phase_calc_get_aa(Sequence *seq, gint pos, gint phase,
                 gboolean has_intron, Translate *translate,
                 Intron_ChainData *intron_chain_data){
    register gchar aa = '\0';
    register gint codon_1, codon_2, codon_3;
    g_assert(has_intron?(seq->alphabet->type == Alphabet_Type_DNA)
                       :TRUE);
    g_assert(seq->alphabet->type != Alphabet_Type_UNKNOWN);
    if(seq->alphabet->type == Alphabet_Type_PROTEIN){
        g_assert(!has_intron);
        aa = Sequence_get_symbol(seq, pos);
    } else { /* seq->type == Alphabet_Type_DNA */
        if(phase == 1){
            if(has_intron){
                g_assert(intron_chain_data->curr_intron_start >= 0);
                codon_1 = intron_chain_data->curr_intron_start - 1;
            } else {
                codon_1 = pos-1;
                }
            codon_2 = pos;
            codon_3 = pos+1;
        } else { /* phase == 2 */
            if(has_intron){
                g_assert(intron_chain_data->curr_intron_start >= 0);
                codon_1 = intron_chain_data->curr_intron_start - 2;
                codon_2 = intron_chain_data->curr_intron_start - 1;
            } else {
                codon_1 = pos-2;
                codon_2 = pos-1;
                }
            codon_3 = pos;
            }
        /**/
        g_assert(codon_1 >= 0);
        g_assert(codon_1 < codon_2);
        g_assert(codon_2 < codon_3);
        g_assert(codon_3 < seq->len);
        aa = Translate_base(translate, Sequence_get_symbol(seq, codon_1),
                                       Sequence_get_symbol(seq, codon_2),
                                       Sequence_get_symbol(seq, codon_3));
        }
    g_assert(aa);
    return aa;
    }

static gboolean Phase_calc_is_valid(Sequence *seq, gint pos, gint phase,
       gboolean has_intron, Intron_ChainData *intron_chain_data){
    if(seq->alphabet->type == Alphabet_Type_PROTEIN)
        return TRUE;
    if(has_intron){
        if(intron_chain_data->curr_intron_start < phase)
            return FALSE;
    } else {
        if(pos < phase)
            return FALSE;
        }
    return TRUE;
    }

static C4_Score Phase_calc_general(gint query_pos, gint target_pos,
                    gpointer user_data, gint phase,
                    gboolean intron_on_query){
    register Ungapped_Data *ud = user_data;
    register Sequence *query = Ungapped_Data_get_query(ud),
                      *target = Ungapped_Data_get_target(ud);
    register Translate *translate = Ungapped_Data_get_translate(ud);
    register Submat *submat = Ungapped_Data_get_protein_submat(ud);
    register gchar qy_aa, tg_aa;
    if((!Phase_calc_is_valid(query, query_pos, phase,
             intron_on_query, ud->intron_data->query_data))
    || (!Phase_calc_is_valid(target, target_pos, phase,
             !intron_on_query, ud->intron_data->target_data)))
        return C4_IMPOSSIBLY_LOW_SCORE;
    qy_aa = Phase_calc_get_aa(query, query_pos, phase,
                intron_on_query, translate,
                ud->intron_data->query_data),
    tg_aa = Phase_calc_get_aa(target, target_pos, phase,
                !intron_on_query, translate,
                ud->intron_data->target_data);
    return Submat_lookup(submat, qy_aa, tg_aa);
    }

/* Phase_{dna,protein}_phase{12}_{query,target}_calc_{func,macro} */

#define Phase_CalcFunc(type, phase, chain, intron_on_query)           \
static C4_Score Phase_##type##_phase##phase##_##chain##_calc_func(    \
                gint query_pos, gint target_pos, gpointer user_data){ \
    return Phase_calc_general(query_pos, target_pos, user_data,       \
                              phase, intron_on_query);                \
    }

/**/
Phase_CalcFunc(protein, 1, query, TRUE, PROTEIN2DNA)
Phase_CalcFunc(dna,     1, query, TRUE, CODON2CODON)
Phase_CalcFunc(protein, 2, query, TRUE, PROTEIN2DNA)
Phase_CalcFunc(dna,     2, query, TRUE, CODON2CODON)
/**/
Phase_CalcFunc(protein, 1, target, FALSE, PROTEIN2DNA)
Phase_CalcFunc(dna,     1, target, FALSE, CODON2CODON)
Phase_CalcFunc(protein, 2, target, FALSE, PROTEIN2DNA)
Phase_CalcFunc(dna,     2, target, FALSE, CODON2CODON)
/**/

#endif /* 0 */

static void Phase_calc_get_positions(Sequence *seq, gint pos, gint phase,
                 gboolean has_intron, Intron_ChainData *intron_chain_data,
                 gint *p1, gint *p2, gint *p3){
    g_assert(has_intron?(seq->alphabet->type == Alphabet_Type_DNA):TRUE);
    g_assert(seq->alphabet->type != Alphabet_Type_UNKNOWN);
    if(seq->alphabet->type == Alphabet_Type_PROTEIN){
        g_assert(!has_intron);
        (*p1) = pos;
    } else { /* seq->type == Alphabet_Type_DNA */
        g_assert(seq->alphabet->type == Alphabet_Type_DNA);
        if(phase == 1){
            if(has_intron){
                g_assert(intron_chain_data->curr_intron_start >= 0);
                (*p1) = intron_chain_data->curr_intron_start - 1;
            } else {
                (*p1) = pos-1;
                }
            (*p2) = pos;
            (*p3) = pos+1;
        } else { /* phase == 2 */
            if(has_intron){
                g_assert(intron_chain_data->curr_intron_start >= 0);
                (*p1) = intron_chain_data->curr_intron_start - 2;
                (*p2) = intron_chain_data->curr_intron_start - 1;
            } else {
                (*p1) = pos-2;
                (*p2) = pos-1;
                }
            (*p3) = pos;
            }
        /**/
        g_assert((*p1) >= 0);
        g_assert((*p1) < (*p2));
        g_assert((*p2) < (*p3));
        g_assert((*p3) < seq->len);
        }
    return;
    }

static gboolean Phase_calc_is_valid(Sequence *seq, gint pos, gint phase,
       gboolean has_intron, Intron_ChainData *intron_chain_data){
    if(seq->alphabet->type == Alphabet_Type_PROTEIN)
        return TRUE;
    if(has_intron){
        if(intron_chain_data->curr_intron_start < phase)
            return FALSE;
    } else {
        if(pos < phase)
            return FALSE;
        }
    return TRUE;
    }

#define Phase_CalcFunc(type, phase, qy_intron, tg_intron)              \
static C4_Score                                                        \
Phase_##phase##_##type##_##qy_intron##_##tg_intron##_calc_func(        \
                gint query_pos, gint target_pos, gpointer user_data){  \
    register Ungapped_Data *ud = user_data;                            \
    register Match *match = ud->match_list[Match_Type_##type];         \
    gint qp1, qp2, qp3, tp1, tp2, tp3;                                 \
    if((!Phase_calc_is_valid(ud->query, query_pos, phase,              \
             qy_intron, ud->intron_data->query_data))                  \
    || (!Phase_calc_is_valid(ud->target, target_pos, phase,            \
             tg_intron, ud->intron_data->target_data)))                \
        return C4_IMPOSSIBLY_LOW_SCORE;                                \
    Phase_calc_get_positions(ud->query, query_pos, phase, qy_intron,   \
                             ud->intron_data->query_data,              \
                             &qp1, &qp2, &qp3);                        \
    Phase_calc_get_positions(ud->target, target_pos, phase, tg_intron, \
                             ud->intron_data->target_data,             \
                             &tp1, &tp2, &tp3);                        \
    return match->split_score_func(match, ud->query, ud->target,       \
                                   qp1, qp2, qp3, tp1, tp2, tp3);      \
    }

Phase_CalcFunc(PROTEIN2DNA, 1, FALSE, TRUE)
Phase_CalcFunc(PROTEIN2DNA, 2, FALSE, TRUE)
Phase_CalcFunc(DNA2PROTEIN, 1, TRUE, FALSE)
Phase_CalcFunc(DNA2PROTEIN, 2, TRUE, FALSE)

Phase_CalcFunc(CODON2CODON, 1, TRUE, TRUE)
Phase_CalcFunc(CODON2CODON, 2, TRUE, TRUE)
Phase_CalcFunc(CODON2CODON, 1, TRUE, FALSE)
Phase_CalcFunc(CODON2CODON, 2, TRUE, FALSE)
Phase_CalcFunc(CODON2CODON, 1, FALSE, TRUE)
Phase_CalcFunc(CODON2CODON, 2, FALSE, TRUE)

static gchar *Phase_get_calc_aa_macro(Alphabet_Type type, gint phase,
                  gboolean is_query, gboolean has_intron){
    register gchar *macro, *codon_1, *codon_2, *codon_3;
    register gchar *chain_name = is_query?"query":"target";
    register gchar chain_symbol = is_query?'Q':'T';
    g_assert(has_intron?(type == Alphabet_Type_DNA):TRUE);
    g_assert(type != Alphabet_Type_UNKNOWN);
    if(type == Alphabet_Type_PROTEIN){
        g_assert(!has_intron);
        /*
        macro = g_strdup_printf("Sequence_get_symbol(ud->%s, %%%cP)",
                    chain_name, chain_symbol);
        */
        macro = g_strdup_printf("%%%cP, 0, 0", chain_symbol);
    } else { /* query is DNA */
        if(phase == 1){
            if(has_intron){
                codon_1 = g_strdup_printf(
                              "(id->%s_data->curr_intron_start - 1)",
                              chain_name);
            } else {
                codon_1 = g_strdup_printf("%%%cP-1", chain_symbol);
                }
            codon_2 = g_strdup_printf("%%%cP", chain_symbol);
            codon_3 = g_strdup_printf("%%%cP+1", chain_symbol);
        } else { /* phase == 2 */
            if(has_intron){
                codon_1 = g_strdup_printf(
                              "(id->%s_data->curr_intron_start - 2)",
                               chain_name);
                codon_2 = g_strdup_printf(
                              "(id->%s_data->curr_intron_start - 1)",
                               chain_name);
            } else {
                codon_1 = g_strdup_printf("%%%cP-2", chain_symbol);
                codon_2 = g_strdup_printf("%%%cP-1", chain_symbol);
                }
            codon_3 = g_strdup_printf("%%%cP", chain_symbol);
            }
        g_assert(codon_1);
        g_assert(codon_2);
        g_assert(codon_3);
        macro = g_strdup_printf("%s, %s, %s", codon_1, codon_2, codon_3);
        /*
        macro = g_strdup_printf("Translate_base(ud->mas->translate,\n"
                                "    Sequence_get_symbol(ud->%s, %s]),\n"
                                "    Sequence_get_symbol(ud->%s, %s]),\n"
                                "    Sequence_get_symbol(ud->%s, %s]))\n",
                                chain_name, codon_1,
                                chain_name, codon_2,
                                chain_name, codon_3);
        */
        g_free(codon_1);
        g_free(codon_2);
        g_free(codon_3);
        }
    return macro;
    }

static gchar *Phase_get_calc_is_valid_macro_chain(Alphabet_Type type,
                  gboolean is_query, gint phase, gboolean has_intron){
    register gchar *macro = NULL;
    if(type == Alphabet_Type_PROTEIN)
        return NULL;
    if(has_intron){
        macro = g_strdup_printf(
                    "(id->%s_data->curr_intron_start >= %d)",
                    is_query?"query":"target", phase);
    } else {
        macro = g_strdup_printf("(%%%cP >= %d)", is_query?'Q':'T', phase);
        }
    return macro;
    }

static gchar *Phase_get_calc_is_valid_macro(Match *match,
                                            gint phase,
                                            gboolean on_query,
                                            gboolean on_target){
    register gchar *query_macro, *target_macro, *macro = NULL;
    query_macro = Phase_get_calc_is_valid_macro_chain(
                      match->query->alphabet->type, TRUE, phase, on_query);
    target_macro = Phase_get_calc_is_valid_macro_chain(
                      match->target->alphabet->type, FALSE, phase, on_target);
    if(query_macro){
        if(target_macro){ /* QT */
            macro = g_strdup_printf("(%s && %s)",
                        query_macro, target_macro);
            g_free(query_macro);
            g_free(target_macro);
        } else { /* Q- */
            macro = query_macro;
            }
    } else {
        if(target_macro){ /* -T */
            macro = target_macro;
        } else { /* -- */
            g_assert_not_reached();
            }
        }
    return macro;
    }

static gchar *Phase_get_calc_macro(Match *match, gint phase,
                                   gboolean on_query, gboolean on_target){
    register gchar *qy_pos_macro, *tg_pos_macro, *macro;
    register gchar *is_valid_macro
        = Phase_get_calc_is_valid_macro(match, phase, on_query, on_target);
    qy_pos_macro = Phase_get_calc_aa_macro(match->query->alphabet->type,
                                           phase, TRUE, on_query);
    tg_pos_macro = Phase_get_calc_aa_macro(match->target->alphabet->type,
                                           phase, FALSE, on_target);
    macro = g_strdup_printf(
              "(%s)?ud->match_list[%d]->split_score_func(\n"
              "ud->match_list[%d], ud->query, ud->target,\n"
              "%s, %s):%d",
              is_valid_macro,
              match->type, match->type,
              qy_pos_macro, tg_pos_macro,
              C4_IMPOSSIBLY_LOW_SCORE);
    g_free(is_valid_macro);
    g_free(qy_pos_macro);
    g_free(tg_pos_macro);
    return macro;
    }

/**/

/*           p2d        d2d
   phase1 0,1...1,2  1,1...2,2
   phase2 0,2...1,1  2,2...1,1
*/

C4_Model *Phase_create(gchar *suffix, Match *match,
                       gboolean on_query, gboolean on_target){
    register gchar *name, *intron_name, *calc_name;
    register C4_Model *model,
        *intron_model_00, *intron_model_12, *intron_model_21;
    register C4_Calc *phase1post_to_dst_calc,
                     *phase2post_to_dst_calc;
    register C4_State *phase1pre_state, *phase1post_state,
                      *phase2pre_state, *phase2post_state;
    register gint pre1_qa = 0, pre1_ta = 0, post2_qa = 0, post2_ta = 0,
                  pre2_qa = 0, pre2_ta = 0, post1_qa = 0, post1_ta = 0;
    register C4_CalcFunc phase1_calc_func = NULL,
                         phase2_calc_func = NULL;
    register C4_Transition *phase1post_transition,
                           *phase2post_transition;
    register C4_Shadow *shadow;
    register gchar *phase1_calc_macro = NULL,
                   *phase2_calc_macro = NULL;
    register gchar *full_suffix = NULL;
    register gboolean against_peptide
        = ((match->query->alphabet->type == Alphabet_Type_PROTEIN)
        || (match->target->alphabet->type == Alphabet_Type_PROTEIN));
    g_assert(on_query || on_target);
    /* Cannot have joint introns vs peptide */
    g_assert((!(on_query && on_target)) || (!against_peptide));
    full_suffix = g_strconcat("phase", suffix?" ":"", suffix?suffix:"",
        suffix?" ":"", on_query?"Q":"-", on_target?"T":"-", NULL);
    model = C4_Model_create(full_suffix);
    /* Prepare intron models */
    intron_name = g_strconcat("0:0 ", full_suffix, NULL);
    intron_model_00 = Intron_create(intron_name,
                                    on_query, on_target, TRUE);
    intron_name[0] = '1';
    intron_name[2] = '2';
    intron_model_12 = Intron_create(intron_name,
                                    on_query, on_target, TRUE);
    intron_name[0] = '2';
    intron_name[2] = '1';
    intron_model_21 = Intron_create(intron_name,
                                    on_query, on_target, TRUE);
    g_free(intron_name);
    if(against_peptide){
        if(on_query){
            pre1_qa  = 1;
            pre1_ta  = 0;
            post1_qa = 2;
            post1_ta = 1;
            pre2_qa  = 2;
            pre2_ta  = 0;
            post2_qa = 1;
            post2_ta = 1;
        } else { /* on_target */
            pre1_qa  = 0;
            pre1_ta  = 1;
            post1_qa = 1;
            post1_ta = 2;
            pre2_qa  = 0;
            pre2_ta  = 2;
            post2_qa = 1;
            post2_ta = 1;
            }
    } else { /* Against DNA */
        pre1_qa  = 1;
        pre1_ta  = 1;
        post1_qa = 2;
        post1_ta = 2;
        pre2_qa  = 2;
        pre2_ta  = 2;
        post2_qa = 1;
        post2_ta = 1;
        }
    /* Select C4 calc functions */
    switch(match->type){
        case Match_Type_PROTEIN2DNA:
            g_assert(!on_query);
            g_assert(on_target);
            phase1_calc_func = Phase_1_PROTEIN2DNA_FALSE_TRUE_calc_func;
            phase2_calc_func = Phase_2_PROTEIN2DNA_FALSE_TRUE_calc_func;
            break;
        case Match_Type_DNA2PROTEIN:
            g_assert(on_query); /* FAILS */
            g_assert(!on_target);
            phase1_calc_func = Phase_1_DNA2PROTEIN_TRUE_FALSE_calc_func;
            phase2_calc_func = Phase_2_DNA2PROTEIN_TRUE_FALSE_calc_func;
            break;
        case Match_Type_CODON2CODON:
            g_assert(on_query || on_target);
            if(on_query){
                if(on_target){ /* Joint */
                    phase1_calc_func = Phase_1_CODON2CODON_TRUE_TRUE_calc_func;
                    phase2_calc_func = Phase_2_CODON2CODON_TRUE_TRUE_calc_func;
                } else { /* Query */
                    phase1_calc_func
                        = Phase_1_CODON2CODON_TRUE_FALSE_calc_func;
                    phase2_calc_func
                        = Phase_2_CODON2CODON_TRUE_FALSE_calc_func;
                    }
            } else { /* Target */
                g_assert(on_target);
                phase1_calc_func = Phase_1_CODON2CODON_FALSE_TRUE_calc_func;
                phase2_calc_func = Phase_2_CODON2CODON_FALSE_TRUE_calc_func;
                }
            break;
        default:
            g_error("Match_Type [%s] unsuitable for phase model",
                    Match_Type_get_name(match->type));
            break;
        }
    phase1_calc_macro = Phase_get_calc_macro(match, 1, on_query, on_target);
    phase2_calc_macro = Phase_get_calc_macro(match, 2, on_query, on_target);
    g_assert(phase1_calc_macro);
    g_assert(phase2_calc_macro);
    /* Add calcs */
    calc_name = g_strconcat("phase1post to dst ", full_suffix, NULL);
    phase1post_to_dst_calc = C4_Model_add_calc(model,
            calc_name, Match_max_score(match),
            phase1_calc_func, phase1_calc_macro,
            NULL, NULL, NULL, NULL, C4_Protect_NONE);
    g_free(calc_name);
    calc_name = g_strconcat("phase2post to dst ", full_suffix, NULL);
    phase2post_to_dst_calc = C4_Model_add_calc(model,
            calc_name, Match_max_score(match),
            phase2_calc_func, phase2_calc_macro,
            NULL, NULL, NULL, NULL, C4_Protect_NONE);
    g_free(calc_name);
    g_free(phase1_calc_macro);
    g_free(phase2_calc_macro);
    /* Add split codon states */
    name = g_strconcat("phase1pre ", full_suffix, NULL);
    phase1pre_state = C4_Model_add_state(model, name);
    g_free(name);
    name = g_strconcat("phase1post ", full_suffix, NULL);
    phase1post_state = C4_Model_add_state(model, name);
    g_free(name);
    name = g_strconcat("phase2pre ", full_suffix, NULL);
    phase2pre_state = C4_Model_add_state(model, name);
    g_free(name);
    name = g_strconcat("phase2post ", full_suffix, NULL);
    phase2post_state = C4_Model_add_state(model, name);
    g_free(name);
    /* Transitions from [START]->phase_pre transitions */
    name = g_strconcat("(START) to ", phase1pre_state->name,
                      NULL);
    C4_Model_add_transition(model, name, NULL, phase1pre_state,
        pre1_qa, pre1_ta, NULL, C4_Label_SPLIT_CODON, NULL);
    g_free(name);
    name = g_strconcat("(START) to ", phase2pre_state->name, NULL);
    C4_Model_add_transition(model, name, NULL, phase2pre_state,
        pre2_qa, pre2_ta, NULL, C4_Label_SPLIT_CODON, NULL);
    g_free(name);
    /* Transitions from post->[END] transitions. */
    name = g_strconcat(phase1post_state->name, " to (END)", NULL);
    phase1post_transition = C4_Model_add_transition(model, name,
            phase1post_state, NULL, post1_qa, post1_ta,
            phase1post_to_dst_calc, C4_Label_SPLIT_CODON, NULL);
    g_free(name);
    name = g_strconcat(phase2post_state->name, " to (END)", NULL);
    phase2post_transition = C4_Model_add_transition(model, name,
            phase2post_state, NULL, post2_qa, post2_ta,
            phase2post_to_dst_calc, C4_Label_SPLIT_CODON, NULL);
    g_free(name);
    /* Insert intron models */
    C4_Model_insert(model, intron_model_00, NULL, NULL);
    C4_Model_insert(model, intron_model_12, phase1pre_state,
                                            phase1post_state);
    C4_Model_insert(model, intron_model_21, phase2pre_state,
                                            phase2post_state);
    /* Add shadows */
    if(on_query && on_target){
        g_assert(model->shadow_list->len == 6);
        shadow = model->shadow_list->pdata[2];
        C4_Shadow_add_dst_transition(shadow, phase1post_transition);
        shadow = model->shadow_list->pdata[3];
        C4_Shadow_add_dst_transition(shadow, phase1post_transition);
        shadow = model->shadow_list->pdata[4];
        C4_Shadow_add_dst_transition(shadow, phase2post_transition);
        shadow = model->shadow_list->pdata[5];
        C4_Shadow_add_dst_transition(shadow, phase2post_transition);
    } else {
        g_assert(model->shadow_list->len == 3);
        shadow = model->shadow_list->pdata[1];
        C4_Shadow_add_dst_transition(shadow, phase1post_transition);
        shadow = model->shadow_list->pdata[2];
        C4_Shadow_add_dst_transition(shadow, phase2post_transition);
        }
    /**/
    C4_Model_destroy(intron_model_00);
    C4_Model_destroy(intron_model_12);
    C4_Model_destroy(intron_model_21);
    /**/
    C4_Model_close(model);
    g_free(full_suffix);
    return model;
    }

/**/

/* FIXME: tidy */

