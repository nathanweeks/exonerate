/****************************************************************\
*                                                                *
*  C4 dynamic programming library - optimal alignment code       *
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

#include <string.h>    /* For strlen() */

#include "slist.h"
#include "optimal.h"

Optimal *Optimal_create(C4_Model *model, gchar *name,
                        Optimal_Type type, gboolean use_codegen){
    register Optimal *optimal = g_new0(Optimal, 1);
    register gchar *viterbi_name;
    g_assert(model);
    g_assert(!model->is_open);
    optimal->ref_count = 1;
    if(name)
        optimal->name = g_strdup(name);
    else
        optimal->name = g_strdup_printf("optimal:%s", model->name);
    optimal->type = type;
    if(type & Optimal_Type_SCORE){
        viterbi_name = g_strconcat(optimal->name, " find score", NULL);
        optimal->find_score = Viterbi_create(model, viterbi_name,
                              Viterbi_Mode_FIND_SCORE, FALSE,
                              use_codegen);
        g_free(viterbi_name);
        }
    if(type & Optimal_Type_PATH){
        viterbi_name = g_strconcat(optimal->name, " find path", NULL);
        optimal->find_path = Viterbi_create(model, viterbi_name,
                             Viterbi_Mode_FIND_PATH, FALSE,
                             use_codegen);
        g_free(viterbi_name);
        }
    if(type & Optimal_Type_REDUCED_SPACE){
        if(!C4_Model_is_global(model)){
            viterbi_name = g_strconcat(optimal->name, " find region",
                                      NULL);
            optimal->find_region = Viterbi_create(model, viterbi_name,
                                      Viterbi_Mode_FIND_REGION, FALSE,
                                      use_codegen);
            g_free(viterbi_name);
            }
        viterbi_name = g_strconcat(optimal->name,
                                  " find checkpoint", NULL);
        optimal->find_checkpoint_continuation
           = Viterbi_create(model, viterbi_name,
                     Viterbi_Mode_FIND_CHECKPOINTS, TRUE, use_codegen);
        g_free(viterbi_name);
        viterbi_name = g_strconcat(optimal->name,
                                  " find path continuation", NULL);
        optimal->find_path_continuation
           = Viterbi_create(model, viterbi_name,
                     Viterbi_Mode_FIND_PATH, TRUE, use_codegen);
        g_free(viterbi_name);
        /**/
        }
    return optimal;
    }
/* FIXME: name should not be model->name + suffix,
 *        but optimal->name + suffix
 */

void Optimal_destroy(Optimal *optimal){
    g_assert(optimal);
    if(--optimal->ref_count)
        return;
    if(optimal->find_score)
        Viterbi_destroy(optimal->find_score);
    if(optimal->find_path)
        Viterbi_destroy(optimal->find_path);
    if(optimal->find_region)
        Viterbi_destroy(optimal->find_region);
    if(optimal->find_checkpoint_continuation)
        Viterbi_destroy(optimal->find_checkpoint_continuation);
    if(optimal->find_path_continuation)
        Viterbi_destroy(optimal->find_path_continuation);
    g_free(optimal->name);
    g_free(optimal);
    return;
    }

Optimal *Optimal_share(Optimal *optimal){
    g_assert(optimal);
    optimal->ref_count++;
    return optimal;
    }

static void Optimal_add_codegen(Viterbi *viterbi,
                                GPtrArray *codegen_list){
    register Codegen *codegen;
    if(viterbi){
        codegen = Viterbi_make_Codegen(viterbi);
        g_ptr_array_add(codegen_list, codegen);
        }
    return;
    }

GPtrArray *Optimal_make_Codegen_list(Optimal *optimal){
    register GPtrArray *codegen_list = g_ptr_array_new();
    Optimal_add_codegen(optimal->find_score, codegen_list);
    Optimal_add_codegen(optimal->find_path, codegen_list);
    Optimal_add_codegen(optimal->find_region, codegen_list);
    Optimal_add_codegen(optimal->find_checkpoint_continuation,
                        codegen_list);
    Optimal_add_codegen(optimal->find_path_continuation, codegen_list);
    g_assert(codegen_list->len);
    return codegen_list;
    }

C4_Score Optimal_find_score(Optimal *optimal, Region *region,
                            gpointer user_data, SubOpt *subopt){
    register C4_Score score;
    register Viterbi_Data *vd;
    g_assert(optimal->find_score);
    vd = Viterbi_Data_create(optimal->find_score, region);
    score = Viterbi_calculate(optimal->find_score, region,
                              vd, user_data, subopt);
    Viterbi_Data_destroy(vd);
    return score;
    }

static Region *Optimal_find_region(Optimal *optimal, Region *region,
                         gpointer user_data, C4_Score threshold,
                         SubOpt *subopt, C4_Score *region_score){
    register Region *alignment_region;
    register Viterbi_Data *vd;
    register C4_Score score;
    g_assert(optimal->type & Optimal_Type_REDUCED_SPACE);
    if(!optimal->find_region) /* Must be a global model */
        return Region_share(region); /* Already know region */
    g_assert(optimal->find_region);
    vd = Viterbi_Data_create(optimal->find_region, region);
    score = Viterbi_calculate(optimal->find_region, region,
                              vd, user_data, subopt);
    if(score < threshold)
        return NULL;
    alignment_region = Region_share(vd->alignment_region);
    (*region_score) = score;
    g_assert(Region_is_valid(alignment_region));
    g_assert(Region_is_within(region, alignment_region));
    Viterbi_Data_destroy(vd);
    return alignment_region;
    }

/**/

static C4_Score Optimal_find_checkpoints_continuation(Optimal *optimal,
                  Region *region, gpointer user_data,
                  SubOpt *subopt, GPtrArray **sub_vsa_list,
                  C4_State *first_state, C4_Score *first_cell,
                  C4_State *final_state, C4_Score *final_cell){
    register C4_Score score = C4_IMPOSSIBLY_LOW_SCORE;
    register Viterbi_Data *vd;
    vd = Viterbi_Data_create(optimal->find_checkpoint_continuation,
                             region);
    g_assert(vd->checkpoint);
    g_assert(!vd->continuation);
    Viterbi_Data_set_continuation(vd, first_state, first_cell,
                                      final_state, final_cell);
    g_assert(vd->continuation);
    score = Viterbi_calculate(optimal->find_checkpoint_continuation,
                              region, vd, user_data, subopt);
    (*sub_vsa_list) = Viterbi_Checkpoint_traceback(
                       optimal->find_checkpoint_continuation,
                       vd, region, first_state, final_cell);
    Viterbi_Data_destroy(vd);
    return score;
    }

static C4_Score Optimal_find_checkpoints_recur(Optimal *optimal,
                  Region *region, gpointer user_data,
                  SubOpt *subopt, SList *vsa_list,
                  C4_State *first_state, C4_Score *first_cell,
                  C4_State *final_state, C4_Score *final_cell){
    register Viterbi_SubAlignment *vsa, *prev_vsa = NULL,
                                        *next_vsa;
    register C4_Score *sub_first_cell;
    register C4_State *sub_final_state;
    register C4_Score score;
    register gint i;
    register gboolean used_vsa = TRUE, prev_used_vsa = TRUE;
    GPtrArray *sub_vsa_list = NULL;
    score = Optimal_find_checkpoints_continuation(
             optimal, region, user_data, subopt, &sub_vsa_list,
             first_state, first_cell, final_state, final_cell);
    g_assert(sub_vsa_list);
    for(i = sub_vsa_list->len-1; i >= 0; i--){
        vsa = sub_vsa_list->pdata[i];
        g_assert(vsa);
        if(Viterbi_use_reduced_space(optimal->find_path, vsa->region)){
            sub_first_cell = prev_vsa
                           ? prev_vsa->final_cell
                           : first_cell;
            if(i){
                next_vsa = sub_vsa_list->pdata[i-1];
                g_assert(next_vsa);
                sub_final_state = next_vsa->first_state;
            } else {
                sub_final_state = final_state;
                }
            Optimal_find_checkpoints_recur(optimal,
                 vsa->region, user_data, subopt, vsa_list,
                 vsa->first_state, sub_first_cell,
                 sub_final_state, vsa->final_cell);
            used_vsa = FALSE;
        } else {
            SList_queue(vsa_list, vsa);
            used_vsa = TRUE;
            }
        if(prev_vsa && (!prev_used_vsa))
            Viterbi_SubAlignment_destroy(prev_vsa);
        prev_vsa = vsa;
        prev_used_vsa = used_vsa;
        }
    g_ptr_array_free(sub_vsa_list, TRUE);
    return score;
    }

static Alignment *Optimal_find_path_quadratic_space_continuation(
           Optimal *optimal, Region *region,
           gpointer user_data, SubOpt *subopt,
           C4_State *first_state, C4_Score *first_cell,
           C4_State *final_state, C4_Score *final_cell){
    register Alignment *alignment;
    register C4_Score score;
    register Viterbi_Data *vd
        = Viterbi_Data_create(optimal->find_path_continuation, region);
    register AlignmentOperation *ao;
    g_assert(!Viterbi_use_reduced_space(optimal->find_path, region));
    Viterbi_Data_set_continuation(vd, first_state, first_cell,
                                      final_state, final_cell);
    score = Viterbi_calculate(optimal->find_path_continuation,
                              region, vd, user_data, subopt);
    g_assert(vd->curr_query_end == region->query_length);
    g_assert(vd->curr_target_end == region->target_length);
    g_assert(vd->continuation->final_cell[0] == score);
    alignment = Viterbi_Data_create_Alignment(vd,
                    optimal->find_path_continuation->model,
                    score, region);
    ao = alignment->operation_list->pdata[0];
    g_assert(ao);
    g_assert(first_state->id == ao->transition->input->id);
    ao = alignment->operation_list->pdata
        [alignment->operation_list->len-1];
    g_assert(ao);
    g_assert(final_state->id == ao->transition->output->id);
    Viterbi_Data_clear_continuation(vd);
    g_assert(score == alignment->score);
    Viterbi_Data_destroy(vd);
    return alignment;
    }

static Alignment *Optimal_compute_subalignments(Optimal *optimal,
                      Region *region, gpointer user_data,
                      SubOpt *subopt,
                      SList *vsa_list, C4_Score score, gint cell_size){
    register gint i;
    register Viterbi_SubAlignment *vsa = NULL;
    register Alignment *sub_alignment, *alignment;
    register AlignmentOperation *ao;
    register C4_State *sub_final_state;
    register C4_Score *sub_first_cell;
    register SListNode *node;
    register C4_Transition *transition;
    register C4_Score *dummy_cell = g_new0(C4_Score, cell_size);
    alignment = Alignment_create(optimal->find_path->model,
                                 region, score);
    SList_for_each(vsa_list, node){
        if(vsa){
            sub_first_cell = vsa->final_cell;
        } else {
            sub_first_cell = dummy_cell;
            }
        if(node->next){
            vsa = node->next->data;
            g_assert(vsa);
            sub_final_state = vsa->first_state;
        } else {
            sub_final_state = optimal->find_checkpoint_continuation
                            ->model->end_state->state;
            }
        vsa = node->data;
        g_assert(vsa);
        g_assert(vsa->region);
        sub_alignment = Optimal_find_path_quadratic_space_continuation(
            optimal, vsa->region, user_data, subopt,
            vsa->first_state, sub_first_cell,
            sub_final_state, vsa->final_cell);
        for(i = 0; i < sub_alignment->operation_list->len; i++){
            ao = sub_alignment->operation_list->pdata[i];
            /* Map transitions back to find_path model */
            transition = optimal->find_path->model->transition_list->pdata
                         [ao->transition->id];
            Alignment_add(alignment, transition, ao->length);
            }
        Alignment_destroy(sub_alignment);
        }
    g_free(dummy_cell);
    return alignment;
    }

static Alignment *Optimal_find_path_reduced_space(Optimal *optimal,
           Region *region, gpointer user_data, SubOpt *subopt){
    register Alignment *alignment;
    register SListSet *slist_set = SListSet_create();
    register SList *vsa_list = SList_create(slist_set);
    register SListNode *node;
    register gint cell_size
             = optimal->find_checkpoint_continuation->cell_size;
    register C4_Score score,
                     *dummy_first_cell = g_new0(C4_Score, cell_size),
                     *dummy_final_cell = g_new0(C4_Score, cell_size);
    register C4_Model *model;
    register Viterbi_SubAlignment *vsa;
    g_assert(Viterbi_use_reduced_space(optimal->find_path, region));
    model = optimal->find_checkpoint_continuation->model;
    score = Optimal_find_checkpoints_recur(optimal,
            region, user_data, subopt, vsa_list,
            model->start_state->state, dummy_first_cell,
            model->end_state->state, dummy_final_cell);
    alignment = Optimal_compute_subalignments(optimal, region,
                       user_data, subopt, vsa_list, score, cell_size);
    SList_for_each(vsa_list, node){
        vsa = node->data;
        Viterbi_SubAlignment_destroy(vsa);
        }
    SList_destroy(vsa_list);
    SListSet_destroy(slist_set);
    g_free(dummy_first_cell);
    g_free(dummy_final_cell);
    return alignment;
    }

/**/

static Alignment *Optimal_find_path_quadratic_space(
           Optimal *optimal, Region *region, gpointer user_data,
           SubOpt *subopt){
    register Alignment *alignment;
    register C4_Score score;
    register Viterbi_Data *vd = Viterbi_Data_create(optimal->find_path,
                                                    region);
    g_assert(!Viterbi_use_reduced_space(optimal->find_path, region));
    score = Viterbi_calculate(optimal->find_path, region,
                              vd, user_data, subopt);
    alignment = Viterbi_Data_create_Alignment(vd,
                     optimal->find_path->model, score, region);
    g_assert(score == alignment->score);
    Viterbi_Data_destroy(vd);
    return alignment;
    }

/**/

Alignment *Optimal_find_path(Optimal *optimal, Region *region,
                             gpointer user_data, C4_Score threshold,
                             SubOpt *subopt){
    register Alignment *alignment = NULL;
    register Region *alignment_region;
    C4_Score region_score = C4_IMPOSSIBLY_LOW_SCORE;
    /**/
    g_assert(!optimal->find_path->model->is_open);
    if(Viterbi_use_reduced_space(optimal->find_path, region)){
        alignment_region = Optimal_find_region(optimal, region,
                               user_data, threshold, subopt,
                               &region_score);
        if(!alignment_region) /* No region when score below threshold */
            return NULL;
        if(Viterbi_use_reduced_space(optimal->find_path,
                                     alignment_region)){
            alignment = Optimal_find_path_reduced_space(optimal,
                          alignment_region, user_data, subopt);
        } else {
            alignment = Optimal_find_path_quadratic_space(
                            optimal, alignment_region,
                            user_data, subopt);
            }
        g_assert(alignment);
        g_assert(Region_is_same(alignment_region, alignment->region));
#ifndef G_DISABLE_ASSERT
        if(alignment->score != region_score){
            g_warning("Region score DIFF have:[%d] expect:[%d]",
                      alignment->score, region_score);
            }
#endif /* G_DISABLE_ASSERT */
        g_assert(alignment->score == region_score);
        Region_destroy(alignment_region);
    } else {
        alignment = Optimal_find_path_quadratic_space(
                        optimal, region, user_data, subopt);
        }
    g_assert(optimal->find_path);
    g_assert(alignment->model == optimal->find_path->model);
    g_assert(Alignment_is_valid(alignment, region, user_data));
    if(alignment->score < threshold){
        Alignment_destroy(alignment);
        return NULL;
        }
    return alignment;
    }

