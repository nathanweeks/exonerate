/****************************************************************\
*                                                                *
*  C4 dynamic programming library - Seeded Dynamic Programming   *
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

#include <stdlib.h> /* For qsort() */

#include "sdp.h"
#include "matrix.h"

/**/

SDP_ArgumentSet *SDP_ArgumentSet_create(Argument *arg){
    register ArgumentSet *as;
    static SDP_ArgumentSet sas = {30};
    if(arg){
        as = ArgumentSet_create("Seeded Dynamic Programming options");
        ArgumentSet_add_option(as, 'x', "extensionthreshold", NULL,
                "Gapped extension threshold", "50",
                Argument_parse_int, &sas.dropoff);
        ArgumentSet_add_option(as, 0, "singlepass", NULL,
                "Generate suboptimal alignment in a single pass", "TRUE",
                Argument_parse_boolean, &sas.single_pass_subopt);
        Argument_absorb_ArgumentSet(arg, as);
        }
    return &sas;
    }

/**/

typedef struct {
     GPtrArray *seed_list;
          gint  seed_pos;
          gint  seed_state_id;
      SDP_Pair *sdp_pair;
} Scheduler_Seed_List;

static Scheduler_Seed_List *Scheduler_Seed_List_create(SDP_Pair *sdp_pair,
                                                       gboolean is_forward){
    register Scheduler_Seed_List *seed_info_list
        = g_new(Scheduler_Seed_List, 1);
    register C4_Model *model = sdp_pair->sdp->model;
    seed_info_list->seed_list = sdp_pair->seed_list;
    seed_info_list->seed_pos = 0;
    if(is_forward){
        seed_info_list->seed_state_id = model->start_state->state->id;
    } else {
        seed_info_list->seed_state_id = model->end_state->state->id;
        }
    seed_info_list->sdp_pair = sdp_pair;
    return seed_info_list;
    }

static void Scheduler_Seed_List_destroy(
            Scheduler_Seed_List *seed_info_list){
    g_free(seed_info_list);
    return;
    }

static void Scheduler_Seed_List_init_forward(gpointer seed_data){
    register Scheduler_Seed_List *seed_info_list = seed_data;
    g_assert(seed_info_list);
    seed_info_list->seed_pos = 0;
    return;
    }

static void Scheduler_Seed_List_init_reverse(gpointer seed_data){
    register Scheduler_Seed_List *seed_info_list = seed_data;
    g_assert(seed_info_list);
    seed_info_list->seed_pos = seed_info_list->seed_list->len - 1;
    return;
    }

static void Scheduler_Seed_List_next_forward(gpointer seed_data){
    register Scheduler_Seed_List *seed_info_list = seed_data;
    g_assert(seed_info_list);
    g_assert(seed_info_list->seed_pos < seed_info_list->seed_list->len);
    seed_info_list->seed_pos++;
    return;
    }

static void Scheduler_Seed_List_next_reverse(gpointer seed_data){
    register Scheduler_Seed_List *seed_info_list = seed_data;
    g_assert(seed_info_list);
    g_assert(seed_info_list->seed_pos >= 0);
    seed_info_list->seed_pos--;
    return;
    }

static gboolean Scheduler_Seed_List_get_forward(gpointer seed_data,
                                      Scheduler_Seed *seed_info){
    register Scheduler_Seed_List *seed_info_list = seed_data;
    register SDP_Seed *seed;
    g_assert(seed_info_list);
    if(seed_info_list->seed_pos >= seed_info_list->seed_list->len)
        return FALSE;
    seed = seed_info_list->seed_list->pdata[seed_info_list->seed_pos];
    seed_info->query_pos = HSP_query_cobs(seed->hsp);
    seed_info->target_pos = HSP_target_cobs(seed->hsp);
    seed_info->seed_id = seed->seed_id;
    seed_info->start_score = seed->max_start->score
                          - (seed->hsp->score >> 1);
    return TRUE;
    }

static gboolean Scheduler_Seed_List_get_reverse(gpointer seed_data,
                                      Scheduler_Seed *seed_info){
    register Scheduler_Seed_List *seed_info_list = seed_data;
    register SDP_Seed *seed;
    g_assert(seed_info_list);
    if(seed_info_list->seed_pos < 0)
        return FALSE;
    seed = seed_info_list->seed_list->pdata[seed_info_list->seed_pos];
    seed_info->query_pos = -HSP_query_cobs(seed->hsp);
    seed_info->target_pos = -HSP_target_cobs(seed->hsp);
    seed_info->seed_id = seed->seed_id;
    seed_info->start_score = (seed->hsp->score >> 1);
    return TRUE;
    }

static void Scheduler_Seed_List_start_func(gpointer seed_data,
                                           gint seed_id, C4_Score score,
                                           gint start_query_pos,
                                           gint start_target_pos,
                                           STraceback_Cell *stcell){
    register Scheduler_Seed_List *seed_info_list = seed_data;
    register SDP_Seed *seed;
    g_assert(seed_info_list);
    g_assert(seed_info_list->seed_list);
    g_assert(seed_id < seed_info_list->sdp_pair->seed_list->len);
    seed = seed_info_list->sdp_pair->seed_list->pdata[seed_id];
    /* Is best start for this seed */
    if(seed->max_start->score < score){
        seed->max_start->score = score;
        seed->max_start->query_pos = start_query_pos;
        seed->max_start->target_pos = start_target_pos;
        if(stcell)
            seed->max_start->cell = STraceback_Cell_share(stcell);
        else
            seed->max_start->cell = NULL;
        }
    return;
    }

static void Scheduler_Seed_List_end_func(gpointer seed_data,
                                         gint seed_id, C4_Score score,
                                         gint end_query_pos,
                                         gint end_target_pos,
                                         STraceback_Cell *stcell){
    register Scheduler_Seed_List *seed_info_list = seed_data;
    register SDP_Seed *seed;
    g_assert(seed_info_list);
    g_assert(seed_info_list->seed_list);
    g_assert(seed_id < seed_info_list->sdp_pair->seed_list->len);
    seed = seed_info_list->sdp_pair->seed_list->pdata[seed_id];
    /* Is best end for this seed */
    if(seed->max_end->score < score){
        seed->max_end->score = score;
        seed->max_end->query_pos = end_query_pos;
        seed->max_end->target_pos = end_target_pos;
        g_assert(stcell);
        seed->max_end->cell = STraceback_Cell_share(stcell);
        }
    return;
    }

/**/

typedef struct {
      Boundary *boundary;
          gint  curr_row;
          gint  curr_interval;
          gint  curr_pos;
      gboolean  is_finished;
          gint  seed_state_id;
      SDP_Pair *sdp_pair;
} Scheduler_Seed_Boundary;

static Scheduler_Seed_Boundary *Scheduler_Seed_Boundary_create(
                                 Boundary *boundary,
                                 SDP_Pair *sdp_pair){
    register Scheduler_Seed_Boundary *seed_info_boundary
        = g_new(Scheduler_Seed_Boundary, 1);
    seed_info_boundary->boundary = Boundary_share(boundary);
    seed_info_boundary->curr_row = 0;
    seed_info_boundary->curr_interval = 0;
    seed_info_boundary->curr_pos = 0;
    seed_info_boundary->is_finished = FALSE;
    seed_info_boundary->seed_state_id
        = sdp_pair->sdp->model->start_state->state->id;
    seed_info_boundary->sdp_pair = sdp_pair;
    return seed_info_boundary;
    }

static void Scheduler_Seed_Boundary_destroy(
            Scheduler_Seed_Boundary *seed_info_boundary){
    Boundary_destroy(seed_info_boundary->boundary);
    g_free(seed_info_boundary);
    return;
    }

static void Scheduler_Seed_Boundary_init_forward(gpointer seed_data){
    register Scheduler_Seed_Boundary *seed_info_boundary = seed_data;
    g_assert(seed_info_boundary);
    seed_info_boundary->curr_row = 0;
    seed_info_boundary->curr_interval = 0;
    seed_info_boundary->curr_pos = 0;
    seed_info_boundary->is_finished = FALSE;
    return;
    }

static void Scheduler_Seed_Boundary_next_forward(gpointer seed_data){
    register Scheduler_Seed_Boundary *seed_info_boundary = seed_data;
    register Boundary_Row *boundary_row
        = seed_info_boundary->boundary->row_list->pdata
         [seed_info_boundary->curr_row];
    register Boundary_Interval *interval
        = boundary_row->interval_list->pdata
         [seed_info_boundary->curr_interval];
    /* If more seeds in interval, move curr_pos */
    if(seed_info_boundary->curr_pos < (interval->length-1)){
        seed_info_boundary->curr_pos++;
        return;
        }
    /* If more intervals in row, move curr_interval */
    if(seed_info_boundary->curr_interval
    <  (boundary_row->interval_list->len-1)){
        seed_info_boundary->curr_interval++;
        seed_info_boundary->curr_pos = 0;
        return;
        }
    /* If more rows in boundary, move curr_row */
    if(seed_info_boundary->curr_row
    <  (seed_info_boundary->boundary->row_list->len-1)){
        seed_info_boundary->curr_row++;
        seed_info_boundary->curr_interval = 0;
        seed_info_boundary->curr_pos = 0;
        return;
        }
    seed_info_boundary->is_finished = TRUE;
    return;
    }

static gboolean Scheduler_Seed_Boundary_get_forward(gpointer seed_data,
                                          Scheduler_Seed *seed_info){
    register Scheduler_Seed_Boundary *seed_info_boundary = seed_data;
    register Boundary_Row *boundary_row
        = seed_info_boundary->boundary->row_list->pdata
         [seed_info_boundary->curr_row];
    /* FIXME: fails when boundary is empty */
    register Boundary_Interval *interval
        = boundary_row->interval_list->pdata
         [seed_info_boundary->curr_interval];
    register SDP_Seed *seed;
    if(seed_info_boundary->is_finished)
        return FALSE;
    g_assert(interval->seed_id < seed_info_boundary->sdp_pair->seed_list->len);
    seed = seed_info_boundary->sdp_pair->seed_list->pdata[interval->seed_id];
    seed_info->query_pos = interval->query_pos
                         + seed_info_boundary->curr_pos;
    seed_info->target_pos = boundary_row->target_pos;
    seed_info->seed_id = interval->seed_id;
    seed_info->start_score = 0;
    return TRUE;
    }

static void Scheduler_Seed_Boundary_end_func(gpointer seed_data,
                                           gint seed_id, C4_Score score,
                                           gint end_query_pos,
                                           gint end_target_pos,
                                           STraceback_Cell *stcell){
    register Scheduler_Seed_Boundary *seed_info_boundary = seed_data;
    register SDP_Seed *seed
        = seed_info_boundary->sdp_pair->seed_list->pdata[seed_id];
    /* Is best end for this seed */
    if(seed->max_end->score < score){
        seed->max_end->score = score;
        seed->max_end->query_pos = end_query_pos;
        seed->max_end->target_pos = end_target_pos;
        g_assert(stcell);
        seed->max_end->cell = STraceback_Cell_share(stcell);
        /* FIXME: free old seed max_end->cell ?? */
        }
    return;
    }

/**/

SDP *SDP_create(C4_Model *model){
    register SDP *sdp = g_new0(SDP, 1);
    register C4_Portal *portal;
    /**/
    sdp->thread_ref = ThreadRef_create();
    sdp->sas = SDP_ArgumentSet_create(NULL);
    g_assert(model);
    g_assert(!model->is_open);
    g_assert(Lookahead_Mask_WIDTH > model->max_query_advance);
    g_assert(Lookahead_Mask_WIDTH > model->max_target_advance);
    sdp->model = C4_Model_share(model);
    /* Perform bidirectional SDP only
     * when there is a single match state and no shadows or spans
     */
    sdp->use_boundary = TRUE;
    if((model->shadow_list->len == 0) && (model->span_list->len == 0)){
        if(model->portal_list->len == 1){
            portal = model->portal_list->pdata[0];
            g_assert(portal);
            if(portal->transition_list->len == 1)
                sdp->use_boundary = FALSE;
            }
        }
    /**/
    if(sdp->use_boundary){
        sdp->find_starts_scheduler = Scheduler_create(model,
                              FALSE, FALSE, TRUE,
                              Scheduler_Seed_List_init_reverse,
                              Scheduler_Seed_List_next_reverse,
                              Scheduler_Seed_List_get_reverse,
                              NULL, NULL,
                              sdp->sas->dropoff);
        sdp->find_ends_scheduler = Scheduler_create(model,
                                 TRUE, TRUE, TRUE,
                                 Scheduler_Seed_Boundary_init_forward,
                                 Scheduler_Seed_Boundary_next_forward,
                                 Scheduler_Seed_Boundary_get_forward,
                                 NULL, Scheduler_Seed_Boundary_end_func,
                                 sdp->sas->dropoff);
    } else {
        sdp->find_starts_scheduler = Scheduler_create(model,
                              FALSE, TRUE, FALSE,
                              Scheduler_Seed_List_init_reverse,
                              Scheduler_Seed_List_next_reverse,
                              Scheduler_Seed_List_get_reverse,
                              Scheduler_Seed_List_start_func, NULL,
                              sdp->sas->dropoff);
        sdp->find_ends_scheduler = Scheduler_create(model,
                                 TRUE, TRUE, FALSE,
                                 Scheduler_Seed_List_init_forward,
                                 Scheduler_Seed_List_next_forward,
                                 Scheduler_Seed_List_get_forward,
                                 NULL, Scheduler_Seed_List_end_func,
                                 sdp->sas->dropoff);
        }
    return sdp;
    }

SDP *SDP_share(SDP *sdp){
    g_assert(sdp);
    ThreadRef_share(sdp->thread_ref);
    return sdp;
    }

void SDP_destroy(SDP *sdp){
    g_assert(sdp);
    if(ThreadRef_destroy(sdp->thread_ref))
        return;
    if(sdp->find_starts_scheduler)
        Scheduler_destroy(sdp->find_starts_scheduler);
    if(sdp->find_ends_scheduler)
        Scheduler_destroy(sdp->find_ends_scheduler);
    C4_Model_destroy(sdp->model);
    g_free(sdp);
    return;
    }

/**/

static void SDP_add_codegen(Scheduler *scheduler,
                            GPtrArray *codegen_list){
    register Codegen *codegen;
    if(scheduler){
        codegen = Scheduler_make_Codegen(scheduler);
        g_ptr_array_add(codegen_list, codegen);
        }
    return;
    }

GPtrArray *SDP_get_codegen_list(SDP *sdp){
    register GPtrArray *codegen_list = g_ptr_array_new();
    SDP_add_codegen(sdp->find_starts_scheduler, codegen_list);
    SDP_add_codegen(sdp->find_ends_scheduler, codegen_list);
    g_assert(codegen_list->len);
    return codegen_list;
    }

/**/

static SDP_Terminal *SDP_Terminal_create(void){
    register SDP_Terminal *terminal = g_new(SDP_Terminal, 1);
    terminal->query_pos = 0;
    terminal->target_pos = 0;
    terminal->score = C4_IMPOSSIBLY_LOW_SCORE;
    terminal->cell = NULL;
    return terminal;
    }

static void SDP_Terminal_destroy(SDP_Terminal *terminal,
                                 STraceback *straceback){
    if(terminal->cell)
        STraceback_Cell_destroy(terminal->cell, straceback);
    g_free(terminal);
    return;
    }

/**/

static SDP_Seed *SDP_Seed_create(HSP *hsp, gint id){
    register SDP_Seed *seed = g_new(SDP_Seed, 1);
    seed->seed_id = id;
    seed->hsp = hsp;
    seed->max_start = SDP_Terminal_create();
    seed->max_end = SDP_Terminal_create();
    seed->pq_node = NULL;
    return seed;
    }
/* FIXME: optimisation: use RecycleBins for SDP_Terminal allocations
 */

static void SDP_Seed_destroy(SDP_Seed *seed, STraceback *fwd_straceback,
                                             STraceback *rev_straceback){
    SDP_Terminal_destroy(seed->max_start, rev_straceback);
    SDP_Terminal_destroy(seed->max_end, fwd_straceback);
    g_free(seed);
    return;
    }

/**/

static int SDP_compare_HSP_cobs_in_forward_dp_order(const void *a,
                                                    const void *b){
    register HSP **hsp_a = (HSP**)a,
                 **hsp_b = (HSP**)b;
    register gint target_diff = HSP_target_cobs(*hsp_a)
                              - HSP_target_cobs(*hsp_b);
    if(!target_diff)
        return HSP_query_cobs(*hsp_a)
             - HSP_query_cobs(*hsp_b);
    return target_diff;
    }

static GPtrArray *SDP_Pair_create_seed_list(Comparison *comparison){
    register gint i;
    register SDP_Seed *seed;
    register HSP *hsp, *prev_hsp = NULL;
    register GPtrArray *hsp_list = g_ptr_array_new(),
                       *seed_list = g_ptr_array_new();
    g_assert(Comparison_has_hsps(comparison));
    /* Build a combined HSP list */
    if(comparison->dna_hspset)
        for(i = 0; i < comparison->dna_hspset->hsp_list->len; i++){
            hsp = comparison->dna_hspset->hsp_list->pdata[i];
            g_ptr_array_add(hsp_list, hsp);
            }
    if(comparison->protein_hspset){
        for(i = 0; i < comparison->protein_hspset->hsp_list->len; i++){
            hsp = comparison->protein_hspset->hsp_list->pdata[i];
            g_ptr_array_add(hsp_list, hsp);
            }
        }
    if(comparison->codon_hspset)
        for(i = 0; i < comparison->codon_hspset->hsp_list->len; i++){
            hsp = comparison->codon_hspset->hsp_list->pdata[i];
            g_ptr_array_add(hsp_list, hsp);
            }
    /* Sort hsps on cobs point in DP order */
    qsort(hsp_list->pdata,
          hsp_list->len, sizeof(gpointer),
          SDP_compare_HSP_cobs_in_forward_dp_order);
    /* Make a seed for each unique HSP */
    g_assert(hsp_list->len);
    for(i = 0; i < hsp_list->len; i++){
        hsp = hsp_list->pdata[i];
        if((!prev_hsp)
        || (HSP_query_cobs(hsp) != HSP_query_cobs(prev_hsp))
        || (HSP_target_cobs(hsp) != HSP_target_cobs(prev_hsp))){
            seed = SDP_Seed_create(hsp, seed_list->len);
            g_ptr_array_add(seed_list, seed);
            }
        prev_hsp = hsp;
        }
    g_ptr_array_free(hsp_list, TRUE);
    g_assert(seed_list->len);
    return seed_list;
    }

SDP_Pair *SDP_Pair_create(SDP *sdp, SubOpt *subopt,
                          Comparison *comparison,
                          gpointer user_data){
    register SDP_Pair *sdp_pair = g_new(SDP_Pair, 1);
    g_assert(sdp);
    sdp_pair->sdp = SDP_share(sdp);
    g_assert(comparison);
    g_assert(Comparison_has_hsps(comparison));
    sdp_pair->comparison = Comparison_share(comparison);
    sdp_pair->user_data = user_data;
    sdp_pair->alignment_count = 0;
    sdp_pair->subopt = SubOpt_share(subopt);
    sdp_pair->seed_list = SDP_Pair_create_seed_list(comparison);
    sdp_pair->seed_list_by_score = NULL;
    sdp_pair->boundary = NULL;
    sdp_pair->last_score = C4_IMPOSSIBLY_LOW_SCORE;
    sdp_pair->fwd_straceback = STraceback_create(sdp->model, TRUE);
    sdp_pair->rev_straceback = STraceback_create(sdp->model, FALSE);
    return sdp_pair;
    }

void SDP_Pair_destroy(SDP_Pair *sdp_pair){
    register gint i;
    register SDP_Seed *seed;
    g_assert(sdp_pair);
    Comparison_destroy(sdp_pair->comparison);
    SDP_destroy(sdp_pair->sdp);
    for(i = 0; i < sdp_pair->seed_list->len; i++){
        seed = sdp_pair->seed_list->pdata[i];
        SDP_Seed_destroy(seed, sdp_pair->fwd_straceback,
                               sdp_pair->rev_straceback);
        }
    g_ptr_array_free(sdp_pair->seed_list, TRUE);
    if(sdp_pair->seed_list_by_score)
        g_free(sdp_pair->seed_list_by_score);
    SubOpt_destroy(sdp_pair->subopt);
    if(sdp_pair->boundary)
        Boundary_destroy(sdp_pair->boundary);
    STraceback_destroy(sdp_pair->fwd_straceback);
    STraceback_destroy(sdp_pair->rev_straceback);
    g_free(sdp_pair);
    return;
    }

/**/

static Boundary *SDP_Pair_find_start_points(SDP_Pair *sdp_pair){
    register Scheduler_Seed_List *seed_info_list
           = Scheduler_Seed_List_create(sdp_pair, FALSE);
    register Boundary *boundary = NULL;
    register Scheduler_Pair *spair;
    if(sdp_pair->sdp->use_boundary)
        boundary = Boundary_create();
    spair = Scheduler_Pair_create(sdp_pair->sdp->find_starts_scheduler,
                              sdp_pair->rev_straceback,
                              sdp_pair->comparison->query->len,
                              sdp_pair->comparison->target->len,
                              sdp_pair->subopt, boundary, -1,
                              seed_info_list, sdp_pair->user_data);
    Scheduler_Pair_calculate(spair);
    Scheduler_Pair_destroy(spair);
    Scheduler_Seed_List_destroy(seed_info_list);
    if(boundary)
        Boundary_reverse(boundary);
    return boundary;
    }

static void SDP_Pair_find_end_points(SDP_Pair *sdp_pair){
    register Scheduler_Pair *spair;
    register Scheduler_Seed_List *seed_info_list = NULL;
    register Scheduler_Seed_Boundary *seed_info_boundary = NULL;
    if(sdp_pair->boundary){
        g_assert(sdp_pair->boundary->row_list->len);
        seed_info_boundary = Scheduler_Seed_Boundary_create(sdp_pair->boundary,
                                                            sdp_pair);
        spair = Scheduler_Pair_create(sdp_pair->sdp->find_ends_scheduler,
                          sdp_pair->fwd_straceback,
                          sdp_pair->comparison->query->len,
                          sdp_pair->comparison->target->len,
                          sdp_pair->subopt,
                          sdp_pair->boundary, -1, seed_info_boundary,
                          sdp_pair->user_data);
    } else {
        seed_info_list = Scheduler_Seed_List_create(sdp_pair, TRUE);
        spair = Scheduler_Pair_create(sdp_pair->sdp->find_ends_scheduler,
                          sdp_pair->fwd_straceback,
                          sdp_pair->comparison->query->len,
                          sdp_pair->comparison->target->len,
                          sdp_pair->subopt,
                          sdp_pair->boundary, -1, seed_info_list,
                          sdp_pair->user_data);
        }
    Scheduler_Pair_calculate(spair);
    Scheduler_Pair_destroy(spair);
    if(sdp_pair->boundary){
        g_assert(seed_info_boundary);
        Scheduler_Seed_Boundary_destroy(seed_info_boundary);
    } else {
        g_assert(seed_info_list);
        Scheduler_Seed_List_destroy(seed_info_list);
        }
    return;
    }

/**/

static Boundary *SDP_Pair_update_starts(SDP_Pair *sdp_pair){
    register gint i;
    register SDP_Seed *seed;
    g_assert(sdp_pair->seed_list->len);
    /* Set max start scores low */
    for(i = 0; i < sdp_pair->seed_list->len; i++){
        seed = sdp_pair->seed_list->pdata[i];
        seed->max_start->score = C4_IMPOSSIBLY_LOW_SCORE;
        if(!sdp_pair->sdp->use_boundary){
            g_assert(seed->max_start->cell);
            STraceback_Cell_destroy(seed->max_start->cell,
                                    sdp_pair->rev_straceback);
            seed->max_start->cell = NULL;
            }
        }
    return SDP_Pair_find_start_points(sdp_pair);
    }

static void SDP_Pair_update_ends(SDP_Pair *sdp_pair){
    register gint i;
    register SDP_Seed *seed;
    g_assert(sdp_pair->seed_list->len);
    /* Set max end scores low */
    for(i = 0; i < sdp_pair->seed_list->len; i++){
        seed = sdp_pair->seed_list->pdata[i];
        seed->max_end->score = C4_IMPOSSIBLY_LOW_SCORE;
        if(seed->max_end->cell){
            STraceback_Cell_destroy(seed->max_end->cell,
                                    sdp_pair->fwd_straceback);
            seed->max_end->cell = NULL;
            }
        }
    SDP_Pair_find_end_points(sdp_pair);
    return;
    }

/**/

static void SDP_Pair_add_traceback(SDP_Pair *sdp_pair, SDP_Seed *best_seed,
                                   gboolean is_forward,
                                   Alignment *alignment){
    register STraceback_List *stlist = NULL;
    register gint i;
    register STraceback_Operation *operation;
    register AlignmentOperation *last;
    register STraceback_Operation *first;
    g_assert(best_seed);
    if(is_forward){
        g_assert(best_seed->max_end);
        g_assert(best_seed->max_end->cell);
        g_assert(best_seed->max_end->cell->transition->output
              == sdp_pair->sdp->model->end_state->state);
        stlist = STraceback_List_create(sdp_pair->fwd_straceback,
                                        best_seed->max_end->cell);
        if(!sdp_pair->sdp->use_boundary){ /* Check join is valid */
            g_assert(alignment->operation_list->len);
            last = alignment->operation_list->pdata
                  [alignment->operation_list->len-1];
            g_assert(stlist->operation_list->len);
            first = stlist->operation_list->pdata[0];
            g_assert(last->transition->output == first->transition->output);
            }
        /* Include 1st operation only when using boundary */
        for(i = sdp_pair->sdp->use_boundary?0:1;
            i < stlist->operation_list->len; i++){
            operation = stlist->operation_list->pdata[i];
            Alignment_add(alignment, operation->transition,
                          operation->length);
            }
    } else {
        g_assert(best_seed->max_start->cell);
        g_assert(best_seed->max_start->cell->transition->input
              == sdp_pair->sdp->model->start_state->state);
        stlist = STraceback_List_create(sdp_pair->rev_straceback,
                                        best_seed->max_start->cell);
        /* Don't include last operation (to end state) */
        for(i = stlist->operation_list->len-1; i >= 1; i--){
            operation = stlist->operation_list->pdata[i];
            Alignment_add(alignment, operation->transition,
                          operation->length);
            }
        }
    STraceback_List_destroy(stlist);
    return;
    }

static void SDP_Seed_find_start(SDP_Seed *best_seed, C4_Model *model){
    register STraceback_Cell *stcell;
    best_seed->max_start->query_pos
        = best_seed->max_end->query_pos;
    best_seed->max_start->target_pos
        = best_seed->max_end->target_pos;
    stcell = best_seed->max_end->cell;
    g_assert(stcell);
    do {
        best_seed->max_start->query_pos
            -= (stcell->transition->advance_query * stcell->length);
        best_seed->max_start->target_pos
            -= (stcell->transition->advance_target * stcell->length);
        stcell = stcell->prev;
        g_assert(stcell);
    } while(stcell->transition->input != model->start_state->state);
    return;
    }

static Alignment *SDP_Pair_find_path(SDP_Pair *sdp_pair,
                                     SDP_Seed *best_seed,
                                     Boundary *boundary){
    register Region *alignment_region;
    register Alignment *alignment;
    /**/
    if(sdp_pair->sdp->use_boundary)
        SDP_Seed_find_start(best_seed, sdp_pair->sdp->model);
    alignment_region = Region_create(best_seed->max_start->query_pos,
                           best_seed->max_start->target_pos,
                           best_seed->max_end->query_pos
                          -best_seed->max_start->query_pos,
                           best_seed->max_end->target_pos
                          -best_seed->max_start->target_pos);
    alignment = Alignment_create(sdp_pair->sdp->model,
                         alignment_region, best_seed->max_end->score);
    if(sdp_pair->sdp->use_boundary){
        SDP_Pair_add_traceback(sdp_pair, best_seed, TRUE, alignment);
    } else {
        SDP_Pair_add_traceback(sdp_pair, best_seed, FALSE, alignment);
        SDP_Pair_add_traceback(sdp_pair, best_seed, TRUE, alignment);
        }
    g_assert(Alignment_is_valid(alignment, alignment_region,
                                sdp_pair->user_data));
    Region_destroy(alignment_region);
    return alignment;
    }

static int SDP_compare_SDP_Seed_by_score(const void *a,
                                         const void *b){
    register SDP_Seed **seed_a = (SDP_Seed**)a,
                      **seed_b = (SDP_Seed**)b;
    return (*seed_b)->max_end->score
         - (*seed_a)->max_end->score;
    }

Alignment *SDP_Pair_next_path(SDP_Pair *sdp_pair, C4_Score threshold){
    register SDP_Seed *seed, *best_seed = NULL;
    register Alignment *alignment = NULL;
    register gint i;
    /* Make the start and end points up to date */
    if(sdp_pair->alignment_count){
        if(!sdp_pair->sdp->sas->single_pass_subopt){ /* multipass */
            if(sdp_pair->boundary)
                Boundary_destroy(sdp_pair->boundary);
            sdp_pair->boundary = SDP_Pair_update_starts(sdp_pair);
            SDP_Pair_update_ends(sdp_pair);
            }
    } else {
        g_assert(!sdp_pair->boundary);
        sdp_pair->boundary = SDP_Pair_find_start_points(sdp_pair);
        /* Boundary_print_gnuplot(sdp_pair->boundary, 1); */
        SDP_Pair_find_end_points(sdp_pair);
        if(sdp_pair->sdp->sas->single_pass_subopt){
            sdp_pair->seed_list_by_score = g_new(gpointer,
                                                 sdp_pair->seed_list->len);
            for(i = 0; i < sdp_pair->seed_list->len; i++)
                sdp_pair->seed_list_by_score[i]
                    = sdp_pair->seed_list->pdata[i];
            qsort(sdp_pair->seed_list_by_score,
                  sdp_pair->seed_list->len, sizeof(gpointer),
                  SDP_compare_SDP_Seed_by_score);
            sdp_pair->single_pass_pos = 0;
            }
        }
/**/
    if(sdp_pair->sdp->sas->single_pass_subopt){ /* singlepass */
        best_seed = NULL;
        while(sdp_pair->single_pass_pos < sdp_pair->seed_list->len){
            best_seed = sdp_pair->seed_list_by_score
                       [sdp_pair->single_pass_pos++];
            if(best_seed->max_end->score < threshold)
                return NULL;
            alignment = SDP_Pair_find_path(sdp_pair, best_seed,
                                           sdp_pair->boundary);
            if(SubOpt_overlaps_alignment(sdp_pair->subopt, alignment)){
                Alignment_destroy(alignment);
                best_seed = NULL;
            } else {
                break;
                }
            }
        if(!best_seed)
            return NULL;
        /**/
    } else { /* multipass */
        g_assert(sdp_pair->seed_list->len);
        best_seed = sdp_pair->seed_list->pdata[0];
        for(i = 1; i < sdp_pair->seed_list->len; i++){
            seed = sdp_pair->seed_list->pdata[i];
            if(best_seed->max_end->score < seed->max_end->score)
                best_seed = seed;
            }
        if(best_seed->max_end->score < threshold)
            return NULL; /* Score below threshold */
        alignment = SDP_Pair_find_path(sdp_pair, best_seed,
                                       sdp_pair->boundary);
        }
    g_assert(best_seed);
    g_assert(alignment);
    /* Check score is less than previous */
    g_assert((sdp_pair->last_score < 0)
          || (best_seed->max_end->score <= sdp_pair->last_score));
    sdp_pair->alignment_count++;
    sdp_pair->last_score = best_seed->max_end->score;
    best_seed->max_end->score = C4_IMPOSSIBLY_LOW_SCORE;
    return alignment;
    }

/**/

