/****************************************************************\
*                                                                *
*  C4 dynamic programming library - code for heuristic pairs     *
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

#include "optimal.h"
#include "hpair.h"
#include "sar.h"
#include "rangetree.h"

/**/

typedef struct {
    Heuristic_Match *match;
               gint  hsp_id;  /* The id of the hsp in the portal */
       SAR_Terminal *sar_start;
       SAR_Terminal *sar_end;
} HPair_NodeData;

#define HPair_NodeData_get_hsp_set(hpair, node_data) \
    ((HSPset*)                                 \
    ((hpair)->portal_data_list->pdata[(node_data)->match->portal->id]))

#define HPair_NodeData_get_hsp(hpair, node_data) \
    ((HSP*)(HPair_NodeData_get_hsp_set(hpair, node_data) \
    ->hsp_list->pdata[(node_data)->hsp_id]))

static HPair_NodeData *HPair_NodeData_create(HPair *hpair,
                             Heuristic_Match *match, gint hsp_id){
    register HPair_NodeData *hnd = g_new(HPair_NodeData, 1);
    register HSP *hsp;
    hnd->match = match;
    hnd->hsp_id = hsp_id;
    hsp = HPair_NodeData_get_hsp(hpair, hnd);
    hnd->sar_start = SAR_Terminal_create(hsp, hpair, match, TRUE);
    hnd->sar_end   = SAR_Terminal_create(hsp, hpair, match, FALSE);
    return hnd;
    }

static void HPair_NodeData_destroy(HPair_NodeData *hnd){
    if(hnd->sar_start)
        SAR_Terminal_destroy(hnd->sar_start);
    if(hnd->sar_end)
        SAR_Terminal_destroy(hnd->sar_end);
    g_free(hnd);
    return;
    }

/**/

typedef struct {
    SAR_Join *sar_join; /* NULL if span edge */
    SAR_Span *sar_span; /* NULL if join edge */
} HPair_EdgeData;

static HPair_EdgeData *HPair_EdgeData_create(SAR_Join *sar_join,
                                             SAR_Span *sar_span){
    register HPair_EdgeData *hed = g_new(HPair_EdgeData, 1);
    g_assert(!(sar_join && sar_span));
    g_assert(sar_join || sar_span);
    hed->sar_join = sar_join;
    hed->sar_span = sar_span;
    return hed;
    }

static void HPair_EdgeData_destroy(HPair_EdgeData *hed){
    if(hed->sar_join)
        SAR_Join_destroy(hed->sar_join);
    if(hed->sar_span)
        SAR_Span_destroy(hed->sar_span);
    g_free(hed);
    return;
    }

/**/

static gboolean HPair_SubOpt_check_diagonal(gint query_pos,
                                            gint target_pos,
                                            gint path_id,
                                            gpointer user_data){
    register HSP *hsp = user_data;
    register gint diagonal = (target_pos*HSP_query_advance(hsp))
                           - (query_pos*HSP_target_advance(hsp));
    if(diagonal == HSP_diagonal(hsp))
        return TRUE;
    return FALSE;
    }
/* FIXME: optimisation: precompute HSP_diagonal() */

static gboolean HPair_check_entry(HPair *hpair, HSP *hsp,
                                  Region *region){
    Region search_region;
    Region_init_static(&search_region,
                       HSP_query_cobs(hsp),
                       HSP_target_cobs(hsp),
                       region->query_start-HSP_query_cobs(hsp),
                       region->target_start-HSP_target_cobs(hsp));
    return SubOpt_find(hpair->subopt, &search_region,
            HPair_SubOpt_check_diagonal, hsp);
    }
/* Returns true if subopt clashes with region -> HSP cobs */

static gboolean HPair_check_exit(HPair *hpair, HSP *hsp,
                                 Region *region){
    Region search_region;
    Region_init_static(&search_region,
          Region_query_end(region),
          Region_target_end(region),
          HSP_query_cobs(hsp)-Region_query_end(region),
          HSP_target_cobs(hsp)-Region_target_end(region));
    return SubOpt_find(hpair->subopt, &search_region,
            HPair_SubOpt_check_diagonal, hsp);
    }

static gboolean HPair_SubOpt_check_region_since(gint query_pos,
                                                gint target_pos,
                                                gint path_id,
                                                gpointer user_data){
    register gint last_updated = GPOINTER_TO_INT(user_data);
    if(path_id >= last_updated)
        return TRUE;
    return FALSE;
    }

static gboolean HPair_check_region_since(HPair *hpair, Region *region,
                                         gint last_updated){
    if(SubOpt_find(hpair->subopt, region,
                   HPair_SubOpt_check_region_since,
                   GINT_TO_POINTER(last_updated)))
        return TRUE;
    return FALSE;
    }

/**/

static C4_Score HPair_confirm_edge_func(gpointer src_data,
                                        gpointer edge_data,
                                        gpointer dst_data,
                                        gpointer user_data){
    register HPair *hpair = user_data;
    register HPair_EdgeData *hpair_edge_data = edge_data;
    register C4_Score score;
    register HPair_NodeData *src_node_data = src_data,
                            *dst_node_data = dst_data;
    register HSP
        *src_hsp = HPair_NodeData_get_hsp(hpair, src_node_data),
        *dst_hsp = HPair_NodeData_get_hsp(hpair, dst_node_data);
    if(hpair_edge_data->sar_join){
        g_assert(!hpair_edge_data->sar_span);
        if((HPair_check_entry(hpair, src_hsp,
                               hpair_edge_data->sar_join->region))
        || (HPair_check_exit(hpair, dst_hsp,
                             hpair_edge_data->sar_join->region)))
            return C4_IMPOSSIBLY_LOW_SCORE;
        score = SAR_Join_find_score(hpair_edge_data->sar_join,
                                    hpair);
    } else {
        if((HPair_check_entry(hpair, src_hsp,
                               hpair_edge_data->sar_span->src_region))
        || (HPair_check_exit(hpair, dst_hsp,
                             hpair_edge_data->sar_span->dst_region)))
            return C4_IMPOSSIBLY_LOW_SCORE;
        g_assert(hpair_edge_data->sar_span);
        score = SAR_Span_find_score(hpair_edge_data->sar_span, hpair);
        }
    return score;
    }

static C4_Score HPair_update_edge_func(gpointer src_data,
                                       gpointer edge_data,
                                       gpointer dst_data,
                                       gpointer user_data,
                                       C4_Score prev_score,
                                       gint last_updated){
    register HPair *hpair = user_data;
    register HPair_EdgeData *hpair_edge_data = edge_data;
    register HPair_NodeData *src_node_data = src_data,
                            *dst_node_data = dst_data;
    register HSP
        *src_hsp = HPair_NodeData_get_hsp(hpair, src_node_data),
        *dst_hsp = HPair_NodeData_get_hsp(hpair, dst_node_data);
    if(hpair_edge_data->sar_join){
        g_assert(!hpair_edge_data->sar_span);
        if((HPair_check_entry(hpair, src_hsp,
                               hpair_edge_data->sar_join->region))
        || (HPair_check_exit(hpair, dst_hsp,
                             hpair_edge_data->sar_join->region)))
            return C4_IMPOSSIBLY_LOW_SCORE;
        if(HPair_check_region_since(hpair,
                                     hpair_edge_data->sar_join->region,
                                     last_updated)){
            return SAR_Join_find_score(hpair_edge_data->sar_join,
                                       hpair);
            }
    } else {
        g_assert(hpair_edge_data->sar_span);
        if((HPair_check_entry(hpair, src_hsp,
                               hpair_edge_data->sar_span->src_region))
        || (HPair_check_exit(hpair, dst_hsp,
                             hpair_edge_data->sar_span->dst_region)))
            return C4_IMPOSSIBLY_LOW_SCORE;
        if((HPair_check_region_since(hpair,
                        hpair_edge_data->sar_span->src_region,
                        last_updated))
        || (HPair_check_region_since(hpair,
                        hpair_edge_data->sar_span->dst_region,
                        last_updated)))
            return SAR_Span_find_score(hpair_edge_data->sar_span,
                                       hpair);
        }
    return prev_score;
    }

static C4_Score HPair_confirm_start_func(gpointer node_data,
                                         gpointer user_data){
    register HPair *hpair = user_data;
    register HPair_NodeData *hpair_node_data = node_data;
    register HSP *hsp = HPair_NodeData_get_hsp(hpair, hpair_node_data);
    g_assert(hpair_node_data->sar_start);
    if(HPair_check_exit(hpair, hsp,
                        hpair_node_data->sar_start->region))
        return C4_IMPOSSIBLY_LOW_SCORE;
    return SAR_Terminal_find_score(hpair_node_data->sar_start,
               hpair_node_data->match->start_terminal->optimal,
               hpair);
    }

static C4_Score HPair_update_start_func(gpointer node_data,
                                        gpointer user_data,
                                        C4_Score prev_score,
                                        gint last_updated){
    register HPair *hpair = user_data;
    register HPair_NodeData *hpair_node_data = node_data;
    register HSP *hsp = HPair_NodeData_get_hsp(hpair, hpair_node_data);
    g_assert(hpair_node_data->sar_start);
    if(HPair_check_exit(hpair, hsp,
                        hpair_node_data->sar_start->region))
        return C4_IMPOSSIBLY_LOW_SCORE;
    if(HPair_check_region_since(hpair,
                                hpair_node_data->sar_start->region,
                                last_updated))
        return SAR_Terminal_find_score(hpair_node_data->sar_start,
                hpair_node_data->match->start_terminal->optimal,
                hpair);
    return prev_score;
    }

static C4_Score HPair_confirm_end_func(gpointer node_data,
                                       gpointer user_data){
    register HPair *hpair = user_data;
    register HPair_NodeData *hpair_node_data = node_data;
    register HSP *hsp = HPair_NodeData_get_hsp(hpair, hpair_node_data);
    if(HPair_check_entry(hpair, hsp,
                          hpair_node_data->sar_end->region))
        return C4_IMPOSSIBLY_LOW_SCORE;
    g_assert(hpair_node_data->sar_end);
    return SAR_Terminal_find_score(hpair_node_data->sar_end,
               hpair_node_data->match->end_terminal->optimal,
               hpair);
    }

static C4_Score HPair_update_end_func(gpointer node_data,
                                      gpointer user_data,
                                      C4_Score prev_score,
                                      gint last_updated){
    register HPair *hpair = user_data;
    register HPair_NodeData *hpair_node_data = node_data;
    register HSP *hsp = HPair_NodeData_get_hsp(hpair, hpair_node_data);
    g_assert(hpair_node_data->sar_end);
    if(HPair_check_entry(hpair, hsp,
                          hpair_node_data->sar_end->region))
        return C4_IMPOSSIBLY_LOW_SCORE;
    if(HPair_check_region_since(hpair,
                                 hpair_node_data->sar_end->region,
                                 last_updated))
        return SAR_Terminal_find_score(hpair_node_data->sar_end,
                hpair_node_data->match->end_terminal->optimal,
                hpair);
    return prev_score;
    }

/**/

static void HPair_destroy_node_data_func(gpointer data){
    register HPair_NodeData *hnd = data;
    HPair_NodeData_destroy(hnd);
    return;
    }

static void HPair_destroy_edge_data_func(gpointer data){
    register HPair_EdgeData *hed = data;
    HPair_EdgeData_destroy(hed);
    return;
    }

/**/

HPair *HPair_create(Heuristic *heuristic, SubOpt *subopt,
                    gint query_length, gint target_length,
                    gint verbosity, gpointer user_data){
    register HPair *hpair = g_new(HPair, 1);
    g_assert(heuristic);
    g_assert(query_length >= 0);
    g_assert(target_length >= 0);
    hpair->heuristic = Heuristic_share(heuristic);
    hpair->query_length = query_length;
    hpair->target_length = target_length;
    hpair->verbosity = verbosity;
    hpair->user_data = user_data;
    hpair->is_finalised = FALSE;
    /**/
    hpair->portal_data_list = g_ptr_array_new();
    g_assert(heuristic->model->portal_list);
    g_assert(heuristic->model->portal_list->len);
    g_ptr_array_set_size(hpair->portal_data_list,
                         heuristic->model->portal_list->len);
    hpair->bsdp_node_offset = g_new0(gint, heuristic->match_total);
    /* FIXME: merge confirm and update funcs */
    hpair->bsdp = BSDP_create(HPair_confirm_edge_func,
                              HPair_confirm_start_func,
                              HPair_confirm_end_func,
                              HPair_update_edge_func,
                              HPair_update_start_func,
                              HPair_update_end_func,
                              HPair_destroy_node_data_func,
                              HPair_destroy_edge_data_func,
                              hpair);
    hpair->subopt = SubOpt_share(subopt);
    return hpair;
    }

static void HPair_destroy_offset_list(HPair *hpair){
    if(hpair->bsdp_node_offset){
        g_free(hpair->bsdp_node_offset);
        hpair->bsdp_node_offset = NULL;
        }
    return;
    }

void HPair_destroy(HPair *hpair){
    register gint i;
    register HSPset *hspset;
    for(i = 0; i < hpair->portal_data_list->len; i++){
        hspset= hpair->portal_data_list->pdata[i];
        if(hspset)
            HSPset_destroy(hspset);
        }
    g_ptr_array_free(hpair->portal_data_list, TRUE);
    HPair_destroy_offset_list(hpair);
    BSDP_destroy(hpair->bsdp);
    Heuristic_destroy(hpair->heuristic);
    SubOpt_destroy(hpair->subopt);
    g_free(hpair);
    return;
    }

/**/

void HPair_add_hspset(HPair *hpair, C4_Portal *portal, HSPset *hsp_set){
    g_assert(hpair);
    g_assert(!hpair->is_finalised);
    g_assert(hsp_set);
    /* Check this is the first hsp_set for this portal */
    g_assert(!hpair->portal_data_list->pdata[portal->id]);
    hpair->portal_data_list->pdata[portal->id] = HSPset_share(hsp_set);
    return;
    }

/**/

static void HPair_initialise_bsdp_nodes(HPair *hpair){
    register gint i, j;
    register Heuristic_Match *match;
    register HSPset *hsp_set;
    register HSP *hsp;
    register gboolean start_ok, end_ok;
    register C4_Score start_bound, end_bound;
    register HPair_NodeData *hpair_node_data;
    register gint node_id, node_total = 0;
    for(i = 0; i < hpair->heuristic->match_total; i++){
        match = hpair->heuristic->match_list[i];
        g_assert(match);
        hsp_set = hpair->portal_data_list->pdata[match->portal->id];
        if(!hsp_set)
            continue;
        for(j = 0; j < hsp_set->hsp_list->len; j++){
            hsp = hsp_set->hsp_list->pdata[j];
            hpair_node_data = HPair_NodeData_create(hpair, match, j);
            if(hpair_node_data->sar_start){
                start_bound
                 = SAR_Terminal_find_bound(hpair_node_data->sar_start,
                                     match->start_terminal->bound);
                start_ok = TRUE;
            } else {
                start_bound = C4_IMPOSSIBLY_LOW_SCORE;
                start_ok = FALSE;
                }
            if(hpair_node_data->sar_end){
                end_bound
                 = SAR_Terminal_find_bound(hpair_node_data->sar_end,
                                     match->end_terminal->bound);
                end_ok = TRUE;
            } else {
                end_bound = C4_IMPOSSIBLY_LOW_SCORE;
                end_ok = FALSE;
                }
            node_id = BSDP_add_node(hpair->bsdp,
                                    hpair_node_data, hsp->score,
                                    start_ok, end_ok,
                                    start_bound, end_bound);
            node_total++;
            if(!hpair->bsdp_node_offset[i])
                hpair->bsdp_node_offset[i] = node_id + 1;
            }
        }
    if(hpair->verbosity > 1)
        g_message("Built HPair using [%d] nodes", node_total);
    return;
    }

static gboolean HPair_hsp_pair_is_valid(HSP *src, HSP *dst){
    if(src == dst) /* If the same */
        return FALSE; /* FIXME: not required with test below */
    if((HSP_query_cobs(src) == HSP_query_cobs(dst))
    && (HSP_target_cobs(src) == HSP_target_cobs(dst)))
        return FALSE; /* Can be on same position with mixed HSPsets */
    /* Check above and to the left */
    if(HSP_query_cobs(src) > HSP_query_cobs(dst))
        return FALSE;
    if(HSP_target_cobs(src) > HSP_target_cobs(dst))
        return FALSE;
    return TRUE;
    }

static void HPair_hsp_pair_calc_emit(HPair *hpair, HSP *src, HSP *dst,
                            gint *query_emit, gint *target_emit){
    register gboolean query_overlapping, target_overlapping;
    g_assert(query_emit);
    g_assert(target_emit);
    (*query_emit) = (*target_emit) = 0;
    query_overlapping = (HSP_query_end(src) > dst->query_start);
    target_overlapping = (HSP_target_end(src) > dst->target_start);
    (*query_emit) = dst->query_start - HSP_query_end(src);
    if(query_overlapping)
        (*query_emit) = (*query_emit)
                      % HSP_query_advance(dst);
    (*target_emit) = dst->target_start - HSP_target_end(src);
    if(target_overlapping)
        (*target_emit) = (*target_emit)
                       % HSP_target_advance(dst);
    /**/
    if(query_overlapping && (!target_overlapping))
        (*target_emit) += ((HSP_query_end(src) - dst->query_start)
                          *(HSP_target_advance(dst)/HSP_query_advance(src)));
    if(target_overlapping && (!query_overlapping))
        (*query_emit) += ((HSP_target_end(src) - dst->target_start)
                          *(HSP_query_advance(dst)/HSP_target_advance(src)));
    /**/
    g_assert((*query_emit) >= 0);
    g_assert((*target_emit) >= 0);
    g_assert((*query_emit) <= hpair->query_length);
    g_assert((*target_emit) <= hpair->target_length);
    /**/
    return;
    }
/* The emit distance is the shortest distance
 * from part of the query to part of the target.
 */

static gboolean HPair_Join_is_valid(Heuristic_Join *join,
                                    gint query_emit, gint target_emit){
    if(query_emit > join->bound->query_range)
        return FALSE;
    if(target_emit > join->bound->target_range)
        return FALSE;
    return TRUE;
    }

static gboolean HPair_Span_is_valid(Heuristic_Span *heuristic_span,
                                    HSP *src, HSP *dst,
                                    gint query_emit, gint target_emit){
    register gint effective_query_limit, effective_target_limit;
    effective_query_limit = heuristic_span->span->max_query
                          + heuristic_span->src_bound->query_range
                          + heuristic_span->dst_bound->query_range;
    effective_target_limit = heuristic_span->span->max_target
                           + heuristic_span->src_bound->target_range
                           + heuristic_span->dst_bound->target_range;
    if(query_emit > effective_query_limit)
        return FALSE;
    if(target_emit > effective_target_limit)
        return FALSE;
    if(query_emit < heuristic_span->span->min_query)
        return FALSE;
    if(target_emit < heuristic_span->span->min_target)
        return FALSE;
    return TRUE;
    }

static void HPair_add_candidate_hsp_pair(HPair *hpair,
                              Heuristic_Pair *pair,
                              HSP *src, HSP *dst,
                              gint src_hsp_id, gint dst_hsp_id,
                              C4_Portal *src_portal,
                              C4_Portal *dst_portal,
                              gint *edge_total){
    register SAR_Join *sar_join;
    register SAR_Span *sar_span;
    register C4_Score bound_score;
    register HPair_EdgeData *hpair_edge_data;
    register Heuristic_Span *heuristic_span;
    register gint i, src_node_id, dst_node_id;
    gint query_emit, target_emit;
    if(!HPair_hsp_pair_is_valid(src, dst))
        return;
    src_node_id = hpair->bsdp_node_offset[pair->src->id]
                + src_hsp_id - 1;
    dst_node_id = hpair->bsdp_node_offset[pair->dst->id]
                + dst_hsp_id - 1;
    HPair_hsp_pair_calc_emit(hpair, src, dst,
                             &query_emit, &target_emit);
    if((HPair_Join_is_valid(pair->join,
                            query_emit, target_emit))
    && (sar_join = SAR_Join_create(src, dst, hpair, pair))){
        bound_score = SAR_Join_find_bound(sar_join);
        hpair_edge_data = HPair_EdgeData_create(sar_join, NULL);
        BSDP_add_edge(hpair->bsdp, hpair_edge_data,
            src_node_id, dst_node_id, bound_score);
        (*edge_total)++;
    } else {
        for(i = 0; i < pair->span_list->len; i++){
            heuristic_span = pair->span_list->pdata[i];
            if(HPair_Span_is_valid(heuristic_span, src, dst,
                                   query_emit, target_emit)){
                sar_span = SAR_Span_create(src, dst, hpair,
                        heuristic_span, src_portal, dst_portal);
                if(sar_span){
                    bound_score = SAR_Span_find_bound(sar_span);
                    if(bound_score <= C4_IMPOSSIBLY_LOW_SCORE){
                        SAR_Span_destroy(sar_span);
                        continue;
                        }
                    /* Build edge */
                    hpair_edge_data = HPair_EdgeData_create(NULL,
                                                        sar_span);
                    /* Add edge to BSDP */
                    BSDP_add_edge(hpair->bsdp, hpair_edge_data,
                        src_node_id, dst_node_id, bound_score);
                    (*edge_total)++;
                    }
                }
            }
        }
    return;
    }

typedef struct {
             HPair *hpair;
    Heuristic_Pair *pair;
         C4_Portal *src_portal;
         C4_Portal *dst_portal;
            HSPset *dst_hsp_set;
               HSP *src_hsp;
              gint  src_hsp_id;
              gint *edge_total;
} HPair_RangeTree_Report_Data;

static gboolean HPair_RangeTree_ReportFunc(gint x, gint y,
                                gpointer info, gpointer user_data){
    register HPair_RangeTree_Report_Data *hrtrd = user_data;
    register gint dst_hsp_id = GPOINTER_TO_INT(info);
    HPair_add_candidate_hsp_pair(hrtrd->hpair, hrtrd->pair,
                      hrtrd->src_hsp,
                      hrtrd->dst_hsp_set->hsp_list->pdata[dst_hsp_id],
                      hrtrd->src_hsp_id, dst_hsp_id,
                      hrtrd->src_portal, hrtrd->dst_portal,
                      hrtrd->edge_total);
    return FALSE;
    }

static void HPair_find_candidate_hsp_pairs(HPair *hpair,
                                           Heuristic_Pair *pair,
                                           gint *edge_total){
    register gint src_portal_id = pair->src->portal->id,
                  dst_portal_id = pair->dst->portal->id;
    register HSPset
        *src_hsp_set = hpair->portal_data_list->pdata[src_portal_id],
        *dst_hsp_set = hpair->portal_data_list->pdata[dst_portal_id];
    register C4_Portal
        *src_portal = hpair->heuristic->model->portal_list->pdata
                                                  [src_portal_id],
        *dst_portal = hpair->heuristic->model->portal_list->pdata
                                                  [dst_portal_id];
    register HSP *src, *dst;
    register gint i;
    register RangeTree *rangetree = RangeTree_create();
    register HSP *max_dst_cobs_hsp;
    g_assert(src_portal == pair->src->portal);
    g_assert(dst_portal == pair->dst->portal);
    gint max_query_join_range, max_target_join_range;
    HPair_RangeTree_Report_Data hrtrd;
    g_assert(src_hsp_set);
    g_assert(dst_hsp_set);
    max_dst_cobs_hsp = dst_hsp_set->hsp_list->pdata[0];
    hrtrd.hpair = hpair;
    hrtrd.pair = pair;
    hrtrd.src_portal = src_portal;
    hrtrd.dst_portal = dst_portal;
    hrtrd.dst_hsp_set = dst_hsp_set;
    hrtrd.edge_total = edge_total;
    Heuristic_Pair_get_max_range(pair, &max_query_join_range,
                                       &max_target_join_range);
    /* Build RangeTree from the dst HSPset */
    for(i = 0; i < dst_hsp_set->hsp_list->len; i++){
        dst = dst_hsp_set->hsp_list->pdata[i];
        RangeTree_add(rangetree,
                      HSP_query_cobs(dst), HSP_target_cobs(dst),
                      GINT_TO_POINTER(i));
        if(max_dst_cobs_hsp->cobs < dst->cobs)
           max_dst_cobs_hsp = dst;
        }
    /* Search RangeTree with each HSP from the src HSPset */
    for(i = 0; i < src_hsp_set->hsp_list->len; i++){
        src = src_hsp_set->hsp_list->pdata[i];
        hrtrd.src_hsp = src;
        hrtrd.src_hsp_id = i;
        RangeTree_find(rangetree,
                       HSP_query_cobs(src),
                       (HSP_query_cobs(src) - src->query_start)
                       + (HSP_query_cobs(max_dst_cobs_hsp) -
                           max_dst_cobs_hsp->query_start)
                       + max_query_join_range,
                       HSP_target_cobs(src),
                       (HSP_target_cobs(src) - src->target_start)
                       + (HSP_target_cobs(max_dst_cobs_hsp) -
                           max_dst_cobs_hsp->target_start)
                       + max_target_join_range,
                       HPair_RangeTree_ReportFunc, &hrtrd);
        }
    RangeTree_destroy(rangetree, NULL, NULL);
    return;
    }

static void HPair_initialise_bsdp_edges(HPair *hpair){
    register gint i, j;
    register Heuristic_Pair *pair;
    gint edge_total = 0;
    for(i = 0; i < hpair->heuristic->match_total; i++){
        for(j = 0; j < hpair->heuristic->match_total; j++){
            pair = hpair->heuristic->pair_matrix[i][j];
            if(pair)
                HPair_find_candidate_hsp_pairs(hpair, pair,
                                               &edge_total);
            }
        }
    if(hpair->verbosity > 1)
        g_message("Built HPair using [%d] edges", edge_total);
    return;
    }

void HPair_finalise(HPair *hpair, C4_Score threshold){
    g_assert(hpair);
    g_assert(!hpair->is_finalised);
    HPair_initialise_bsdp_nodes(hpair);
    HPair_initialise_bsdp_edges(hpair);
    HPair_destroy_offset_list(hpair);
    BSDP_initialise(hpair->bsdp, threshold);
    hpair->is_finalised = TRUE;
    return;
    }

/**/

#if 0
static Alignment *HPair_refine_alignment(HPair *hpair,
                                         Alignment *alignment){
    register Alignment *refined_alignment = NULL;
    register Region *region;
    register gint query_region_start, target_region_start;
    g_assert(hpair->heuristic->optimal);
    if(hpair->verbosity > 1)
        g_message("Refining alignment ... (%d)", alignment->score);
    switch(hpair->heuristic->has->refinement){
        case Heuristic_Refinement_FULL:
            region = Region_create(0, 0, hpair->query_length,
                                         hpair->target_length);
            refined_alignment = Optimal_find_path(
                    hpair->heuristic->optimal, region,
                    hpair->user_data, alignment->score, hpair->subopt);
            g_assert(refined_alignment);
            Region_destroy(region);
            break;
        case Heuristic_Refinement_REGION:
            query_region_start
                = MAX(0, alignment->region->query_start
                       - hpair->heuristic->has->refinement_boundary);
            target_region_start
                = MAX(0, alignment->region->target_start
                       - hpair->heuristic->has->refinement_boundary);
            region = Region_create(query_region_start,
                                   target_region_start,
                MIN(hpair->query_length,
                    Region_query_end(alignment->region)
                  + hpair->heuristic->has->refinement_boundary)
                - query_region_start,
                MIN(hpair->target_length,
                    Region_target_end(alignment->region)
                  + hpair->heuristic->has->refinement_boundary)
                - target_region_start);
            refined_alignment = Optimal_find_path(
                    hpair->heuristic->optimal, region,
                    hpair->user_data, alignment->score, hpair->subopt);
            g_assert(refined_alignment);
            Region_destroy(region);
            break;
        default:
            g_error("Bad type for refinement [%s]",
                Heuristic_Refinement_to_string(
                       hpair->heuristic->has->refinement));
            break;
        }
    g_assert(refined_alignment);
    g_assert(refined_alignment->score >= alignment->score);
    return refined_alignment;
    }
#endif /* 0 */

Alignment *HPair_next_path(HPair *hpair, C4_Score threshold){
    register Alignment *alignment;
    register SAR_Alignment *sar_alignment;
    register BSDP_Path *bsdp_path;
    register BSDP_Node *first_node, *last_node,
                       *src_node, *dst_node;
    register gint i;
    register HPair_EdgeData *hpair_edge_data;
    register HPair_NodeData *first_node_data, *last_node_data,
                            *src_node_data, *dst_node_data;
    g_assert(hpair->is_finalised);
    bsdp_path = BSDP_next_path(hpair->bsdp, threshold);
    if(!bsdp_path)
        return NULL;
    g_assert(bsdp_path->node_list->len); /* Check path not empty */
    first_node = bsdp_path->node_list->pdata[0];
    last_node = bsdp_path->node_list->pdata
               [bsdp_path->node_list->len-1];
    /**/
    first_node_data = first_node->node_data;
    last_node_data = last_node->node_data;
    /**/
    sar_alignment = SAR_Alignment_create(first_node_data->sar_start,
                                         last_node_data->sar_end,
                                         first_node_data->match,
                                         last_node_data->match,
                                         hpair, bsdp_path->score);
    /* Add first HSP traceback */
    SAR_Alignment_add_HSP(sar_alignment,
                      HPair_NodeData_get_hsp(hpair, first_node_data),
                      first_node_data->match);
    for(i = 1; i < bsdp_path->node_list->len; i++){
        src_node = bsdp_path->node_list->pdata[i-1];
        dst_node = bsdp_path->node_list->pdata[i];
        g_assert(src_node->edge.used->mailbox != -1);
        hpair_edge_data = src_node->edge.used->edge_data;
        src_node_data = src_node->node_data;
        dst_node_data = dst_node->node_data;
        if(hpair_edge_data->sar_join){
            SAR_Alignment_add_SAR_Join(sar_alignment,
                                       hpair_edge_data->sar_join);
        } else {
            SAR_Alignment_add_SAR_Span(sar_alignment,
                                       hpair_edge_data->sar_span);
            }
        SAR_Alignment_add_HSP(sar_alignment,
                  HPair_NodeData_get_hsp(hpair, dst_node_data),
                  dst_node_data->match);
        }
    SAR_Alignment_finalise(sar_alignment);
    BSDP_Path_destroy(bsdp_path);
    alignment = Alignment_share(sar_alignment->alignment);
    SAR_Alignment_destroy(sar_alignment);
#if 0
    if(hpair->heuristic->has->refinement != Heuristic_Refinement_NONE){
        refined_alignment = HPair_refine_alignment(hpair, alignment);
        Alignment_destroy(alignment);
        alignment = refined_alignment;
        /* Refinement may put alignment below threshold */
        if(alignment->score < threshold){
            Alignment_destroy(alignment);
            return NULL;
            }
        }
#endif /* 0 */
    g_assert(alignment->score >= threshold);
#if 0
    SubOpt_add_alignment(hpair->subopt, alignment);
#endif /* 0 */
    return alignment;
    }

