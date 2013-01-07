/****************************************************************\
*                                                                *
*  C4 dynamic programming library - sub-alignment regions        *
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

#include "sar.h"

/**/

SAR_ArgumentSet *SAR_ArgumentSet_create(Argument *arg){
    register ArgumentSet *as;
    static SAR_ArgumentSet sas;
    if(arg){
        as = ArgumentSet_create("SAR Options");
        /**/
        ArgumentSet_add_option(as, '\0', "quality", "percent",
            "HSP quality threshold", "0",
            Argument_parse_float, &sas.hsp_quality);
        /**/
        Argument_absorb_ArgumentSet(arg, as);
        }
    return &sas;
    }

/**/

static gboolean SAR_region_valid_hsp_exit(Region *region, HSP *hsp){
#ifndef G_DISABLE_ASSERT
    register gint query_width = Region_query_end(region)
                              - hsp->query_start;
    register gint target_width = Region_target_end(region)
                               - hsp->target_start;
    register gint query_prefix = Region_query_end(region)
                               - hsp->query_start;
#endif /* G_DISABLE_ASSERT */
    /* Check region corner meets the HSP */
    g_assert(Region_is_valid(region));
    g_assert((query_width * HSP_target_advance(hsp))
          == (target_width * HSP_query_advance(hsp)));
    /* Check region corner is at HSP word boundary */
    g_assert(!(query_width % HSP_query_advance(hsp)));
    g_assert(!(target_width % HSP_target_advance(hsp)));
    /* Check region boundary is between cobs and HSP end */
    g_assert(query_prefix >= 0);
    g_assert(query_prefix <= HSP_query_cobs(hsp));
    return TRUE;
    }

static gboolean SAR_region_valid_hsp_entry(Region *region, HSP *hsp){
#ifndef G_DISABLE_ASSERT
    register gint query_width = HSP_query_end(hsp)
                              - region->query_start;
    register gint target_width = HSP_target_end(hsp)
                               - region->target_start;
    register gint query_suffix = HSP_query_end(hsp)
                               - region->query_start;
#endif /* G_DISABLE_ASSERT */
    g_assert(Region_is_valid(region));
    /* Check region corner meets the HSP */
    g_assert((query_width * HSP_target_advance(hsp))
          == (target_width * HSP_query_advance(hsp)));
    /* Check region corner is at HSP word boundary */
    g_assert(!(query_width % HSP_query_advance(hsp)));
    g_assert(!(target_width % HSP_target_advance(hsp)));
    /* Check region boundary is between HSP start and cobs */
    g_assert(query_suffix >= 0);
    g_assert(query_suffix <= (HSP_query_end(hsp)-HSP_query_cobs(hsp)));
    return TRUE;
    }

/**/

static gboolean SAR_Terminal_calculate_start_region(HSP *hsp,
            HPair *hpair, Heuristic_Range *range, C4_Scope scope,
            Region *region){
    register gint to_shrink;
    register gboolean at_query_edge = FALSE, at_target_edge = FALSE;
    Region outer_box;
    /* Find outer box (cobs -> corner) */
    outer_box.query_start = 0;
    outer_box.target_start = 0;
    outer_box.query_length = HSP_query_cobs(hsp);
    outer_box.target_length = HSP_target_cobs(hsp);
    /* Find inner_box (end point) */
    region->query_start = hsp->query_start;
    region->target_start = hsp->target_start;
    region->query_length = 0;
    region->target_length = 0;
    /* Grow inner_box by external range */
    region->query_start -= range->external_query;
    region->target_start -= range->external_target;
    region->query_length += range->external_query;
    region->target_length += range->external_target;
    /* Grow inner_box on HSPs by internal range */
    region->query_length += range->internal_query;
    region->target_length += range->internal_target;
    /* Trim inner_box to intersection (inner,outer) */
    if(region->query_start < outer_box.query_start){
        region->query_length -= (outer_box.query_start
                                 - region->query_start);
        region->query_start = outer_box.query_start;
        }
    if(region->target_start < outer_box.target_start){
        region->target_length -= (outer_box.target_start
                                  - region->target_start);
        region->target_start = outer_box.target_start;
        }
    to_shrink = Region_query_end(region)
              - Region_query_end(&outer_box);
    if(to_shrink > 0)
        region->query_length -= to_shrink;
    to_shrink = Region_target_end(region)
              - Region_target_end(&outer_box);
    if(to_shrink > 0)
        region->target_length -= to_shrink;
    /* Fail when region is empty */
    if(region->query_length <= 0)
        return FALSE;
    if(region->target_length <= 0)
        return FALSE;
    /* Check region is still on HSP boundary */
    g_assert(SAR_region_valid_hsp_exit(region, hsp));
    /* Check validity for non-local models */
    at_query_edge = (region->query_start == 0);
    at_target_edge = (region->target_start == 0);
    switch(scope){
        case C4_Scope_ANYWHERE:
            return TRUE;
            break;
        case C4_Scope_CORNER:
            if(at_query_edge && at_target_edge)
                return TRUE;
            break;
        case C4_Scope_EDGE:
            if(at_query_edge || at_target_edge)
                return TRUE;
            break;
        case C4_Scope_QUERY:
            if(at_query_edge)
                return TRUE;
            break;
        case C4_Scope_TARGET:
            if(at_target_edge)
                return TRUE;
            break;
        }
    return FALSE;
    }

static gboolean SAR_Terminal_calculate_end_region(HSP *hsp,
            HPair *hpair, Heuristic_Range *range, C4_Scope scope,
            Region *region){
    register gint to_shrink;
    register gboolean at_query_edge = FALSE, at_target_edge = FALSE;
    Region outer_box;
    /* Find outer box (cobs -> corner) */
    outer_box.query_start = HSP_query_cobs(hsp);
    outer_box.target_start = HSP_target_cobs(hsp);
    outer_box.query_length = hsp->hsp_set->query->len
                           - HSP_query_cobs(hsp);
    outer_box.target_length = hsp->hsp_set->target->len
                           - HSP_target_cobs(hsp);
    /**/
    /* Find inner_box (end point) */
    region->query_start = HSP_query_end(hsp);
    region->target_start = HSP_target_end(hsp);
    region->query_length = 0;
    region->target_length = 0;
    /* Grow inner_box by external range */
    region->query_length += range->external_query;
    region->target_length += range->external_target;
    /* Grow inner_box on HSPs by internal range */
    region->query_start -= range->internal_query;
    region->query_length += range->internal_query;
    region->target_start -= range->internal_target;
    region->target_length += range->internal_target;
    /* Trim inner_box to intersection (inner,outer) */
    if(Region_query_end(region) > Region_query_end(&outer_box)){
        region->query_length -= (Region_query_end(region)
                                - Region_query_end(&outer_box));
        }
    if(Region_target_end(region) > Region_target_end(&outer_box)){
        region->target_length -= (Region_target_end(region)
                                  - Region_target_end(&outer_box));
        }
    to_shrink = outer_box.query_start
              - region->query_start;
    if(to_shrink > 0){
        region->query_start += to_shrink;
        region->query_length -= to_shrink;
        }
    to_shrink = outer_box.target_start
              - region->target_start;
    if(to_shrink > 0){
        region->target_start += to_shrink;
        region->target_length -= to_shrink;
        }
    /* Fail when region is empty */
    if(region->query_length <= 0)
        return FALSE;
    if(region->target_length <= 0)
        return FALSE;
    /* Check region is still on HSP boundary */
    g_assert(SAR_region_valid_hsp_entry(region, hsp));
    /* Check validity for non-local models */
    at_query_edge = (Region_query_end(region)
                     == hsp->hsp_set->query->len);
    at_target_edge = (Region_target_end(region)
                     == hsp->hsp_set->target->len);
    switch(scope){
        case C4_Scope_ANYWHERE:
            return TRUE;
            break;
        case C4_Scope_CORNER:
            if(at_query_edge && at_target_edge)
                return TRUE;
            break;
        case C4_Scope_EDGE:
            if(at_query_edge || at_target_edge)
                return TRUE;
            break;
        case C4_Scope_QUERY:
            if(at_query_edge)
                return TRUE;
            break;
        case C4_Scope_TARGET:
            if(at_target_edge)
                return TRUE;
            break;
        }
    return FALSE;
    }

/**/

static C4_Score SAR_find_start_component(HPair *hpair,
                Region *region, HSP *hsp, C4_Calc *calc, gint *prefix){
    register C4_Score component = 0;
    register gint i, query_pos = hsp->query_start,
                     target_pos = hsp->target_start;
    register Region *calc_region;
    g_assert(prefix);
    (*prefix) = (Region_query_end(region) - hsp->query_start)
              / HSP_query_advance(hsp);
    g_assert((*prefix) >= 0);
    calc_region = Region_create(hsp->query_start, hsp->target_start,
                     (*prefix) * HSP_query_advance(hsp),
                     (*prefix) * HSP_target_advance(hsp));
    C4_Calc_init(calc, calc_region, hpair->user_data);
    for(i = 0; i < (*prefix); i++){
        component += C4_Calc_score(calc, query_pos, target_pos,
                                   hpair->user_data);
        query_pos += HSP_query_advance(hsp);
        target_pos += HSP_target_advance(hsp);
        }
    g_assert(component <= hsp->score);
    C4_Calc_exit(calc, calc_region, hpair->user_data);
    Region_destroy(calc_region);
    return component;
    }

static C4_Score SAR_find_end_component(HPair *hpair,
                Region *region, HSP *hsp, C4_Calc *calc, gint *suffix){
    register C4_Score component = 0;
    register gint i, query_pos, target_pos;
    register Region *calc_region;
    g_assert(suffix >= 0);
    g_assert(suffix);
    (*suffix) = (HSP_query_end(hsp) - region->query_start)
              / HSP_query_advance(hsp);
    calc_region = Region_create(region->query_start,
                                region->target_start,
                    (*suffix) * HSP_query_advance(hsp),
                    (*suffix) * HSP_target_advance(hsp));
    C4_Calc_init(calc, calc_region, hpair->user_data);
    query_pos = region->query_start;
    target_pos = region->target_start;
    for(i = 0; i < (*suffix); i++){
        component += C4_Calc_score(calc, query_pos, target_pos,
                                   hpair->user_data);
        query_pos += HSP_query_advance(hsp);
        target_pos += HSP_target_advance(hsp);
        }
    g_assert(component <= hsp->score);
    C4_Calc_exit(calc, calc_region, hpair->user_data);
    Region_destroy(calc_region);
    return component;
    }

/**/

static void SAR_HSP_quality(HPair *hpair, C4_Calc *calc,
                            HSP *hsp, C4_Score component,
                            gint start, gint length,
                            gint *half_score, gint *max_score){
    register gint i, query_pos, target_pos;
    g_assert(half_score);
    g_assert(max_score);
    query_pos = hsp->query_start + (start * HSP_query_advance(hsp));
    target_pos = hsp->target_start + (start * HSP_target_advance(hsp));
    (*half_score) = (*max_score) = 0;
    for(i = 0; i < length; i++){
        (*half_score) += HSP_get_score(hsp, query_pos, target_pos);
        (*max_score) += HSP_query_self(hsp, query_pos);
        query_pos += HSP_query_advance(hsp);
        target_pos += HSP_target_advance(hsp);
        }
    return;
    }

SAR_Terminal *SAR_Terminal_create(HSP *hsp, HPair *hpair,
                                  Heuristic_Match *match,
                                  gboolean is_start){
    register SAR_Terminal *sar_terminal;
    Region region;
    register C4_Score component;
    gint prefix, suffix;
    register SAR_ArgumentSet *sas = SAR_ArgumentSet_create(NULL);
    register gboolean is_valid = FALSE;
    gint half_score, max_score;
    register gint start, length;
    prefix = suffix = 0;
    if(is_start){
        is_valid = SAR_Terminal_calculate_start_region(hsp, hpair,
                       match->start_terminal->range,
                       hpair->heuristic->model->start_state->scope,
                       &region);
    } else { /* is end */
        is_valid = SAR_Terminal_calculate_end_region(hsp, hpair,
                       match->end_terminal->range,
                       hpair->heuristic->model->end_state->scope,
                       &region);
        }
    if(!is_valid)
        return NULL;
    if(is_start){
        component = SAR_find_start_component(hpair,
                                   &region, hsp, match->portal->calc,
                                   &prefix);
        start = prefix;
        length = hsp->cobs - prefix;
    } else {
        component = SAR_find_end_component(hpair,
                                   &region, hsp, match->portal->calc,
                                   &suffix);
        start = hsp->cobs;
        length = hsp->length - hsp->cobs - suffix;
        }
    if(length && (sas->hsp_quality > 0.0)){
        SAR_HSP_quality(hpair, match->portal->calc,
                        hsp, component, start, length,
                        &half_score, &max_score);
        if((((gdouble)half_score/(gdouble)max_score) * 100.0)
          < sas->hsp_quality)
            return NULL;
        }
    sar_terminal = g_new(SAR_Terminal, 1);
    sar_terminal->region = Region_copy(&region);
    sar_terminal->component = component;
    return sar_terminal;
    }

void SAR_Terminal_destroy(SAR_Terminal *sar_terminal){
    g_assert(sar_terminal);
    Region_destroy(sar_terminal->region);
    g_free(sar_terminal);
    return;
    }

C4_Score SAR_Terminal_find_bound(SAR_Terminal *sar_terminal,
                                 Heuristic_Bound *heuristic_bound){
    g_assert(sar_terminal->region->query_length >= 0);
    g_assert(sar_terminal->region->target_length >= 0);
    g_assert(sar_terminal->region->query_length
            <= heuristic_bound->query_range);
    g_assert(sar_terminal->region->target_length
            <= heuristic_bound->target_range);
    return heuristic_bound->matrix[sar_terminal->region->query_length]
                                  [sar_terminal->region->target_length]
           - sar_terminal->component;
    }

C4_Score SAR_Terminal_find_score(SAR_Terminal *sar_terminal,
                                 Optimal *optimal, HPair *hpair){
    return Optimal_find_score(optimal, sar_terminal->region,
                              hpair->user_data, hpair->subopt)
          - sar_terminal->component;
    }

/**/

static void SAR_reduce_mid_overlap(HPair *hpair, HSP *src, HSP *dst,
            Region *region, C4_Calc *src_calc, C4_Calc *dst_calc){
    register C4_Score src_total, dst_total, max_total;
    register gint src_query_pos, src_target_pos,
                  dst_query_pos, dst_target_pos;
    register gint max_src_query_pos, max_src_target_pos,
                  max_dst_query_pos, max_dst_target_pos,
                  max_dist_from_centre;
    /* Reverse up dst region, counting dst_total */
    if(!(region->query_length + region->target_length))
        return;
    src_total = dst_total = 0;
    dst_query_pos = Region_query_end(region)
                  - HSP_query_advance(dst);
    dst_target_pos = Region_target_end(region)
                   - HSP_target_advance(dst);
    while((dst_query_pos >= region->query_start)
       && (dst_target_pos >= region->target_start)
       && (dst_query_pos >= dst->query_start)
       && (dst_target_pos >= dst->target_start)){
        dst_total += C4_Calc_score(dst_calc,
                                   dst_query_pos, dst_target_pos,
                                   hpair->user_data);
        dst_query_pos -= HSP_query_advance(dst);
        dst_target_pos -= HSP_target_advance(dst);
        }
    dst_query_pos += HSP_query_advance(dst);
    dst_target_pos += HSP_target_advance(dst);
    /* Go down src region, storing max total */
    src_query_pos = region->query_start;
    src_target_pos = region->target_start;
    max_total = dst_total;
    max_src_query_pos = src_query_pos;
    max_src_target_pos = src_target_pos;
    max_dst_query_pos = dst_query_pos;
    max_dst_target_pos = dst_target_pos;
    max_dist_from_centre = Region_query_end(region)
                         - src_query_pos;
    while((src_query_pos < Region_query_end(region))
       && (src_target_pos < Region_target_end(region))
       && (src_query_pos < HSP_query_end(src))
       && (src_target_pos < HSP_target_end(src))){
        src_total += C4_Calc_score(src_calc,
                            src_query_pos, src_target_pos,
                            hpair->user_data);
        while((src_query_pos >= dst_query_pos)
           || (src_target_pos >= dst_target_pos)){
            dst_total -= C4_Calc_score(dst_calc,
                                       dst_query_pos, dst_target_pos,
                                       hpair->user_data);
            dst_query_pos += HSP_query_advance(dst);
            dst_target_pos += HSP_target_advance(dst);
            }
        if(max_total <= (src_total + dst_total)){
            /* If the score better
             * or the distance from the overlap centre is less
             */
            if((max_total < (src_total + dst_total))
            || (ABS(Region_query_end(region)-src_query_pos)
                < max_dist_from_centre)){
                max_dist_from_centre = ABS(Region_query_end(region)
                                           - src_query_pos);
                max_total = (src_total + dst_total);
                max_src_query_pos = src_query_pos;
                max_src_target_pos = src_target_pos;
                max_dst_query_pos = dst_query_pos;
                max_dst_target_pos = dst_target_pos;
                }
            }
        src_query_pos += HSP_query_advance(src);
        src_target_pos += HSP_target_advance(src);
        }
    /**/
    region->query_start = max_src_query_pos;
    region->target_start = max_src_target_pos;
    region->query_length = max_dst_query_pos
                         - max_src_query_pos;
    region->target_length = max_dst_target_pos
                          - max_src_target_pos;
    g_assert(Region_is_valid(region));
    return;
    }

static void SAR_find_end_box(HPair *hpair, HSP *src, HSP *dst,
                             Region *region, Region *cobs_box,
                             C4_Calc *src_calc, C4_Calc *dst_calc){
    register gint
        query_overlap = (HSP_query_end(src) - dst->query_start),
        target_overlap = (HSP_target_end(src) - dst->target_start);
    register gint src_query_move = 0, src_target_move = 0,
                  dst_query_move = 0, dst_target_move = 0;
    region->query_start = MIN(HSP_query_end(src), dst->query_start);
    region->target_start = MIN(HSP_target_end(src), dst->target_start);
    region->query_length = MAX(HSP_query_end(src), dst->query_start)
                         - region->query_start;
    region->target_length = MAX(HSP_target_end(src), dst->target_start)
                          - region->target_start;
    if((query_overlap > 0) || (target_overlap > 0)){
        /* Snap end_box to HSPs within cobs */
        src_query_move = region->query_start
                   - cobs_box->query_start;
        src_target_move = region->target_start
                    - cobs_box->target_start;
        if((src_query_move <= 0) || (src_target_move <= 0)){
            src_query_move = src_target_move = 0;
        } else {
            src_query_move -= (src_query_move % HSP_query_advance(src));
            src_target_move -= (src_target_move
                                % HSP_target_advance(src));
            if((src_query_move / HSP_query_advance(src))
             < (src_target_move / HSP_target_advance(src))){
                src_target_move = (src_query_move
                                   / HSP_query_advance(src))
                                * HSP_target_advance(src);
            } else {
                src_query_move = (src_target_move
                                  / HSP_target_advance(src))
                               * HSP_query_advance(src);
                }
            }
        /**/
        dst_query_move = Region_query_end(cobs_box)
                   - Region_query_end(region);
        dst_target_move = Region_target_end(cobs_box)
                    - Region_target_end(region);
        if((dst_query_move <= 0) || (dst_target_move <= 0)){
            dst_query_move = dst_target_move = 0;
        } else {
            dst_query_move -= (dst_query_move
                               % HSP_query_advance(dst));
            dst_target_move -= (dst_target_move
                                % HSP_target_advance(dst));
            if((dst_query_move / HSP_query_advance(dst))
             < (dst_target_move / HSP_target_advance(dst))){
                dst_target_move = (dst_query_move
                                   / HSP_query_advance(dst))
                            * HSP_target_advance(dst);
            } else {
                dst_query_move = (dst_target_move
                                  / HSP_target_advance(dst))
                               * HSP_query_advance(dst);
                }
            }
        region->query_start = cobs_box->query_start
                            + src_query_move;
        region->target_start = cobs_box->target_start
                             + src_target_move;
        region->query_length = Region_query_end(cobs_box)
                             - dst_query_move
                             - region->query_start;
        region->target_length = Region_target_end(cobs_box)
                             - dst_target_move
                             - region->target_start;
        g_assert(SAR_region_valid_hsp_entry(region, src));
        g_assert(SAR_region_valid_hsp_exit(region, dst));
        SAR_reduce_mid_overlap(hpair, src, dst, region,
                               src_calc, dst_calc);
        }
    g_assert(SAR_region_valid_hsp_entry(region, src));
    g_assert(SAR_region_valid_hsp_exit(region, dst));
    return;
    }

static gboolean SAR_find_cobs_box(HSP *src, HSP *dst, Region *region){
    region->query_start = HSP_query_cobs(src);
    region->target_start = HSP_target_cobs(src);
    region->query_length = HSP_query_cobs(dst)
                         - region->query_start;
    region->target_length = HSP_target_cobs(dst)
                          - region->target_start;
    if(region->query_length <= 0)
        return FALSE;
    if(region->target_length <= 0)
        return FALSE;
    return TRUE;
    }

static gboolean SAR_Join_calculate_region(HPair *hpair,
                                          HSP *src, HSP *dst,
                                          Region *region,
                                          Heuristic_Pair *pair){
    register gint to_shrink;
    Region outer_box;
    /* Find outer_box (cobs box) */
    if(!SAR_find_cobs_box(src, dst, &outer_box))
        return FALSE;
    /* Find inner_box (end_box) */
    SAR_find_end_box(hpair, src, dst, region, &outer_box,
                     pair->src->portal->calc, pair->dst->portal->calc);
    /* Check within external range */
    if(region->query_length > (pair->join->src_range->external_query
                             + pair->join->dst_range->external_query))
        return FALSE;
    if(region->target_length > (pair->join->src_range->external_target
                              + pair->join->dst_range->external_target))
        return FALSE;
    /* Grow inner_region by internal_range */
    region->query_start -= pair->join->src_range->internal_query;
    region->query_length += (pair->join->src_range->internal_query
                           + pair->join->dst_range->internal_query);
    region->target_start -= pair->join->src_range->internal_target;
    region->target_length += (pair->join->src_range->internal_target
                           + pair->join->dst_range->internal_target);
    /* Trim inner_box to intersection (inner,outer) */
    to_shrink = outer_box.query_start
              - region->query_start;
    if(to_shrink > 0){
        region->query_start += to_shrink;
        region->query_length -= to_shrink;
        }
    to_shrink = outer_box.target_start
              - region->target_start;
    if(to_shrink > 0){
        region->target_start += to_shrink;
        region->target_length -= to_shrink;
        }
    to_shrink = Region_query_end(region)
              - Region_query_end(&outer_box);
    if(to_shrink > 0)
        region->query_length -= to_shrink;
    to_shrink = Region_target_end(region)
              - Region_target_end(&outer_box);
    if(to_shrink > 0)
        region->target_length -= to_shrink;
    /* Check region not empty */
    if(region->query_length < 1)
        return FALSE;
    if(region->target_length < 1)
        return FALSE;
    /* Check region is still on HSP boundaries */
    g_assert(SAR_region_valid_hsp_entry(region, src));
    g_assert(SAR_region_valid_hsp_exit(region, dst));
    return TRUE;
    }

SAR_Join *SAR_Join_create(HSP *src_hsp, HSP *dst_hsp, HPair *hpair,
                          Heuristic_Pair *pair){
    register SAR_Join *sar_join;
    register C4_Score src_component, dst_component;
    register gint src_length, dst_length;
    gint prefix, suffix, src_half_score, dst_half_score,
                         src_max_score, dst_max_score;
    Region region;
    register SAR_ArgumentSet *sas = SAR_ArgumentSet_create(NULL);
    if(!SAR_Join_calculate_region(hpair, src_hsp, dst_hsp, &region,
                                  pair))
        return NULL;
    src_component = SAR_find_end_component(hpair,
              &region, src_hsp, pair->src->portal->calc, &suffix);
    dst_component = SAR_find_start_component(hpair,
              &region, dst_hsp, pair->dst->portal->calc, &prefix);
    src_length = src_hsp->length - src_hsp->cobs - suffix;
    dst_length = dst_hsp->cobs - prefix;
    if((src_length + dst_length) && (sas->hsp_quality > 0.0)){
        SAR_HSP_quality(hpair, pair->src->portal->calc,
                        src_hsp, src_component,
                        src_hsp->cobs, src_length,
                        &src_half_score, &src_max_score);
        SAR_HSP_quality(hpair, pair->dst->portal->calc,
                        dst_hsp, dst_component,
                        prefix, dst_length,
                        &dst_half_score, &dst_max_score);
        if(((((gdouble)(src_half_score+dst_half_score))
           / ((gdouble)(src_max_score+dst_max_score))) * 100.0)
          < sas->hsp_quality)
            return NULL;
        }
    sar_join = g_new(SAR_Join, 1);
    sar_join->region = Region_copy(&region);
    sar_join->src_component = src_component;
    sar_join->dst_component = dst_component;
    sar_join->pair = pair;
    return sar_join;
    }

void SAR_Join_destroy(SAR_Join *sar_join){
    g_assert(sar_join);
    Region_destroy(sar_join->region);
    g_free(sar_join);
    return;
    }

C4_Score SAR_Join_find_bound(SAR_Join *sar_join){
    g_assert(sar_join->region->query_length >= 0);
    g_assert(sar_join->region->target_length >= 0);
    g_assert(sar_join->region->query_length
            <= sar_join->pair->join->bound->query_range);
    g_assert(sar_join->region->target_length
            <= sar_join->pair->join->bound->target_range);
    return sar_join->pair->join->bound->matrix
               [sar_join->region->query_length]
               [sar_join->region->target_length]
          - (sar_join->src_component + sar_join->dst_component);
    }

C4_Score SAR_Join_find_score(SAR_Join *sar_join, HPair *hpair){
    return Optimal_find_score(sar_join->pair->join->optimal,
                              sar_join->region,
                              hpair->user_data, hpair->subopt)
          - (sar_join->src_component + sar_join->dst_component);
    }

/**/

static gboolean SAR_Span_calculate_regions(HPair *hpair,
                                           HSP *src, HSP *dst,
                                           Heuristic_Span *span,
                                           Region *src_region,
                                           Region *dst_region,
                                           C4_Calc *src_calc,
                                           C4_Calc *dst_calc){
    register gint to_shrink;
    Region outer_box, end_box;
    /* Find outer_box (cobs box) */
    if(!SAR_find_cobs_box(src, dst, &outer_box))
        return FALSE;
    /* Find inner_box (end_box) */
    SAR_find_end_box(hpair, src, dst, &end_box, &outer_box,
                     src_calc, dst_calc);
    /* Set src_region to entry points to end_box */
    src_region->query_start = end_box.query_start;
    src_region->target_start = end_box.target_start;
    src_region->query_length = 0;
    src_region->target_length = 0;
    /* Set dst_region to exit points from end_box */
    dst_region->query_start = Region_query_end(&end_box);
    dst_region->target_start = Region_target_end(&end_box);
    dst_region->query_length = 0;
    dst_region->target_length = 0;
    /* Grow src_region by external range */
    src_region->query_length += span->src_range->external_query;
    src_region->target_length += span->src_range->external_target;
    /* Grow dst_region by external range */
    dst_region->query_start -= span->dst_range->external_query;
    dst_region->target_start -= span->dst_range->external_target;
    dst_region->query_length += span->dst_range->external_query;
    dst_region->target_length += span->dst_range->external_target;
    /* Grow inner_box on HSPs by internal range */
    src_region->query_start -= span->src_range->internal_query;
    src_region->query_length += span->src_range->internal_query;
    src_region->target_start -= span->src_range->internal_target;
    src_region->target_length += span->src_range->internal_target;
    /**/
    /* Grow inner_box on HSPs by internal range */
    dst_region->query_length += span->dst_range->internal_query;
    dst_region->target_length += span->dst_range->internal_target;
    /**/
    /* Trim inner_box to intersection (inner,outer) */
    if(Region_query_end(src_region) > Region_query_end(&outer_box)){
        src_region->query_length -= (Region_query_end(src_region)
                                - Region_query_end(&outer_box));
        }
    if(Region_target_end(src_region) > Region_target_end(&outer_box)){
        src_region->target_length -= (Region_target_end(src_region)
                                  - Region_target_end(&outer_box));
        }
    to_shrink = outer_box.query_start
              - src_region->query_start;
    if(to_shrink > 0){
        src_region->query_start += to_shrink;
        src_region->query_length -= to_shrink;
        }
    to_shrink = outer_box.target_start
              - src_region->target_start;
    if(to_shrink > 0){
        src_region->target_start += to_shrink;
        src_region->target_length -= to_shrink;
        }
    /* Trim inner_box to intersection (inner,outer) */
    if(dst_region->query_start < outer_box.query_start){
        dst_region->query_length -= (outer_box.query_start
                                 - dst_region->query_start);
        dst_region->query_start = outer_box.query_start;
        }
    if(dst_region->target_start < outer_box.target_start){
        dst_region->target_length -= (outer_box.target_start
                                  - dst_region->target_start);
        dst_region->target_start = outer_box.target_start;
        }
    to_shrink = Region_query_end(dst_region)
              - Region_query_end(&outer_box);
    if(to_shrink > 0)
        dst_region->query_length -= to_shrink;
    to_shrink = Region_target_end(dst_region)
              - Region_target_end(&outer_box);
    if(to_shrink > 0)
        dst_region->target_length -= to_shrink;
    /* Check size not zero */
    if(src_region->query_length < 1)
        return FALSE;
    if(src_region->target_length < 1)
        return FALSE;
    if(dst_region->query_length < 1)
        return FALSE;
    if(dst_region->target_length < 1)
        return FALSE;
    /* Check a path is possible */
    if((dst_region->query_start - Region_query_end(src_region))
       > span->span->max_query)
        return FALSE;
    if((dst_region->target_start - Region_target_end(src_region))
       > span->span->max_target)
        return FALSE;
    /* Check region are still on HSP boundaries */
    g_assert(SAR_region_valid_hsp_entry(src_region, src));
    g_assert(SAR_region_valid_hsp_exit(dst_region, dst));
    return TRUE;
    }

SAR_Span *SAR_Span_create(HSP *src_hsp, HSP *dst_hsp, HPair *hpair,
                          Heuristic_Span *span,
                          C4_Portal *src_portal, C4_Portal *dst_portal){
    register SAR_Span *sar_span;
    register SAR_ArgumentSet *sas = SAR_ArgumentSet_create(NULL);
    register C4_Score src_component, dst_component;
    register gint src_length, dst_length;
    gint suffix, prefix, src_half_score, dst_half_score,
                         src_max_score, dst_max_score;
    Region src_region, dst_region;
    if(!SAR_Span_calculate_regions(hpair, src_hsp, dst_hsp, span,
                                   &src_region, &dst_region,
                                   src_portal->calc, dst_portal->calc))
        return NULL;
    src_component = SAR_find_end_component(hpair,
              &src_region, src_hsp, src_portal->calc, &suffix);
    dst_component = SAR_find_start_component(hpair,
              &dst_region, dst_hsp, dst_portal->calc, &prefix);
    src_length = src_hsp->length - src_hsp->cobs - suffix;
    dst_length = dst_hsp->cobs - prefix;
    if((src_length + dst_length) && (sas->hsp_quality > 0.0)){
        SAR_HSP_quality(hpair, src_portal->calc,
                        src_hsp, src_component,
                        src_hsp->cobs, src_length,
                        &src_half_score, &src_max_score);
        SAR_HSP_quality(hpair, dst_portal->calc,
                        dst_hsp, dst_component,
                        prefix, dst_length,
                        &dst_half_score, &dst_max_score);
        if(((((gdouble)(src_half_score+dst_half_score))
           / ((gdouble)(src_max_score+dst_max_score))) * 100.0)
          < sas->hsp_quality)
            return NULL;
        }
    if((100.0) < sas->hsp_quality)
        return NULL;
    sar_span = g_new(SAR_Span, 1);
    sar_span->src_region = Region_copy(&src_region);
    sar_span->dst_region = Region_copy(&dst_region);
    sar_span->src_component = src_component;
    sar_span->dst_component = dst_component;
    sar_span->span = span;
    return sar_span;
    }

void SAR_Span_destroy(SAR_Span *sar_span){
    g_assert(sar_span);
    Region_destroy(sar_span->src_region);
    Region_destroy(sar_span->dst_region);
    g_free(sar_span);
    return;
    }

C4_Score SAR_Span_find_bound(SAR_Span *sar_span){
    register C4_Score src_raw_score, dst_raw_score;
    register gint query_overlap, target_overlap;
    g_assert(sar_span->src_region->query_length
            <= sar_span->span->src_bound->query_range);
    g_assert(sar_span->src_region->target_length
            <= sar_span->span->src_bound->target_range);
    g_assert(sar_span->dst_region->query_length
            <= sar_span->span->dst_bound->query_range);
    g_assert(sar_span->dst_region->target_length
            <= sar_span->span->dst_bound->target_range);
    query_overlap = Region_query_end(sar_span->src_region)
                  - sar_span->dst_region->query_start;
    target_overlap = Region_target_end(sar_span->src_region)
                   - sar_span->dst_region->target_start;
    if(query_overlap < 0)
        query_overlap = 0;
    if(target_overlap < 0)
        target_overlap = 0;
    /**/
    /* Get the source bounds */
    src_raw_score = sar_span->span->src_bound->matrix
      [sar_span->src_region->query_length-(query_overlap>>1)]
      [sar_span->src_region->target_length-(target_overlap>>1)];
    dst_raw_score = sar_span->span->dst_bound->matrix
      [sar_span->dst_region->query_length-(query_overlap>>1)
                                         -(query_overlap & 1)]
      [sar_span->dst_region->target_length-(target_overlap>>1)
                                          -(target_overlap & 1)];
    /* Calculate the bound score */
    return (src_raw_score - sar_span->src_component)
         + (dst_raw_score - sar_span->dst_component);
    }

C4_Score SAR_Span_find_score(SAR_Span *sar_span, HPair *hpair){
    register C4_Score score;
    register Heuristic_Data *heuristic_data;
    Heuristic_Span_register(sar_span->span, sar_span->src_region,
                                            sar_span->dst_region);
    /* Calculate src optimal, reporting to the integration matrix */
    heuristic_data = hpair->user_data;
    heuristic_data->heuristic_span = sar_span->span;
    score = Optimal_find_score(sar_span->span->src_optimal,
                               sar_span->src_region, hpair->user_data,
                               hpair->subopt);
    Heuristic_Span_integrate(sar_span->span, sar_span->src_region,
                                             sar_span->dst_region);
    /* Calculate dst optimal, starting from the integration matrix */
    score = Optimal_find_score(sar_span->span->dst_optimal,
                               sar_span->dst_region, hpair->user_data,
                               hpair->subopt);
    heuristic_data->heuristic_span = NULL;
    return score - (sar_span->src_component + sar_span->dst_component);
    }

/**/

SAR_Alignment *SAR_Alignment_create(SAR_Terminal *sar_start,
                                    SAR_Terminal *sar_end,
                                    Heuristic_Match *start_match,
                                    Heuristic_Match *end_match,
                                    HPair *hpair, C4_Score score){
    register SAR_Alignment *sar_alignment = g_new(SAR_Alignment, 1);
    register Region *align_region;
    register Alignment *start_alignment
      = Optimal_find_path(start_match->start_terminal->optimal,
                          sar_start->region, hpair->user_data,
                          C4_IMPOSSIBLY_LOW_SCORE, hpair->subopt);
    g_assert(start_alignment);
    sar_alignment->end_alignment
      = Optimal_find_path(end_match->end_terminal->optimal,
                          sar_end->region, hpair->user_data,
                          C4_IMPOSSIBLY_LOW_SCORE, hpair->subopt);
    sar_alignment->end_region = Region_share(sar_end->region);
    sar_alignment->end_match = end_match;
    g_assert(sar_alignment->end_alignment);
    align_region = Region_create(start_alignment->region->query_start,
                                 start_alignment->region->target_start,
        Region_query_end(sar_alignment->end_alignment->region)
        - start_alignment->region->query_start,
        Region_target_end(sar_alignment->end_alignment->region)
        - start_alignment->region->target_start);
    sar_alignment->alignment = Alignment_create(hpair->heuristic->model,
                                                align_region, score);
    Alignment_import_derived(sar_alignment->alignment, start_alignment,
                           start_match->start_terminal->terminal_model);
    sar_alignment->hpair = hpair;
    sar_alignment->last_region = Region_share(start_alignment->region);
    sar_alignment->last_hsp = NULL;
    sar_alignment->last_match = NULL;
    Alignment_destroy(start_alignment);
    Region_destroy(align_region);
    return sar_alignment;
    }

void SAR_Alignment_destroy(SAR_Alignment *sar_alignment){
    g_assert(sar_alignment);
    g_assert(sar_alignment->alignment);
    Alignment_destroy(sar_alignment->alignment);
    if(sar_alignment->end_alignment)
        Alignment_destroy(sar_alignment->end_alignment);
    if(sar_alignment->end_region)
        Region_destroy(sar_alignment->end_region);
    if(sar_alignment->last_region)
        Region_destroy(sar_alignment->last_region);
    g_free(sar_alignment);
    return;
    }

static void SAR_Alignment_add_region(SAR_Alignment *sar_alignment,
                            Region *src_region, Region *dst_region){
    register gint suffix;
    g_assert(sar_alignment->last_hsp);
    g_assert(sar_alignment->last_match);
    g_assert(!sar_alignment->last_region);
    suffix = (HSP_query_end(sar_alignment->last_hsp)
              - src_region->query_start)
           / HSP_query_advance(sar_alignment->last_hsp);
    Alignment_add(sar_alignment->alignment,
                  sar_alignment->last_match->transition,
                  -suffix);
    sar_alignment->last_hsp = NULL;
    sar_alignment->last_match = NULL;
    sar_alignment->last_region = Region_share(dst_region);
    return;
    }

void SAR_Alignment_finalise(SAR_Alignment *sar_alignment){
    Region seq_region;
    g_assert(sar_alignment->end_alignment);
    SAR_Alignment_add_region(sar_alignment, sar_alignment->end_region,
                                            sar_alignment->end_region);
    Alignment_import_derived(sar_alignment->alignment,
              sar_alignment->end_alignment,
              sar_alignment->end_match->end_terminal->terminal_model);
    Alignment_destroy(sar_alignment->end_alignment);
    sar_alignment->end_alignment = NULL;
    Region_destroy(sar_alignment->end_region);
    sar_alignment->end_region = NULL;
    Region_init_static(&seq_region, 0, 0,
                       sar_alignment->hpair->query_length,
                       sar_alignment->hpair->target_length);
    g_assert(Alignment_is_valid(sar_alignment->alignment, &seq_region,
                                sar_alignment->hpair->user_data));
    return;
    }

void SAR_Alignment_add_SAR_Join(SAR_Alignment *sar_alignment,
                                SAR_Join *sar_join){
    register Alignment *edge_alignment;
    edge_alignment = Optimal_find_path(
            sar_join->pair->join->optimal,
            sar_join->region, sar_alignment->hpair->user_data,
            C4_IMPOSSIBLY_LOW_SCORE, sar_alignment->hpair->subopt);
    SAR_Alignment_add_region(sar_alignment, sar_join->region,
                                            sar_join->region);
    Alignment_import_derived(sar_alignment->alignment, edge_alignment,
            sar_join->pair->join->join_model);
    Alignment_destroy(edge_alignment);
    return;
    }

void SAR_Alignment_add_SAR_Span(SAR_Alignment *sar_alignment,
                                SAR_Span *sar_span){
    register Alignment *src_alignment, *dst_alignment;
    register gint query_span_start, target_span_start,
                  query_span_end, target_span_end;
    register Region *src_align_region;
    register Heuristic_Data *heuristic_data
                           = sar_alignment->hpair->user_data;
    heuristic_data->heuristic_span = sar_span->span;
    Heuristic_Span_register(sar_span->span, sar_span->src_region,
                                            sar_span->dst_region);
    Optimal_find_score(sar_span->span->src_optimal,
        sar_span->src_region, sar_alignment->hpair->user_data,
        sar_alignment->hpair->subopt);
    Heuristic_Span_integrate(sar_span->span, sar_span->src_region,
                                             sar_span->dst_region);
    dst_alignment = Optimal_find_path(sar_span->span->dst_optimal,
        sar_span->dst_region, sar_alignment->hpair->user_data,
        C4_IMPOSSIBLY_LOW_SCORE, sar_alignment->hpair->subopt);
    /* Merge edge traceback */
    query_span_end = dst_alignment->region->query_start
                   - sar_span->dst_region->query_start;
    target_span_end = dst_alignment->region->target_start
                    - sar_span->dst_region->target_start;
    query_span_start
        = sar_span->span->dst_integration_matrix
          [query_span_end][target_span_end].query_pos
        - sar_span->src_region->query_start;
    target_span_start
        = sar_span->span->dst_integration_matrix
          [query_span_end][target_span_end].target_pos
        - sar_span->src_region->target_start;
    src_align_region = Region_create(sar_span->src_region->query_start,
                                     sar_span->src_region->target_start,
                                     query_span_start,
                                     target_span_start);
    src_alignment = Optimal_find_path(
        sar_span->span->src_traceback_optimal, src_align_region,
        sar_alignment->hpair->user_data, C4_IMPOSSIBLY_LOW_SCORE,
        sar_alignment->hpair->subopt);
    Region_destroy(src_align_region);
    SAR_Alignment_add_region(sar_alignment, sar_span->src_region,
                                            sar_span->dst_region);
    Alignment_import_derived(sar_alignment->alignment, src_alignment,
              sar_span->span->src_traceback_model);
    /* Add traceback for the span */
    Heuristic_Span_add_traceback(sar_span->span,
               sar_alignment->alignment,
               dst_alignment->region->query_start
               -Region_query_end(src_alignment->region),
               dst_alignment->region->target_start
               -Region_target_end(src_alignment->region));
    heuristic_data->heuristic_span = NULL;
    Alignment_import_derived(sar_alignment->alignment, dst_alignment,
                sar_span->span->dst_model);
    Alignment_destroy(src_alignment);
    Alignment_destroy(dst_alignment);
    return;
    }

void SAR_Alignment_add_HSP(SAR_Alignment *sar_alignment, HSP *hsp,
                           Heuristic_Match *match){
    register gint prefix;
    g_assert(sar_alignment->last_region);
    g_assert(!sar_alignment->last_hsp);
    g_assert(!sar_alignment->last_match);
    prefix = (Region_query_end(sar_alignment->last_region)
              - hsp->query_start) / HSP_query_advance(hsp);
    Alignment_add(sar_alignment->alignment, match->transition,
                  hsp->length - prefix);
    Region_destroy(sar_alignment->last_region);
    sar_alignment->last_region = NULL;
    sar_alignment->last_hsp = hsp;
    sar_alignment->last_match = match;
    return;
    }

