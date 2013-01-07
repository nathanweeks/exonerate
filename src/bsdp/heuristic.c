/****************************************************************\
*                                                                *
*  C4 dynamic programming library - code for heuristics          *
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

#include "heuristic.h"
#include "matrix.h"
#include "optimal.h"
#include "region.h"

/**/

#if 0
gchar *Heuristic_Refinement_to_string(Heuristic_Refinement refinement){
    register gchar *name = NULL;
    switch(refinement){
        case Heuristic_Refinement_NONE:
            name = "none";
            break;
        case Heuristic_Refinement_FULL:
            name = "full";
            break;
        case Heuristic_Refinement_REGION:
            name = "region";
            break;
        default:
            g_error("Unknown Heuristic Refinement [%d]", refinement);
            break;
        }
    return name;
    }

Heuristic_Refinement Heuristic_Refinement_from_string(gchar *str){
    gchar *name[Heuristic_Refinement_TOTAL] =
        {"none", "full", "region"};
    Heuristic_Refinement refinement[Heuristic_Refinement_TOTAL] =
       {Heuristic_Refinement_NONE,
        Heuristic_Refinement_FULL,
        Heuristic_Refinement_REGION};
    register gint i;
    for(i = 0; i < Heuristic_Refinement_TOTAL; i++)
        if(!g_strcasecmp(name[i], str))
            return refinement[i];
    g_error("Unknown refinement type [%s]", str);
    return Heuristic_Refinement_NONE; /* Not reached */
    }

static gchar *Heuristic_Argument_parse_refinement(gchar *arg_string,
                                                  gpointer data){
    register Heuristic_Refinement *refinement
          = (Heuristic_Refinement*)data;
    gchar *refine_str;
    register gchar *ret_val = Argument_parse_string(arg_string,
                                                    &refine_str);
    if(ret_val)
        return ret_val;
    (*refinement) = Heuristic_Refinement_from_string(refine_str);
    return NULL;
    }
#endif /* 0 */

Heuristic_ArgumentSet *Heuristic_ArgumentSet_create(Argument *arg){
    register ArgumentSet *as;
    static Heuristic_ArgumentSet has;
    if(arg){
        as = ArgumentSet_create("Heuristic Options");
        /**/
        ArgumentSet_add_option(as, '\0', "terminalrangeint", NULL,
            "Internal terminal range", "12",
            Argument_parse_int, &has.terminal_range_internal);
        ArgumentSet_add_option(as, '\0', "terminalrangeext", NULL,
            "External terminal range", "12",
            Argument_parse_int, &has.terminal_range_external);
        /**/
        ArgumentSet_add_option(as, '\0', "joinrangeint", NULL,
            "Internal join range", "12",
            Argument_parse_int, &has.join_range_internal);
        ArgumentSet_add_option(as, '\0', "joinrangeext", NULL,
            "External join range", "12",
            Argument_parse_int, &has.join_range_external);
        /**/
        ArgumentSet_add_option(as, '\0', "spanrangeint", NULL,
            "Internal span range", "12",
            Argument_parse_int, &has.span_range_internal);
        ArgumentSet_add_option(as, '\0', "spanrangeext", NULL,
            "External span range", "12",
            Argument_parse_int, &has.span_range_external);
        /**/
#if 0
        ArgumentSet_add_option(as, '\0', "refine", NULL,
            "Alignment refinement strategy [none|full|region]", "none",
            Heuristic_Argument_parse_refinement, &has.refinement);
        ArgumentSet_add_option(as, '\0', "refineboundary", NULL,
            "Refinement region boundary", "32",
            Argument_parse_int, &has.refinement_boundary);
#endif /* 0 */
        /**/
        Argument_absorb_ArgumentSet(arg, as);
        }
    return &has;
    }

/**/

static Heuristic_Range *Heuristic_Range_create(gint internal,
                                               gint external,
                                               C4_Portal *portal){
    register Heuristic_Range *heuristic_range = g_new(Heuristic_Range,
                                                      1);
    g_assert(internal);
    g_assert(external);
    heuristic_range->internal_query = internal
                                    * portal->advance_query;
    heuristic_range->internal_target = internal
                                     * portal->advance_target;
    heuristic_range->external_query = external
                                    * portal->advance_query;
    heuristic_range->external_target = external
                                     * portal->advance_target;
    return heuristic_range;
    }

static void Heuristic_Range_destroy(Heuristic_Range *heuristic_range){
    g_free(heuristic_range);
    return;
    }

/**/

static void Heuristic_Bound_report_end_func(C4_Score *cell,
        gint cell_size, gint query_pos, gint target_pos,
        gpointer user_data){
    register C4_Score **matrix = user_data;
    matrix[query_pos][target_pos] = cell[0];
    return;
    }

static gchar *Heuristic_Bound_report_end_macro =
    "matrix[%QP][%TP] = %C[0]";

static Heuristic_Bound *Heuristic_Bound_create(C4_Model *model,
                        gint query_range, gint target_range){
    register Heuristic_Bound *heuristic_bound
           = g_new(Heuristic_Bound, 1);
    register gint i, j;
    register Optimal *optimal;
    register C4_Model *bound_model = C4_Model_copy(model);
    register gchar *optimal_name;
    register Region *region;
    register C4_Calc *calc;
    C4_Model_open(bound_model);
    heuristic_bound->query_range = query_range;
    heuristic_bound->target_range = target_range;
    heuristic_bound->matrix = (C4_Score**)Matrix2d_create(
                                  query_range+1, target_range+1,
                                  sizeof(C4_Score));
    /* Ensure bounds set to -INF for correct behaviour with:
     * when later calling Heuristic_Bound_max_region_convert()
     */
    for(i = 0; i <= query_range; i++)
        for(j = 0; j <= target_range; j++)
            heuristic_bound->matrix[i][j] = C4_IMPOSSIBLY_LOW_SCORE;
    /* Empty all calcs */
    for(i = 0; i < bound_model->calc_list->len; i++){
        calc = bound_model->calc_list->pdata[i];
        if(calc){
            calc->calc_func = NULL;
            calc->init_func = NULL;
            calc->exit_func = NULL;
            }
        }
    /* Remove all shadows */
    C4_Model_remove_all_shadows(bound_model);
    /* Strip of any extras, (no codegen anyway) */
    C4_Model_clear_codegen(bound_model);
    /* Strip any start_func, add an end_func */
    C4_Model_configure_extra(bound_model, NULL, NULL, NULL, NULL);
    C4_Model_configure_start_state(bound_model,
                       model->start_state->scope, NULL, NULL);
    C4_Model_configure_end_state(bound_model,
                       C4_Scope_ANYWHERE,
                       Heuristic_Bound_report_end_func,
                       Heuristic_Bound_report_end_macro);
    /* Collect the bound scores */
    C4_Model_close(bound_model);
    optimal_name = g_strdup_printf("Bound:%s", bound_model->name);
    optimal = Optimal_create(bound_model, optimal_name,
                             Optimal_Type_SCORE, FALSE);
    /* No codegen for bound models (as only used at startup) */
    g_free(optimal_name);
    region = Region_create(0, 0, query_range, target_range);
    Optimal_find_score(optimal, region, heuristic_bound->matrix, NULL);
    Region_destroy(region);
    Optimal_destroy(optimal);
    C4_Model_destroy(bound_model);
    return heuristic_bound;
    }

static void Heuristic_Bound_destroy(Heuristic_Bound *heuristic_bound){
    g_free(heuristic_bound->matrix);
    g_free(heuristic_bound);
    return;
    }

static void Heuristic_Bound_max_region_convert(
            Heuristic_Bound *heuristic_bound){
    register gint i, j;
    for(i = 1; i <= heuristic_bound->query_range; i++){
        for(j = 1; j <= heuristic_bound->target_range; j++){
            if(heuristic_bound->matrix[i][j]
             < heuristic_bound->matrix[i-1][j-1])
               heuristic_bound->matrix[i][j]
             = heuristic_bound->matrix[i-1][j-1];
            if(heuristic_bound->matrix[i][j]
             < heuristic_bound->matrix[i-1][j])
               heuristic_bound->matrix[i][j]
             = heuristic_bound->matrix[i-1][j];
            if(heuristic_bound->matrix[i][j]
             < heuristic_bound->matrix[i][j-1])
               heuristic_bound->matrix[i][j]
             = heuristic_bound->matrix[i][j-1];
            }
        }
    return;
    }
/* Converts the scores in the Heuristic_Bound,
 * so that each cell contains the max of all preceeding cells.
 */

/**/

static Heuristic_Join *Heuristic_Join_create(C4_Model *model,
                       Heuristic_Match *src, Heuristic_Match *dst,
                       Heuristic_ArgumentSet *has){
    register Heuristic_Join *heuristic_join;
    register C4_Model *derived;
    register gchar *optimal_name;
    heuristic_join = g_new(Heuristic_Join, 1);
    heuristic_join->src_range = Heuristic_Range_create(
            has->join_range_internal, has->join_range_external,
            src->portal);
    heuristic_join->dst_range = Heuristic_Range_create(
            has->join_range_internal, has->join_range_external,
            dst->portal);
    heuristic_join->join_model = C4_DerivedModel_create(model,
            src->transition->output, dst->transition->output,
            C4_Scope_CORNER, NULL, NULL,
            C4_Scope_CORNER, NULL, NULL);
    derived = heuristic_join->join_model->derived;
    heuristic_join->bound = Heuristic_Bound_create(derived,
            heuristic_join->src_range->internal_query
          + heuristic_join->src_range->external_query
          + heuristic_join->src_range->internal_query
          + heuristic_join->src_range->external_query,
            heuristic_join->dst_range->internal_target
          + heuristic_join->dst_range->external_target
          + heuristic_join->dst_range->internal_target
          + heuristic_join->dst_range->external_target);
    optimal_name = g_strdup_printf("HeuristicJoin:[%s]", derived->name);
    heuristic_join->optimal = Optimal_create(derived, optimal_name,
                                             Optimal_Type_SCORE
                                            |Optimal_Type_PATH, TRUE);
    g_free(optimal_name);
    return heuristic_join;
    }

static void Heuristic_Join_destroy(Heuristic_Join *heuristic_join){
    Heuristic_Range_destroy(heuristic_join->src_range);
    Heuristic_Range_destroy(heuristic_join->dst_range);
    C4_DerivedModel_destroy(heuristic_join->join_model);
    Heuristic_Bound_destroy(heuristic_join->bound);
    Optimal_destroy(heuristic_join->optimal);
    g_free(heuristic_join);
    return;
    }

/**/

static Heuristic_Terminal *Heuristic_Terminal_create(C4_Model *model,
      C4_Portal *portal, C4_Transition *transition, gboolean is_start,
      Heuristic_ArgumentSet *has){
    register Heuristic_Terminal *heuristic_terminal
            = g_new(Heuristic_Terminal, 1);
    register C4_Model *derived;
    register gchar *optimal_name;
    heuristic_terminal->range = Heuristic_Range_create(
            has->terminal_range_internal,
            has->terminal_range_external,
            portal);
    if(is_start){
        heuristic_terminal->terminal_model = C4_DerivedModel_create(
           model, model->start_state->state, transition->output,
           model->start_state->scope, NULL, NULL,
           C4_Scope_CORNER, NULL, NULL);
    } else { /* end terminal */
        heuristic_terminal->terminal_model = C4_DerivedModel_create(
           model, transition->output, model->end_state->state,
           C4_Scope_CORNER, NULL, NULL,
           model->end_state->scope, NULL, NULL);
        }
    derived = heuristic_terminal->terminal_model->derived;
    heuristic_terminal->bound = Heuristic_Bound_create(derived,
        (heuristic_terminal->range->internal_query
       + heuristic_terminal->range->external_query),
        (heuristic_terminal->range->internal_target
       + heuristic_terminal->range->external_target));
    Heuristic_Bound_max_region_convert(heuristic_terminal->bound);
    optimal_name = g_strdup_printf("Terminal [%s]", derived->name);
    heuristic_terminal->optimal = Optimal_create(derived, optimal_name,
                                                 Optimal_Type_SCORE
                                                |Optimal_Type_PATH,
                                                TRUE);
    g_free(optimal_name);
    return heuristic_terminal;
    }

static void Heuristic_Terminal_destroy(
            Heuristic_Terminal *heuristic_terminal){
    Heuristic_Range_destroy(heuristic_terminal->range);
    C4_DerivedModel_destroy(heuristic_terminal->terminal_model);
    Heuristic_Bound_destroy(heuristic_terminal->bound);
    Optimal_destroy(heuristic_terminal->optimal);
    g_free(heuristic_terminal);
    return;
    }

/**/

static Heuristic_Match *Heuristic_Match_create(C4_Model *model,
                C4_Portal *portal, C4_Transition *transition, gint id,
                Heuristic_ArgumentSet *has){
    register Heuristic_Match *match = g_new(Heuristic_Match, 1);
    match->id = id;
    match->portal = portal;
    match->transition = transition;
    match->start_terminal = Heuristic_Terminal_create(model, portal,
                                transition, TRUE, has);
    match->end_terminal = Heuristic_Terminal_create(model, portal,
                                transition, FALSE, has);
    return match;
    }

static void Heuristic_Match_destroy(Heuristic_Match *match){
    Heuristic_Terminal_destroy(match->start_terminal);
    Heuristic_Terminal_destroy(match->end_terminal);
    g_free(match);
    return;
    }

/**/

static C4_Score Heuristic_Span_score(Heuristic_Span *heuristic_span,
                              gint query_length, gint target_length){
    /* FIXME: need to implement properly for non-zero span scores */
    return 0;
    }

void Heuristic_Span_add_traceback(Heuristic_Span *heuristic_span,
         Alignment *alignment, gint query_length, gint target_length){
    if(query_length){
        g_assert(heuristic_span->span->query_loop);
        Alignment_add(alignment, heuristic_span->span->query_loop,
            query_length
            / heuristic_span->span->query_loop->advance_query);
        }
    if(target_length){
        g_assert(heuristic_span->span->target_loop);
        Alignment_add(alignment, heuristic_span->span->target_loop,
            target_length
            / heuristic_span->span->target_loop->advance_target);
        }
    return;
    }

static void Heuristic_Span_src_report_end_func(C4_Score *cell,
        gint cell_size, gint query_pos, gint target_pos,
        gpointer user_data){
    register Heuristic_Data *heuristic_data = user_data;
    /* This requires heuristic_data to be the 1st member of user_data */
    register Heuristic_Span *heuristic_span
                           = heuristic_data->heuristic_span;
    register gint real_query_pos, real_target_pos;
    register gint i;
    g_assert(heuristic_span);
    g_assert(heuristic_span->curr_src_region);
    real_query_pos = query_pos
                   - heuristic_span->curr_src_region->query_start;
    real_target_pos = target_pos
                    - heuristic_span->curr_src_region->target_start;
    g_assert(real_query_pos >= 0);
    g_assert(real_target_pos >= 0);
    g_assert(real_query_pos
            <= heuristic_span->src_bound->query_range);
    g_assert(real_target_pos
             <= heuristic_span->src_bound->target_range);
    for(i = 0; i < cell_size; i++)
        heuristic_span->src_integration_matrix
            [real_query_pos][real_target_pos][i] = cell[i];
    return;
    }

static C4_Score *Heuristic_Span_dst_init_start_func(
                     gint query_pos, gint target_pos,
                     gpointer user_data){
    register Heuristic_Data *heuristic_data = user_data;
    /* This requires heuristic_data to be the 1st member of user_data */
    register Heuristic_Span *heuristic_span
                           = heuristic_data->heuristic_span;
    register gint real_query_pos, real_target_pos;
    register Heuristic_Span_Cell *span_cell;
    g_assert(heuristic_span);
    g_assert(heuristic_span->curr_dst_region);
    real_query_pos = query_pos
                   - heuristic_span->curr_dst_region->query_start;
    real_target_pos = target_pos
                    - heuristic_span->curr_dst_region->target_start;
    g_assert(real_query_pos >= 0);
    g_assert(real_target_pos >= 0);
    g_assert(real_query_pos
            <= heuristic_span->dst_bound->query_range);
    g_assert(real_target_pos
            <= heuristic_span->dst_bound->target_range);
    span_cell = &heuristic_span->dst_integration_matrix
                 [real_query_pos][real_target_pos];
    if((span_cell->query_pos == -1)
    || (span_cell->target_pos == -1))
        return heuristic_span->dummy_cell;
    return heuristic_span->src_integration_matrix
           [span_cell->query_pos
           - heuristic_span->curr_src_region->query_start]
           [span_cell->target_pos
           - heuristic_span->curr_src_region->target_start];
    }

static Heuristic_Span *Heuristic_Span_create(C4_Model *model,
                       C4_State *src, C4_State *dst,
                       C4_Portal *src_portal, C4_Portal *dst_portal,
                       C4_Span *span, Heuristic_ArgumentSet *has){
    register Heuristic_Span *heuristic_span = g_new(Heuristic_Span, 1);
    register gchar *optimal_name;
    register gint cell_size = 1 + model->total_shadow_designations;
    heuristic_span->span = span;
    /**/
    heuristic_span->src_range = Heuristic_Range_create(
            has->span_range_internal, has->span_range_external,
            src_portal);
    heuristic_span->dst_range = Heuristic_Range_create(
            has->span_range_internal, has->span_range_external,
            dst_portal);
    /**/
    heuristic_span->src_model = C4_DerivedModel_create(model,
                                src, span->span_state,
        C4_Scope_CORNER, NULL, NULL, C4_Scope_ANYWHERE,
        Heuristic_Span_src_report_end_func, NULL);
    heuristic_span->dst_model = C4_DerivedModel_create(model,
                                span->span_state, dst,
        C4_Scope_ANYWHERE, Heuristic_Span_dst_init_start_func,
        NULL, C4_Scope_CORNER, NULL, NULL);
    heuristic_span->src_traceback_model = C4_DerivedModel_create(model,
                                src, span->span_state,
        C4_Scope_CORNER, NULL, NULL,
        C4_Scope_CORNER, NULL, NULL);
    /**/
    heuristic_span->src_bound = Heuristic_Bound_create(
            heuristic_span->src_model->derived,
        (heuristic_span->src_range->internal_query
       + heuristic_span->src_range->external_query),
        (heuristic_span->src_range->internal_target
       + heuristic_span->src_range->external_target));
    heuristic_span->dst_bound = Heuristic_Bound_create(
            heuristic_span->dst_model->derived,
        (heuristic_span->dst_range->internal_query
       + heuristic_span->dst_range->external_query),
        (heuristic_span->dst_range->internal_target
       + heuristic_span->dst_range->external_target));
    Heuristic_Bound_max_region_convert(heuristic_span->src_bound);
    Heuristic_Bound_max_region_convert(heuristic_span->dst_bound);
    /**/
    optimal_name = g_strdup_printf("Src span [%s]",
                             heuristic_span->src_model->derived->name);
    heuristic_span->src_optimal = Optimal_create(
                        heuristic_span->src_model->derived,
                        optimal_name,
                        Optimal_Type_SCORE|Optimal_Type_PATH, TRUE);
    g_free(optimal_name);
    /**/
    optimal_name = g_strdup_printf("Dst span [%s]",
                             heuristic_span->dst_model->derived->name);
    heuristic_span->dst_optimal = Optimal_create(
                         heuristic_span->dst_model->derived,
                         optimal_name,
                         Optimal_Type_SCORE|Optimal_Type_PATH, TRUE);
    g_free(optimal_name);
    /**/
    optimal_name = g_strdup_printf("Src span traceback [%s]",
                             heuristic_span->src_model->derived->name);
    heuristic_span->src_traceback_optimal = Optimal_create(
                         heuristic_span->src_traceback_model->derived,
                         optimal_name,
                         Optimal_Type_PATH, TRUE);
    g_free(optimal_name);
    /**/
    heuristic_span->src_integration_matrix
                     = (C4_Score***)Matrix3d_create(
                    heuristic_span->src_bound->query_range+1,
                    heuristic_span->src_bound->target_range+1,
                    cell_size, sizeof(C4_Score));
    heuristic_span->dst_integration_matrix
                     = (Heuristic_Span_Cell**)Matrix2d_create(
                    heuristic_span->dst_bound->query_range+1,
                    heuristic_span->dst_bound->target_range+1,
                    sizeof(Heuristic_Span_Cell));
    heuristic_span->curr_src_region = NULL;
    heuristic_span->curr_dst_region = NULL;
    heuristic_span->dummy_cell = g_new0(C4_Score, cell_size);
    heuristic_span->dummy_cell[0] = C4_IMPOSSIBLY_LOW_SCORE;
    return heuristic_span;
    }
/* FIXME: tidy : can use automatic optimal names instead ?
 */

static void Heuristic_Span_destroy(Heuristic_Span *heuristic_span){
    if(heuristic_span->curr_src_region)
        Region_destroy(heuristic_span->curr_src_region);
    if(heuristic_span->curr_dst_region)
        Region_destroy(heuristic_span->curr_dst_region);
    Heuristic_Range_destroy(heuristic_span->src_range);
    Heuristic_Range_destroy(heuristic_span->dst_range);
    C4_DerivedModel_destroy(heuristic_span->src_model);
    C4_DerivedModel_destroy(heuristic_span->src_traceback_model);
    C4_DerivedModel_destroy(heuristic_span->dst_model);
    Heuristic_Bound_destroy(heuristic_span->src_bound);
    Heuristic_Bound_destroy(heuristic_span->dst_bound);
    Optimal_destroy(heuristic_span->src_optimal);
    Optimal_destroy(heuristic_span->src_traceback_optimal);
    Optimal_destroy(heuristic_span->dst_optimal);
    g_free(heuristic_span->src_integration_matrix);
    g_free(heuristic_span->dst_integration_matrix);
    g_free(heuristic_span->dummy_cell);
    g_free(heuristic_span);
    return;
    }

static void Heuristic_Span_clear(Heuristic_Span *heuristic_span,
                          gint query_length, gint target_length){
    register gint i, j;
    g_assert(query_length <= heuristic_span->src_bound->query_range);
    g_assert(target_length <= heuristic_span->src_bound->target_range);
    for(i = 0; i <= heuristic_span->src_bound->query_range; i++)
        for(j = 0; j <= heuristic_span->src_bound->target_range; j++)
            heuristic_span->src_integration_matrix[i][j][0]
                                           = C4_IMPOSSIBLY_LOW_SCORE;
    return;
    }

void Heuristic_Span_register(Heuristic_Span *heuristic_span,
                             Region *src_region, Region *dst_region){
    g_assert(heuristic_span);
    g_assert(src_region);
    g_assert(dst_region);
    if(heuristic_span->curr_src_region)
        Region_destroy(heuristic_span->curr_src_region);
    if(heuristic_span->curr_dst_region)
        Region_destroy(heuristic_span->curr_dst_region);
    heuristic_span->curr_src_region = Region_share(src_region);
    heuristic_span->curr_dst_region = Region_share(dst_region);
    Heuristic_Span_clear(heuristic_span,
                         src_region->query_length,
                         src_region->target_length);
    return;
    }
/* Algorithm:
 *     for each dst matrix pos
 *         set dst matrix -INF
 *         for each valid src matrix pos
 *            challenge with score
 */

void Heuristic_Span_integrate(Heuristic_Span *heuristic_span,
                              Region *src_region, Region *dst_region){
    register gint i, j, x, y;
    register gint init_query, finish_query,
                  init_target, finish_target;
    register gint prev_init_query = -1, prev_finish_query = -1,
                  prev_init_target = -1, prev_finish_target = -1;
    register C4_Score candidate_score, top_score = 0;
    register gint top_query_pos = -1, top_target_pos = -1;
    /**/
    g_assert(heuristic_span);
    g_assert(src_region);
    g_assert(dst_region);
    for(i = 0; i <= dst_region->query_length; i++){
        for(j = 0; j <= dst_region->target_length; j++){
            heuristic_span->dst_integration_matrix[i][j].query_pos
                                           = -1;
            heuristic_span->dst_integration_matrix[i][j].target_pos
                                           = -1;
            init_query = MAX(src_region->query_start,
                            (dst_region->query_start+i
                             - heuristic_span->span->max_query));
            init_target = MAX(src_region->target_start,
                             (dst_region->target_start+j
                              - heuristic_span->span->max_target));
            finish_query = MIN((src_region->query_start
                               +src_region->query_length),
                            (dst_region->query_start+i
                             - heuristic_span->span->min_query));
            finish_target = MIN((src_region->target_start
                               +src_region->target_length),
                            (dst_region->target_start+j
                             - heuristic_span->span->min_target));
            if((init_query != prev_init_query)
            || (init_target != prev_init_target)
            || (finish_query != prev_finish_query)
            || (finish_target != prev_finish_target)){
                top_score = C4_IMPOSSIBLY_LOW_SCORE;
                top_query_pos = -1;
                top_target_pos = -1;
                for(x = init_query; x <= finish_query; x++){
                    for(y = init_target; y <= finish_target; y++){
                        candidate_score
                            = heuristic_span->src_integration_matrix
                              [x-src_region->query_start]
                              [y-src_region->target_start][0]
                            + Heuristic_Span_score(heuristic_span,
                                    dst_region->query_start+i-x,
                                    dst_region->target_start+j-y);
                        if(top_score < candidate_score){
                            top_score = candidate_score;
                            top_query_pos = x;
                            top_target_pos = y;
                            }
                        }
                    }
                }
            heuristic_span->dst_integration_matrix[i][j].query_pos
                      = top_query_pos;
            heuristic_span->dst_integration_matrix[i][j].target_pos
                      = top_target_pos;
            prev_init_query = init_query;
            prev_init_target = init_target;
            prev_finish_query = finish_query;
            prev_finish_target = finish_target;
            }
        }
#if 0
    g_message("SRC int matrix");
    for(i = 0; i <= src_region->query_length; i++){
        for(j = 0; j <= src_region->target_length; j++){
            g_print("[%d],",
                    heuristic_span->src_integration_matrix[i][j][0]);
            }
        g_print("\n");
        }
    g_message("DST int matrix");
    for(i = 0; i <= dst_region->query_length; i++){
        for(j = 0; j <= dst_region->target_length; j++){
            g_print("[%d],",
                 Heuristic_Span_dst_init_start_func(
                     dst_region->query_start+i,
                     dst_region->target_start+j,
                     NULL, heuristic_span)[0]);
            }
        g_print("\n");
        }
#endif /* 0 */
    return;
    }
/* FIXME: optimisation: could use a faster algorithm for this
 *                      with precomputed ranges or something.
 */

static gint Heuristic_Span_get_max_query_range(
            Heuristic_Span *heuristic_span){
    return heuristic_span->src_range->external_query
         + heuristic_span->dst_range->external_query
         + heuristic_span->span->max_query;
    }

static gint Heuristic_Span_get_max_target_range(
            Heuristic_Span *heuristic_span){
    return heuristic_span->src_range->external_target
         + heuristic_span->dst_range->external_target
         + heuristic_span->span->max_target;
    }

/**/

static Heuristic_Pair *Heuristic_Pair_create(C4_Model *model,
                       Heuristic_Match *src, Heuristic_Match *dst,
                       Heuristic_ArgumentSet *has){
    register Heuristic_Pair *pair;
    register gint i;
    register C4_State *src_state = src->transition->output,
                      *dst_state = dst->transition->output;
    register C4_Span *span;
    register Heuristic_Span *heuristic_span;
    if(!C4_Model_path_is_possible(model,
                                  src->transition->output,
                                  dst->transition->output))
        return NULL;
    pair = g_new(Heuristic_Pair, 1);
    pair->src = src;
    pair->dst = dst;
    pair->join = Heuristic_Join_create(model, src, dst, has);
    pair->span_list = g_ptr_array_new();
    for(i = 0; i < model->span_list->len; i++){
        span = model->span_list->pdata[i];
        if(C4_Model_path_is_possible(model, src_state,
                                     span->span_state)
        && C4_Model_path_is_possible(model, span->span_state,
                                     dst_state)){
            heuristic_span = Heuristic_Span_create(model,
                                     src_state, dst_state,
                                     src->portal, dst->portal,
                                     span, has);
            g_ptr_array_add(pair->span_list, heuristic_span);
            }
        }
    return pair;
    }

static void Heuristic_Pair_destroy(
            Heuristic_Pair *match_pair){
    register gint i;
    Heuristic_Join_destroy(match_pair->join);
    for(i = 0; i < match_pair->span_list->len; i++)
        Heuristic_Span_destroy(match_pair->span_list->pdata[i]);
    g_ptr_array_free(match_pair->span_list, TRUE);
    g_free(match_pair);
    return;
    }


void Heuristic_Pair_get_max_range(Heuristic_Pair *pair,
                                  gint *max_query_join_range,
                                  gint *max_target_join_range){
    register gint i, range;
    register Heuristic_Span *heuristic_span;
    g_assert(pair);
    g_assert(pair->join);
    (*max_query_join_range) = pair->join->src_range->external_query
                            + pair->join->dst_range->external_query;
    (*max_target_join_range) = pair->join->src_range->external_target
                             + pair->join->dst_range->external_target;
    for(i = 0; i < pair->span_list->len; i++){
        heuristic_span = pair->span_list->pdata[i];
        range = Heuristic_Span_get_max_query_range(heuristic_span);
        if((*max_query_join_range) < range);
           (*max_query_join_range) = range;
        range = Heuristic_Span_get_max_target_range(heuristic_span);
        if((*max_target_join_range) < range);
           (*max_target_join_range) = range;
        }
    return;
    }

/**/

Heuristic *Heuristic_create(C4_Model *model){
    register Heuristic *heuristic = g_new(Heuristic, 1);
    register gint i, j, counter;
    register C4_Transition *transition;
    register C4_Portal *portal;
    register Heuristic_Match *src, *dst;
    g_assert(model);
    g_assert(!model->is_open);
    g_assert(model->portal_list->len); /* Must have a portal */
    heuristic->name = g_strdup_printf("heuristic:%s", model->name);
    heuristic->ref_count = 1;
    heuristic->model = C4_Model_share(model);
    heuristic->has = Heuristic_ArgumentSet_create(NULL);
#if 0
    if(heuristic->has->refinement != Heuristic_Refinement_NONE){
        heuristic->optimal = Optimal_create(model, NULL,
                             Optimal_Type_SCORE
                            |Optimal_Type_PATH
                            |Optimal_Type_REDUCED_SPACE, TRUE);
    } else {
        heuristic->optimal = NULL;
        }
#endif /* 0 */
    /**/
    heuristic->match_total = 0;
    for(i = 0; i < model->portal_list->len; i++){
        portal = model->portal_list->pdata[i];
        heuristic->match_total += portal->transition_list->len;
        }
    g_assert(heuristic->match_total); /* Must have some matches */
    /* Prepare match_list */
    heuristic->match_list = g_new(Heuristic_Match*,
                                  heuristic->match_total);
    for(i = counter = 0; i < model->portal_list->len; i++){
        portal = model->portal_list->pdata[i];
        for(j = 0; j < portal->transition_list->len; j++){
            transition = portal->transition_list->pdata[j];
            heuristic->match_list[counter]
                = Heuristic_Match_create(model, portal, transition,
                                         counter, heuristic->has);
            counter++;
            }
        }
    g_assert(counter == heuristic->match_total);
    /* Prepare pair_matrix */
    heuristic->pair_matrix = (Heuristic_Pair***)
        Matrix2d_create(heuristic->match_total,
                        heuristic->match_total,
                        sizeof(Heuristic_Pair*));
    for(i = 0; i < heuristic->match_total; i++){
        src = heuristic->match_list[i];
        for(j = 0; j < heuristic->match_total; j++){
            dst = heuristic->match_list[j];
            heuristic->pair_matrix[i][j] =
                Heuristic_Pair_create(model, src, dst, heuristic->has);
            }
        }
    return heuristic;
    }

void Heuristic_destroy(Heuristic *heuristic){
    register gint i, j;
    g_assert(heuristic);
    if(--heuristic->ref_count)
        return;
    for(i = 0; i < heuristic->match_total; i++)
        Heuristic_Match_destroy(heuristic->match_list[i]);
    for(i = 0; i < heuristic->match_total; i++)
        for(j = 0; j < heuristic->match_total; j++)
            if(heuristic->pair_matrix[i][j])
                Heuristic_Pair_destroy(
                          heuristic->pair_matrix[i][j]);
#if 0
    if(heuristic->optimal)
        Optimal_destroy(heuristic->optimal);
#endif /* 0 */
    C4_Model_destroy(heuristic->model);
    g_free(heuristic->match_list);
    g_free(heuristic->pair_matrix);
    g_free(heuristic->name);
    g_free(heuristic);
    return;
    }

Heuristic *Heuristic_share(Heuristic *heuristic){
    g_assert(heuristic);
    heuristic->ref_count++;
    return heuristic;
    }

/**/

GPtrArray *Heuristic_get_Optimal_list(Heuristic *heuristic){
    register gint i, j, k;
    register Heuristic_Match *match;
    register Heuristic_Terminal *heuristic_terminal;
    register Heuristic_Pair *match_pair;
    register Heuristic_Join *heuristic_join;
    register Heuristic_Span *heuristic_span;
    register GPtrArray *optimal_list = g_ptr_array_new();
    g_assert(heuristic);
#if 0
    /* Refinement optimal is not included here,
     * as will be included for exhaustive alignment anyway
     */
#endif /* 0 */
    for(i = 0; i < heuristic->match_total; i++){
        match = heuristic->match_list[i];
        heuristic_terminal = match->start_terminal;
        g_ptr_array_add(optimal_list, heuristic_terminal->optimal);
        heuristic_terminal = match->end_terminal;
        g_ptr_array_add(optimal_list, heuristic_terminal->optimal);
        }
    for(i = 0; i < heuristic->match_total; i++){
        for(j = 0; j < heuristic->match_total; j++){
            if(heuristic->pair_matrix[i][j]){
                match_pair = heuristic->pair_matrix[i][j];
                heuristic_join = match_pair->join;
                g_ptr_array_add(optimal_list, heuristic_join->optimal);
                for(k = 0; k < match_pair->span_list->len; k++){
                    heuristic_span = match_pair->span_list->pdata[k];
                    g_ptr_array_add(optimal_list,
                          heuristic_span->src_optimal);
                    g_ptr_array_add(optimal_list,
                          heuristic_span->src_traceback_optimal);
                    g_ptr_array_add(optimal_list,
                          heuristic_span->dst_optimal);

                    }
                }
            }
        }
    return optimal_list;
    }

/**/

