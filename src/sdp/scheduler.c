/****************************************************************\
*                                                                *
*  C4 dynamic programming library - DP Scheduler                 *
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

#include "scheduler.h"
#include "matrix.h"
#include "cgutil.h"

#include <string.h> /* For memset()   */

static gboolean Scheduler_Pair_is_valid(Scheduler_Pair *spair);
static void Scheduler_Pair_align_rows(Scheduler_Pair *spair);

#ifdef USE_COMPILED_MODELS
#include "c4_model_archive.h"
#endif /* USE_COMPILED_MODELS */

static C4_Span **Scheduler_get_span_map(C4_Model *model){
    register gint i, j;
    register C4_Span *span,
                **span_map = g_new0(C4_Span*, model->state_list->len);
    register C4_State *state;
    for(i = 0; i < model->state_list->len; i++){
        state = model->state_list->pdata[i];
        for(j = 0; j < model->span_list->len; j++){
            span = model->span_list->pdata[j];
            if(state == span->span_state)
                span_map[state->id] = span;
            }
        }
    return span_map;
    }

Scheduler *Scheduler_create(C4_Model *model,
                      gboolean is_forward, gboolean has_traceback,
                      gboolean use_boundary,
                      Scheduler_Seed_init_Func init_func,
                      Scheduler_Seed_next_Func next_func,
                      Scheduler_Seed_get_Func get_func,
                      Scheduler_start_Func start_func,
                      Scheduler_end_Func end_func,
                      C4_Score dropoff){
    register Scheduler *scheduler = g_new(Scheduler, 1);
    register gchar *raw_name;
    register Codegen_ArgumentSet *cas
           = Codegen_ArgumentSet_create(NULL);
    g_assert(!(is_forward && start_func));
    g_assert(!(use_boundary && start_func));
    g_assert(!((!is_forward) && end_func));
    scheduler->is_forward = is_forward;
    scheduler->has_traceback = has_traceback;
    scheduler->use_boundary = use_boundary;
    /**/
    scheduler->shadow_start = 3; /* score, max, seed */
    scheduler->cell_score_offset = sizeof(Scheduler_Cell);
    scheduler->cell_traceback_offset = scheduler->cell_score_offset
                        + Matrix2d_size(model->state_list->len,
                                        scheduler->shadow_start
                                        + model->shadow_list->len,
                                        sizeof(C4_Score));
    scheduler->cell_size = scheduler->cell_traceback_offset;
    if(has_traceback)
        scheduler->cell_size += (sizeof(STraceback_Cell*)
                                 * model->state_list->len);
    scheduler->init_func = init_func;
    scheduler->next_func = next_func;
    scheduler->get_func = get_func;
    scheduler->start_func = start_func;
    scheduler->end_func = end_func;
    scheduler->model = C4_Model_share(model);
    scheduler->span_map = Scheduler_get_span_map(model);
    scheduler->dropoff = dropoff;
    /**/
    raw_name = g_strdup_printf("Scheduler_Cell_process_%s:%s:%s:%s",
                               model->name,
                               is_forward?"forward":"reverse",
                               has_traceback?"traceback":"notraceback",
                               use_boundary?"boundary":"seeded");
    scheduler->name = Codegen_clean_path_component(raw_name);
    g_free(raw_name);
    if(cas->use_compiled){
#ifdef USE_COMPILED_MODELS
        scheduler->cell_func
            = (Scheduler_DP_Func)Bootstrapper_lookup(scheduler->name);
        if(!scheduler->cell_func){
            g_warning("Could not find compiled implementation of [%s]",
                    scheduler->name);
            scheduler->cell_func = Scheduler_Cell_process;
            }
#else /* USE_COMPILED_MODELS */
        g_warning("No compiled models - will be very slow");
        scheduler->cell_func = Scheduler_Cell_process;
#endif /* USE_COMPILED_MODELS */
    } else {
        scheduler->cell_func = Scheduler_Cell_process;
        }
    g_assert(scheduler->cell_func);
    return scheduler;
    }

void Scheduler_destroy(Scheduler *scheduler){
    C4_Model_destroy(scheduler->model);
    if(scheduler->span_map)
        g_free(scheduler->span_map);
    g_free(scheduler->name);
    g_free(scheduler);
    return;
    }

static void Scheduler_end_transition(Codegen *codegen,
                                     gint transition_id){
    register gint sorry = 'o';
    Codegen_indent(codegen, 1);
    Codegen_printf(codegen, "%c%c%c%c end_of_transition_%d;\n",
                             sorry-8, sorry, sorry+5, sorry, transition_id);
    Codegen_indent(codegen, -1);
    return;
    }

static void Scheduler_implement_DP(Scheduler *scheduler,
                                   Codegen *codegen){
    register gint i, j, input_pos, output_pos;
    register C4_Transition *transition;
    register C4_Span *span;
    register C4_Shadow *shadow;
    Codegen_printf(codegen, "void %s(Scheduler_Cell *cell,\n"
                            "        Scheduler_Row *row,\n"
                            "        Scheduler_Pair *spair){\n",
                   scheduler->name);
    Codegen_indent(codegen, 1);
    Codegen_printf(codegen, "register gint seed_id = 0;\n");
    Codegen_printf(codegen, "register gint src_query_pos,\n"
                            "              src_target_pos,\n"
                            "              dst_query_pos,\n"
                            "              dst_target_pos;\n");
    Codegen_printf(codegen, "register C4_Score src_score, dst_score,\n"
                            "         transition_score, max_score;\n");
    Codegen_printf(codegen, "register C4_Transition *transition;\n");
    if(scheduler->model->span_list->len){
        Codegen_printf(codegen, "Scheduler_SpanSeed span_seed;\n");
        Codegen_printf(codegen,
                   "register Scheduler_SpanData *span_data;\n");
        }
    if(scheduler->model->shadow_list->len)
        Codegen_printf(codegen, "register C4_Shadow *shadow;\n");
    Codegen_printf(codegen, "register Scheduler_Row *dst_row;\n");
    Codegen_printf(codegen, "register Scheduler_Cell *dst_cell;\n");
    Codegen_printf(codegen, "src_query_pos = %scell->query_pos;\n",
                   scheduler->is_forward?"":"-");
    Codegen_printf(codegen, "src_target_pos = %srow->target_pos;\n",
                   scheduler->is_forward?"":"-");
    /* Calculate transitions in reverse order for forward DP */
    for(i = scheduler->model->transition_list->len-1; i >= 0; i--){
        transition = scheduler->model->transition_list->pdata[i];
        if(C4_Transition_is_span(transition)){
            if(scheduler->is_forward && scheduler->use_boundary){
                /* Copy output to span history */
                span = scheduler->span_map[transition->output->id];
                if(span){
                    input_pos = transition->input->id;
                    Codegen_printf(codegen,
                        "span_data = spair->span_history->span_data[%d];\n",
                        span->id);
                    Codegen_printf(codegen,
                        "span_seed.score = cell->score[%d][0];\n",
                        input_pos);
                    Codegen_printf(codegen, "if(span_seed.score >= 0){\n");
                    Codegen_indent(codegen, 1);
                    Codegen_printf(codegen,
                            "span_seed.max = cell->score[%d][1];\n"
                            "span_seed.seed_id = cell->score[%d][2];\n"
                            "span_seed.cell = cell->traceback[%d];\n"
                            "span_seed.query_entry = src_query_pos;\n"
                            "span_seed.target_entry = src_target_pos;\n"
                            "span_seed.shadow_data = &cell->score[%d][%d];\n",
                            input_pos, input_pos, input_pos,
                            input_pos, scheduler->shadow_start);
                    /* FIXME: func call */
                    Codegen_printf(codegen,
                        "Scheduler_SpanData_submit(span_data,\n"
                        "    &span_seed, spair->span_seed_cache, spair);\n");
                    /**/
                    Codegen_printf(codegen, "}");
                    Codegen_indent(codegen, -1);
                    }
                }
            continue;
            }
        Codegen_printf(codegen, "transition = spair->scheduler->model"
                                 "->transition_list->pdata[%d];\n", i);
        if(scheduler->is_forward){
            Codegen_printf(codegen,
                "dst_query_pos = src_query_pos + %d;\n",
                transition->advance_query);
            Codegen_printf(codegen,
                "dst_target_pos = src_target_pos + %d;\n",
                transition->advance_target);
            Codegen_printf(codegen,
          "if((dst_query_pos > Region_query_end(spair->region))\n"
          "|| (dst_target_pos > Region_target_end(spair->region)))\n");
            Scheduler_end_transition(codegen, i);
            input_pos = transition->input->id;
            output_pos = transition->output->id;
            if(scheduler->use_boundary){
                /* Copy curr best span to input */
                span = scheduler->span_map[transition->input->id];
                if(span){
                    Codegen_printf(codegen, "if(cell->permit_span_thaw){\n");
                    Codegen_indent(codegen, 1);
                    Codegen_printf(codegen,
                        "span_data =\n"
                        "    spair->span_history->span_data[%d];\n",
                        span->id);
                    /* FIXME: func call */
                    Codegen_printf(codegen,
                        "Scheduler_SpanData_get_curr(\n"
                        "  span_data, spair->span_seed_cache,\n"
                        "  cell->query_pos, row->target_pos,\n"
                        "  spair->straceback);\n");
                    /* Protect good present score from overwriting */
                    Codegen_printf(codegen,
                        "if(span_data->curr_span_seed\n"
                        "&& (cell->score[%d][0]\n"
                        "   < span_data->curr_span_seed->score)){\n",
                        input_pos);
                    Codegen_indent(codegen, 1);
                    Codegen_printf(codegen,
                            "cell->score[%d][0] = %s;\n"
                            "cell->score[%d][1] = %s;\n"
                            "cell->score[%d][2] = %s;\n",
                        input_pos,
                        "span_data->curr_span_seed->score",
                        input_pos,
                        "span_data->curr_span_seed->max",
                        input_pos,
                        "span_data->curr_span_seed->seed_id");
                        /* FIXME: func call */
                        Codegen_printf(codegen,
                        "cell->traceback[%d] = Scheduler_Cell_add_span(\n"
                        "    span_data->curr_span_seed->cell,\n"
                        "    spair->straceback,\n"
                        "    spair->scheduler->span_map[%d],\n"
                        "    src_query_pos\n"
                        "    -span_data->curr_span_seed->query_entry,\n"
                        "    src_target_pos\n"
                        "    -span_data->curr_span_seed->target_entry);\n",
                        input_pos, transition->input->id);
                    /* Thaw shadows */
                    for(j = 0;
                        j < scheduler->model->total_shadow_designations;
                        j++){
                        Codegen_printf(codegen,
                           "cell->score[%d][%d] = span_data"
                           "->curr_span_seed->shadow_data[%d];\n",
                           input_pos, scheduler->shadow_start+j, j);
                        }
                    Codegen_printf(codegen, "}\n");
                    Codegen_indent(codegen, -1);
                    Codegen_printf(codegen, "}\n");
                    Codegen_indent(codegen, -1);
                    }
                }
            /* End Shadows */
            /* This replaces Scheduler_Cell_shadow_end() */
            for(j = 0; j < transition->dst_shadow_list->len; j++){
                shadow = transition->dst_shadow_list->pdata[j];
                Codegen_printf(codegen,
                  "shadow = transition->dst_shadow_list->pdata[%d];\n",
                  j);
                Codegen_printf(codegen,
                    "shadow->end_func(cell->score[%d][%d],\n"
                    "    dst_query_pos, dst_target_pos,\n"
                    "    spair->user_data);\n",
                    input_pos,
                    scheduler->shadow_start+shadow->designation);
                }
            /* FIXME: func call */
            Codegen_printf(codegen,
                "transition_score = C4_Calc_score(transition->calc,\n"
                "   src_query_pos, src_target_pos,\n"
                "   spair->user_data);\n");
        } else { /* reverse pass */
            Codegen_printf(codegen,
                "dst_query_pos = src_query_pos - %d;\n"
                "dst_target_pos = src_target_pos - %d;\n",
                    transition->advance_query,
                    transition->advance_target);
            Codegen_printf(codegen,
                "if((dst_query_pos < spair->region->query_start)\n"
                "|| (dst_target_pos < spair->region->target_start))\n");
            Scheduler_end_transition(codegen, i);
            input_pos = transition->output->id;
            output_pos = transition->input->id;
            /* Allow extension through shadowed transitions */
            if(transition->dst_shadow_list->len)
                Codegen_printf(codegen, "transition_score = 0;\n");
            else /* FIXME: func call */
                Codegen_printf(codegen,
                    "transition_score = C4_Calc_score(transition->calc,\n"
                    "    dst_query_pos, dst_target_pos,\n"
                    "    spair->user_data);\n");
            }
        Codegen_printf(codegen,
            "src_score = cell->score[%d][0];\n"
            "max_score = cell->score[%d][1];\n"
            "seed_id   = cell->score[%d][2];\n"
            "dst_score = src_score + transition_score;\n",
            input_pos, input_pos, input_pos);
        Codegen_printf(codegen, "if(dst_score < 0)\n");
        Scheduler_end_transition(codegen, i);
        Codegen_printf(codegen, "if((max_score - dst_score)"
                                "    > spair->scheduler->dropoff)\n");
        Scheduler_end_transition(codegen, i);
        /**/
        if(C4_Transition_is_match(transition)){
            Codegen_printf(codegen, "if(spair->subopt_index){\n");
            Codegen_indent(codegen, 1);
            /* FIXME: func call */
            Codegen_printf(codegen, "if(SubOpt_Index_is_blocked(\n"
                    "spair->subopt_index,\n"
                    "src_query_pos-spair->region->query_start))\n");
            Scheduler_end_transition(codegen, i);
            Codegen_printf(codegen, "}\n");
            Codegen_indent(codegen, -1);
            }
        /* Fetch/create and set the dst cell */
        /* FIXME: func call */
        Codegen_printf(codegen,
            "dst_row = Lookahead_get(spair->row_index, %d);\n",
            transition->advance_target);
        Codegen_printf(codegen, "if(!dst_row){\n");
        Codegen_indent(codegen, 1);
        /* FIXME: func call */
        Codegen_printf(codegen,
        "dst_row = Scheduler_Row_create(%sdst_target_pos, spair);\n",
            scheduler->is_forward?"":"-");
        /* FIXME: func call */
        Codegen_printf(codegen,
                "Lookahead_set(spair->row_index, %d, dst_row);\n",
                transition->advance_target);
        Codegen_printf(codegen,
            "Lookahead_move(dst_row->cell_index, cell->query_pos);\n");
        Codegen_printf(codegen, "}\n");
        Codegen_indent(codegen, -1);
        /* FIXME: func call */
        Codegen_printf(codegen,
            "dst_cell = Lookahead_get(dst_row->cell_index, %d);\n",
            transition->advance_query);
        Codegen_printf(codegen, "if(dst_cell){\n");
        Codegen_indent(codegen, 1);
        /* Keep score if higher (don't tie break on max_score) */
        Codegen_printf(codegen,
            "if(dst_score <= dst_cell->score[%d][0])\n", output_pos);
        Scheduler_end_transition(codegen, i);
        Codegen_indent(codegen, -1);
        Codegen_printf(codegen, "} else {\n");
        Codegen_indent(codegen, 1);
        /* New cell */
        /* FIXME: func call */
        Codegen_printf(codegen,
          "dst_cell = Scheduler_Cell_create(%sdst_query_pos, FALSE, spair);\n",
           scheduler->is_forward?"":"-");
        /* FIXME: func call */
        Codegen_printf(codegen,
           "Lookahead_set(dst_row->cell_index, %d, dst_cell);\n",
           transition->advance_query);
        Codegen_printf(codegen, "}\n");
        Codegen_indent(codegen, -1);
        /* FIXME: func call */
        Codegen_printf(codegen, "Scheduler_Cell_assign(spair,\n"
                                "    cell, %d,\n"
                                "    dst_cell, %d,\n"
                                "    dst_score, max_score,\n"
                                "    transition, seed_id,\n"
                                "    dst_query_pos, dst_target_pos);\n",
                                input_pos, output_pos);
        Codegen_printf(codegen, "end_of_transition_%d:\n", i);
        }
    /**/
    Codegen_printf(codegen, "return;\n");
    Codegen_printf(codegen, "}\n\n");
    Codegen_indent(codegen, -1);
    return;
    }
/* Equivalent to Scheduler_Cell_process()
 */

Codegen *Scheduler_make_Codegen(Scheduler *scheduler){
    register Codegen *codegen = Codegen_create(NULL, scheduler->name);
    CGUtil_print_header(codegen, scheduler->model);
    Codegen_printf(codegen,
            "/* Seeded Viterbi DP implementation. */\n"
            "\n"
            "#include \"scheduler.h\"\n"
            "\n");
    Scheduler_implement_DP(scheduler, codegen);
    CGUtil_print_footer(codegen);
    CGUtil_compile(codegen, scheduler->model);
    return codegen;
    }
/* This only implements the Scheduler_Cell_process() function,
 * as is is the most time consuming part of the SDP,
 * and that which is mostly likely to benefit
 * from the codegen optimisations.
 */

/**/

static gpointer Scheduler_Cache_get(SparseCache *cache,
                                    gint slot, gint pos){
    register gpointer *cache_data = SparseCache_get(cache, pos);
    return cache_data[slot];
    }

static void Scheduler_Cache_set(SparseCache *cache,
                                gint slot, gint pos, gpointer data){
    register gpointer *cache_data = SparseCache_get(cache, pos);
    cache_data[slot] = data;
    return;
    }

/**/

static Scheduler_SpanSeed *Scheduler_SpanSeed_create(C4_Score score,
                         C4_Score max, gint seed_id,
                         gint query_entry, gint target_entry,
                         STraceback_Cell *cell,
                         C4_Score *shadow_data,
                         Scheduler *scheduler){
    register Scheduler_SpanSeed *span_seed
     = g_new(Scheduler_SpanSeed, 1);
    register gint i;
    register C4_Model *model = scheduler->model;
    span_seed->score = score;
    span_seed->max = max;
    span_seed->seed_id = seed_id;
    span_seed->query_entry = query_entry;
    span_seed->target_entry = target_entry;
    span_seed->cell = STraceback_Cell_share(cell);
    if(scheduler->is_forward
    && model->total_shadow_designations){
        span_seed->shadow_data = g_new0(C4_Score,
                                 model->total_shadow_designations);
        if(shadow_data)
            for(i = 0; i < model->total_shadow_designations; i++)
                span_seed->shadow_data[i] = shadow_data[i];
    } else {
        span_seed->shadow_data = NULL;
        }
    return span_seed;
    }

static void Scheduler_SpanSeed_destroy(Scheduler_SpanSeed *span_seed,
                                       STraceback *straceback){
    if(span_seed->shadow_data)
        g_free(span_seed->shadow_data);
    STraceback_Cell_destroy(span_seed->cell, straceback);
    g_free(span_seed);
    return;
    }

static void Scheduler_SpanSeed_copy(Scheduler_SpanSeed *src,
                                    Scheduler_SpanSeed *dst,
                                    Scheduler_Pair *spair){
    register gint i;
    register C4_Model *model = spair->scheduler->model;
    g_assert(src);
    g_assert(dst);
    dst->score        = src->score;
    dst->max          = src->max;
    dst->seed_id      = src->seed_id;
    dst->query_entry  = src->query_entry;
    dst->target_entry = src->target_entry;
    STraceback_Cell_destroy(dst->cell, spair->straceback);
    dst->cell = STraceback_Cell_share(src->cell);
    /* Copy shadows */
    if(spair->scheduler->is_forward)
        for(i = 0; i < model->total_shadow_designations; i++)
            dst->shadow_data[i] = src->shadow_data[i];
    return;
    }

static Scheduler_SpanSeed *Scheduler_SpanSeed_duplicate(
                           Scheduler_SpanSeed *span_seed,
                           Scheduler_Pair *spair){
    return Scheduler_SpanSeed_create(span_seed->score,
                               span_seed->max,
                               span_seed->seed_id,
                               span_seed->query_entry,
                               span_seed->target_entry,
                               span_seed->cell,
                               span_seed->shadow_data,
                               spair->scheduler);
    }

/**/

static gpointer Scheduler_Cache_get_func(gint pos, gpointer page_data,
                                         gpointer user_data){
    register Scheduler_SpanSeed **seed_matrix = page_data;
    return seed_matrix[pos];
    }

static SparseCache_Page *Scheduler_Cache_fill_func(gint start,
                                                   gpointer user_data){
    register Scheduler_Pair *spair = user_data;
    register Scheduler_SpanSeed ***seed_matrix
          = (Scheduler_SpanSeed***)Matrix2d_create(
                         SparseCache_PAGE_SIZE,
                         spair->scheduler->model->span_list->len,
                         sizeof(Scheduler_SpanSeed*));
    register SparseCache_Page *page = g_new(SparseCache_Page, 1);
    page->data = seed_matrix;
    page->data_size = Matrix2d_size(SparseCache_PAGE_SIZE,
                                    spair->scheduler->model->span_list->len,
                                    sizeof(Scheduler_SpanSeed*));
    page->get_func = Scheduler_Cache_get_func;
    page->copy_func = NULL;
    return page;
    }

static void Scheduler_Cache_empty_func(SparseCache_Page *page,
                                       gpointer user_data){
    register Scheduler_Pair *spair = user_data;
    register Scheduler_SpanSeed *seed, ***seed_matrix = page->data;
    register gint i, j;
    for(i = 0; i < SparseCache_PAGE_SIZE; i++)
        for(j = 0; j < spair->scheduler->model->span_list->len; j++){
            seed = seed_matrix[i][j];
            if(seed)
                Scheduler_SpanSeed_destroy(seed, spair->straceback);
            }
    g_free(seed_matrix);
    g_free(page);
    return;
    }

/**/

static Scheduler_SpanData *Scheduler_SpanData_create(
                           SListSet *slist_set, C4_Span *span){
    register Scheduler_SpanData *span_data
     = g_new(Scheduler_SpanData, 1);
    span_data->span = span;
    span_data->curr_span_seed = NULL;
    return span_data;
    }

static void Scheduler_SpanData_destroy(Scheduler_SpanData *span_data){
    g_free(span_data);
    return;
    }

void Scheduler_SpanData_get_curr(Scheduler_SpanData *span_data,
                           SparseCache *cache,
                           gint query_pos, gint target_pos,
                           STraceback *straceback){
    register Scheduler_SpanSeed *stored_seed;
    g_assert(span_data);
    g_assert(cache);
    /* Remove curr if expired */
    if(span_data->curr_span_seed){
        if((span_data->curr_span_seed->query_entry > query_pos)
        || ((span_data->curr_span_seed->query_entry
          + span_data->span->max_query) < query_pos)
        || ((span_data->curr_span_seed->target_entry
          + span_data->span->max_target) < target_pos)){
            span_data->curr_span_seed = NULL;
            }
        }
    /* Challenge curr with stored */
    stored_seed = Scheduler_Cache_get(cache,
                                      span_data->span->id, query_pos);
    if(stored_seed){
       g_assert(stored_seed->query_entry == query_pos);
       if((stored_seed->target_entry + span_data->span->max_target)
        >= target_pos){ /* Still valid */
            if(span_data->curr_span_seed){
                if(span_data->curr_span_seed->score
                  < stored_seed->score){
                    span_data->curr_span_seed = stored_seed;
                    }
            } else {
                span_data->curr_span_seed = stored_seed;
                }
        } else { /* Stored has expired */
            Scheduler_SpanSeed_destroy(stored_seed, straceback);
            Scheduler_Cache_set(cache, span_data->span->id, query_pos, NULL);
            }
        }
    return;
    }
/* Scheduler_SpanData_get_curr()
 * is called with transition from a span state
 * sets span_data->curr_span_seed to the best span_seed in scope
 */

/* Algorithm:
 * ---------
 * if(qy_span) challenging curr
 * if(ty_span) challenging cached
 * if(qy_span && ty_span)
 *    challenging curr and stored
 */

void Scheduler_SpanData_submit(Scheduler_SpanData *span_data,
                               Scheduler_SpanSeed *span_seed,
                               SparseCache *cache, Scheduler_Pair *spair){
    register Scheduler_SpanSeed *stored_seed;
    if(span_data->span->max_target){ /* Challenge cache */
        stored_seed = Scheduler_Cache_get(cache,
                                          span_data->span->id,
                                          span_seed->query_entry);
        if(stored_seed){
            g_assert(stored_seed->query_entry
                  == span_seed->query_entry);
            if(stored_seed->score <= span_seed->score){
                Scheduler_SpanSeed_copy(span_seed, stored_seed, spair);
                }
        } else {
            Scheduler_Cache_set(cache, span_data->span->id,
                            span_seed->query_entry,
                            Scheduler_SpanSeed_duplicate(span_seed, spair));
            }
        }
    return;
    }
/* Scheduler_SpanData_submit()
 * is called with transition to a span state
 */

/**/

static Scheduler_SpanHistory *Scheduler_SpanHistory_create(
                                        SListSet *slist_set,
                                        C4_Model *model){
    register gint i;
    register Scheduler_SpanHistory *span_history
     = g_new(Scheduler_SpanHistory, 1);
    register Scheduler_SpanData *span_data;
    register C4_Span *span;
    span_history->model = C4_Model_share(model);
    span_history->span_data = g_new(Scheduler_SpanData*,
                                    model->span_list->len);
    for(i = 0; i < model->span_list->len; i++){
        span = model->span_list->pdata[i];
        span_data = Scheduler_SpanData_create(slist_set, span);
        span_history->span_data[i] = span_data;
        }
    return span_history;
    }

static void Scheduler_SpanHistory_destroy(
            Scheduler_SpanHistory *span_history){
    register gint i;
    register Scheduler_SpanData *span_data;
    for(i = 0; i < span_history->model->span_list->len; i++){
        span_data = span_history->span_data[i];
        Scheduler_SpanData_destroy(span_data);
        }
    C4_Model_destroy(span_history->model);
    g_free(span_history->span_data);
    g_free(span_history);
    return;
    }

/**/



static void Scheduler_Cell_init(Scheduler_Cell *cell,
                                gint query_pos,
                                gboolean permit_span_thaw,
                                Scheduler_Pair *spair){
    register gint i;
    register C4_Model *model = spair->scheduler->model;
    memset(cell, 0, spair->scheduler->cell_size);
    cell->score = (C4_Score**)
                  ((gchar*)cell + spair->scheduler->cell_score_offset);
    Matrix2d_init((gchar**)cell->score,
            model->state_list->len,
            spair->scheduler->shadow_start + model->shadow_list->len,
            sizeof(C4_Score));
    /**/
    if(spair->scheduler->has_traceback)
        cell->traceback = (STraceback_Cell**)
                   ((gchar*)cell + spair->scheduler->cell_traceback_offset);
    else
        cell->traceback = NULL;
    if(spair->scheduler->is_forward){
        g_assert(query_pos >= 0);
    } else {
        g_assert(query_pos <= 0);
        }
    cell->query_pos = query_pos;
    cell->permit_span_thaw = permit_span_thaw;
    for(i = 0; i < spair->scheduler->model->state_list->len; i++)
        cell->score[i][0] = C4_IMPOSSIBLY_LOW_SCORE;
    if(spair->scheduler->has_traceback)
        for(i = 0; i < spair->scheduler->model->state_list->len; i++)
            cell->traceback[i] = NULL;
    return;
    }

Scheduler_Cell *Scheduler_Cell_create(gint query_pos,
                                      gboolean permit_span_thaw,
                                      Scheduler_Pair *spair){
    register Scheduler_Cell *cell;
    cell = RecycleBin_alloc(spair->cell_recycle);
    Scheduler_Cell_init(cell, query_pos, permit_span_thaw, spair);
    return cell;
    }

static void Scheduler_Cell_destroy(Scheduler_Cell *cell,
                                   Scheduler_Pair *spair){
    RecycleBin_recycle(spair->cell_recycle, cell);
    return;
    }

static void Scheduler_Cell_shadow_start(C4_Transition *transition,
                              C4_Score *dst_cell, gint shadow_start,
                              gint dst_query_pos, gint dst_target_pos,
                              gpointer user_data){
    register gint i;
    register C4_Shadow *shadow;
    for(i = 0; i < transition->input->src_shadow_list->len; i++){
        shadow = transition->input->src_shadow_list->pdata[i];
        dst_cell[shadow_start+shadow->designation]
              = shadow->start_func(
                dst_query_pos, dst_target_pos, user_data);
        }
    return;
    }
/* FIXME: do in codegen */

static void Scheduler_Cell_shadow_end(C4_Transition *transition,
                              C4_Score *src_cell, gint shadow_start,
                              gint dst_query_pos, gint dst_target_pos,
                              gpointer user_data){
    register gint i;
    register C4_Shadow *shadow;
    for(i = 0; i < transition->dst_shadow_list->len; i++){
        shadow = transition->dst_shadow_list->pdata[i];
        shadow->end_func(src_cell[shadow_start+shadow->designation],
                dst_query_pos, dst_target_pos, user_data);
        }
    return;
    }

void Scheduler_Cell_assign(Scheduler_Pair *spair,
                           Scheduler_Cell *src_cell, gint input_pos,
                           Scheduler_Cell *dst_cell, gint output_pos,
                           C4_Score dst_score, C4_Score max_score,
                           C4_Transition *transition, gint seed_id,
                           gint dst_query_pos, gint dst_target_pos){
    register C4_Model *model = spair->scheduler->model;
    register gint i;
    /*
    g_message("Assign (%d,%d)->(%d,%d) dst_score[%d] seed[%d] [%s]",
            dst_query_pos-transition->advance_query,
            dst_target_pos-transition->advance_target,
            dst_query_pos,
            dst_target_pos,
            dst_score, seed_id, transition->name);
    */
    dst_cell->score[output_pos][0] = dst_score; /* Assign score */
    dst_cell->score[output_pos][2] = seed_id;
    if(spair->scheduler->has_traceback){
        if(dst_cell->traceback[output_pos])
            STraceback_Cell_destroy(dst_cell->traceback[output_pos],
                                    spair->straceback);
        dst_cell->traceback[output_pos]
                = STraceback_add(spair->straceback, transition, 1,
                                 src_cell->traceback[input_pos]);
        /* FIXME: optimisation: overwrite prev cell directly ? */
        }
    if(spair->scheduler->is_forward){
        /* Start shadows */
        Scheduler_Cell_shadow_start(transition, src_cell->score[input_pos],
                              spair->scheduler->shadow_start,
                              dst_query_pos-transition->advance_query,
                              dst_target_pos-transition->advance_target,
                              spair->user_data);
        /* Transport shadows */
        for(i = 0; i < model->total_shadow_designations; i++){
            dst_cell->score[output_pos][spair->scheduler->shadow_start+i]
               = src_cell->score[input_pos][spair->scheduler->shadow_start+i];
            }
        /* FIXME: is this necessary ? or write shadows directly to dst ? */
        }
    if(dst_score < max_score){
        dst_cell->score[output_pos][1] = max_score;
        /*
        g_message("MAX [%d] on (%d,%d) seed=%d opp=%d (%s)",
                dst_score,
                dst_query_pos-transition->advance_query,
                dst_target_pos-transition->advance_target,
                seed_id, output_pos, transition->name);
        */
    } else { /* Is the best score on this path */
        dst_cell->score[output_pos][1] = dst_score;
        if(spair->scheduler->start_func){
            if(transition->input == model->start_state->state)
                spair->scheduler->start_func(spair->seed_data,
                                      seed_id, dst_score,
                                      dst_query_pos, dst_target_pos,
                                      spair->scheduler->has_traceback
                                      ?dst_cell->traceback[output_pos]
                                      :NULL);
            }
        if(spair->scheduler->end_func){
            if(transition->output == model->end_state->state){
                g_assert(dst_cell->traceback[output_pos]);
                spair->scheduler->end_func(spair->seed_data,
                                    seed_id, dst_score,
                                    dst_query_pos, dst_target_pos,
                                    dst_cell->traceback[output_pos]);
                }
            }
        }
    return;
    }

/**/

STraceback_Cell *Scheduler_Cell_add_span(STraceback_Cell *prev,
                                         STraceback *straceback,
                                         C4_Span *span,
                                         gint query_len, gint target_len){
    g_assert(span);
    g_assert(span->query_loop || (!query_len));
    g_assert(span->target_loop || (!target_len));
    if(query_len){
        g_assert(query_len >= 0);
        g_assert(span->query_loop);
        prev = STraceback_add(straceback, span->query_loop, query_len, prev);
        }
    if(target_len){
        g_assert(target_len >= 0);
        g_assert(span->target_loop);
        prev = STraceback_add(straceback, span->target_loop, target_len, prev);
        }
    return prev;
    }

void Scheduler_Cell_process(Scheduler_Cell *cell,
                            Scheduler_Row *row,
                            Scheduler_Pair *spair){
    register C4_Transition *transition;
    register C4_Model *model = spair->scheduler->model;
    register C4_Score src_score, dst_score, transition_score, max_score;
    register Scheduler_Row *dst_row;
    register Scheduler_Cell *dst_cell;
    register Scheduler_SpanData *span_data;
    register C4_Span *span;
    register gint i, j, seed_id = 0,
                  src_query_pos, src_target_pos,
                  dst_query_pos, dst_target_pos,
                  rel_dst_query_pos, rel_dst_target_pos,
                  input_pos, output_pos;
    Scheduler_SpanSeed span_seed;
    /*
    g_print("%s %s %d %d\n", __FUNCTION__,
            spair->scheduler->is_forward?"fwd":"rev",
            cell->query_pos, row->target_pos);
    */
    if(spair->scheduler->is_forward){
        src_query_pos = cell->query_pos;
        src_target_pos = row->target_pos;
    } else {
        src_query_pos = -cell->query_pos;
        src_target_pos = -row->target_pos;
        }
    /* Calculate transitions in reverse order for forward DP */
    for(i = model->transition_list->len-1; i >= 0; i--){
        transition = model->transition_list->pdata[i];
        if(C4_Transition_is_span(transition)){
            if(spair->scheduler->is_forward
            && spair->scheduler->use_boundary){
            /* Copy output to span history */
                span = spair->scheduler->span_map[transition->output->id];
                if(span){
                    input_pos = transition->input->id;
                    span_data = spair->span_history->span_data[span->id];
                    span_seed.score        = cell->score[input_pos][0];
                    if(span_seed.score >= 0){
                        span_seed.max          = cell->score[input_pos][1];
                        span_seed.seed_id      = cell->score[input_pos][2];
                        span_seed.cell         = cell->traceback[input_pos];
                        /* Ref count does not need to be taken here */
                        span_seed.query_entry  = src_query_pos;
                        span_seed.target_entry = src_target_pos;
                        /*
                        g_message("FREEZING [%d,%d] [%d] entry (%d,%d)",
                           span_seed.score,
                           span_seed.max,
                           span_seed.seed_id,
                           span_seed.query_entry,
                           span_seed.target_entry);
                        */
                        span_seed.shadow_data = &cell->score[input_pos]
                                       [spair->scheduler->shadow_start];
                        Scheduler_SpanData_submit(span_data, &span_seed,
                                              spair->span_seed_cache, spair);
                        }
                    }
                }
            continue;
            }
        if(spair->scheduler->is_forward){
            dst_query_pos = src_query_pos
                          + transition->advance_query;
            dst_target_pos = src_target_pos
                           + transition->advance_target;
            if((dst_query_pos > Region_query_end(spair->region))
            || (dst_target_pos > Region_target_end(spair->region)))
                continue;
            input_pos = transition->input->id;
            output_pos = transition->output->id;
            rel_dst_query_pos = dst_query_pos;
            rel_dst_target_pos = dst_target_pos;
            if(cell->permit_span_thaw){
                /* Copy curr best span to input */
                span = spair->scheduler
                     ->span_map[transition->input->id];
                if(span){
                    span_data = spair
                              ->span_history->span_data[span->id];
                    Scheduler_SpanData_get_curr(span_data,
                                                spair->span_seed_cache,
                                                cell->query_pos,
                                                row->target_pos,
                                                spair->straceback);
                    /* Protect good present score from overwriting */
                    if((span_data->curr_span_seed)
                    && (cell->score[input_pos][0]
                        < span_data->curr_span_seed->score)){
                        /*
                        g_message("THAWING [%d,%d] [%d] entry(%d,%d) [%d,%d]",
                           span_data->curr_span_seed->score,
                           span_data->curr_span_seed->max,
                           span_data->curr_span_seed->seed_id,
                           span_data->curr_span_seed->query_entry,
                           span_data->curr_span_seed->target_entry,
                           src_query_pos, src_target_pos);
                           */
                        /**/
                        cell->score[input_pos][0]
                           = span_data->curr_span_seed->score;
                        cell->score[input_pos][1]
                           = span_data->curr_span_seed->max;
                        cell->score[input_pos][2]
                           = span_data->curr_span_seed->seed_id;
                        g_assert(spair->scheduler->has_traceback);
                        cell->traceback[input_pos]
                            = Scheduler_Cell_add_span(
                                span_data->curr_span_seed->cell,
                                spair->straceback, span,
                                src_query_pos
                                -span_data->curr_span_seed->query_entry,
                                src_target_pos
                                -span_data->curr_span_seed->target_entry);
                        /**/
                        /* Thaw shadows */
                        for(j = 0;
                            j < model->total_shadow_designations;
                            j++)
                            cell->score[input_pos]
                                  [spair->scheduler->shadow_start+j]
                            = span_data->curr_span_seed->shadow_data[j];
                        }
                    }
                }
            /* End shadows */
            Scheduler_Cell_shadow_end(transition,
                                cell->score[input_pos],
                                spair->scheduler->shadow_start,
                                dst_query_pos, dst_target_pos,
                                spair->user_data);
            transition_score = C4_Calc_score(transition->calc,
                        src_query_pos, src_target_pos,
                        spair->user_data);
        } else { /* reverse pass */
            dst_query_pos = src_query_pos
                          - transition->advance_query;
            dst_target_pos = src_target_pos
                           - transition->advance_target;
            if((dst_query_pos < spair->region->query_start)
            || (dst_target_pos < spair->region->target_start))
                continue;
            rel_dst_query_pos = -dst_query_pos;
            rel_dst_target_pos = -dst_target_pos;
            input_pos = transition->output->id;
            output_pos = transition->input->id;
            /* Allow extension through shadowed transitions */
            if(transition->dst_shadow_list->len)
                transition_score = 0;
            else
                transition_score = C4_Calc_score(transition->calc,
                        dst_query_pos, dst_target_pos,
                        spair->user_data);
            }
        src_score = cell->score[input_pos][0];
        max_score = cell->score[input_pos][1];
        seed_id   = cell->score[input_pos][2];
        dst_score = src_score + transition_score;
        if(spair->scheduler->is_forward)
        if(dst_score < 0)
            continue;
        if((max_score - dst_score) > spair->scheduler->dropoff)
            continue;
/**/
        if(C4_Transition_is_match(transition)){
            if(spair->subopt_index){
                if(SubOpt_Index_is_blocked(spair->subopt_index,
                   src_query_pos-spair->region->query_start))
                    continue;
                }
            }
/**/
        /* Fetch/create and set the dst cell */
        dst_row = Lookahead_get(spair->row_index,
                                transition->advance_target);
        if(!dst_row){
            dst_row = Scheduler_Row_create(rel_dst_target_pos, spair);
            Lookahead_set(spair->row_index,
                          transition->advance_target, dst_row);
            Lookahead_move(dst_row->cell_index, cell->query_pos);
            }
        g_assert(dst_row);
        g_assert(dst_row->target_pos == rel_dst_target_pos);
        g_assert(dst_row->cell_index->pos == cell->query_pos);
        dst_cell = Lookahead_get(dst_row->cell_index,
                                 transition->advance_query);
        if(dst_cell){
            /* Keep score if higher (don't tie break on max_score) */
            if(dst_score <= dst_cell->score[output_pos][0])
                continue;
            g_assert(dst_cell->query_pos == rel_dst_query_pos);
        } else { /* New cell */
            dst_cell = Scheduler_Cell_create(rel_dst_query_pos, FALSE, spair);
            Lookahead_set(dst_row->cell_index,
                          transition->advance_query, dst_cell);
            }
        g_assert(dst_cell);
        Scheduler_Cell_assign(spair, cell, input_pos,
                              dst_cell, output_pos,
                              dst_score, max_score, transition, seed_id,
                              dst_query_pos, dst_target_pos);
        }
    return;
    }
/* FIXME: tidy and optimise */

static void Scheduler_Cell_seed(Scheduler_Cell *cell,
                                Scheduler_Seed *seed,
                                Scheduler_Pair *spair){
    register gboolean seed_state_id = spair->scheduler->is_forward
         ?  spair->scheduler->model->start_state->state->id
         :  spair->scheduler->model->end_state->state->id;
    cell->score[seed_state_id][0] = seed->start_score;
    cell->score[seed_state_id][1] = seed->start_score;
    cell->score[seed_state_id][2] = seed->seed_id;
    if(cell->traceback)
        cell->traceback[seed_state_id] = NULL;
    return;
    }

/**/

static void Scheduler_Row_Lookahead_free_func(gpointer data,
                                        gpointer user_data){
    register Scheduler_Cell *cell = data;
    register Scheduler_Row *row = user_data;
    g_assert(cell);
    /* Do not destroy cells - just add to the used_list */
    SList_queue(row->used_list, cell);
    return;
    }

Scheduler_Row *Scheduler_Row_create(gint target_pos,
                                    Scheduler_Pair *spair){
    register Scheduler_Row *row = g_new0(Scheduler_Row, 1);
    if(spair->scheduler->is_forward){
        g_assert(target_pos >= 0);
    } else {
        g_assert(target_pos <= 0);
        }
    row->target_pos = target_pos;
    row->cell_index = Lookahead_create(
          spair->scheduler->is_forward
          ?0:-Region_query_end(spair->region),
          spair->scheduler->model->max_query_advance,
          Scheduler_Row_Lookahead_free_func, row);
    row->used_list = SList_create(spair->slist_set);
    row->unused_list = SList_create(spair->slist_set);
    return row;
    }

static void Scheduler_Row_process(Scheduler_Row *row,
                                  Scheduler_Pair *spair){
    register Scheduler_Cell *cell, *next_cell;
    register gint advance;
    g_assert(Scheduler_Pair_is_valid(spair));
    if(spair->scheduler->is_forward)
        SubOpt_Index_set_row(spair->subopt_index,
                     row->target_pos-spair->region->target_start);
    else
        SubOpt_Index_set_row(spair->subopt_index,
                   (-row->target_pos)-spair->region->target_start);
    do {
        /* Add first cell if row index is empty */
        if(!row->cell_index->mask){
            if(SList_isempty(row->unused_list))
                break; /* No more cells */
            cell = SList_pop(row->unused_list);
            g_assert(cell);
            Lookahead_move(row->cell_index, cell->query_pos);
            Lookahead_set(row->cell_index, 0, cell);
            }
        /* Populate initial cell positions with seeds
         * (no more than max_query_advance columns)
         */
        cell = Lookahead_get(row->cell_index, 0);
        g_assert(cell);
        do {
            g_assert(row);
            g_assert(row->unused_list);
            next_cell = SList_first(row->unused_list);
            if(!next_cell)
                break;
            advance = next_cell->query_pos - cell->query_pos;
            g_assert(advance >= 0);
            if(advance > row->cell_index->max_advance)
                break;
            /* Cell cannot be already created */
            g_assert(!Lookahead_get(row->cell_index, advance));
            next_cell = SList_pop(row->unused_list);
            Lookahead_set(row->cell_index, advance, next_cell);
        } while(TRUE);
        /* Process cell */
        Scheduler_Pair_align_rows(spair);
        /* Will use codegen when available ... */
        spair->scheduler->cell_func(cell, row, spair);
        Lookahead_next(row->cell_index);
    } while(TRUE);
    return;
    }

typedef struct {
    Scheduler_Pair *spair;
    gint target_pos;
} Scheduler_TraverseData;

static void Scheduler_Row_traverse_cell_destroy(gpointer data,
                                                gpointer user_data){
    register Scheduler_Cell *cell = data;
    register Scheduler_TraverseData *std = user_data;
    register Scheduler_Pair *spair = std->spair;
    register gint i, state_id;
    register Boundary_Row *boundary_row;
    register C4_Model *model = spair->scheduler->model;
    register STraceback_Cell *stcell, *prev_stcell;
    register C4_Span *span;
    if(spair->straceback){
        for(i = 0; i < model->state_list->len; i++){
            stcell = cell->traceback[i];
            if(stcell){
                if((stcell->prev)
                && (stcell->prev->ref_count == 1)
                && (stcell->prev->transition == stcell->transition)){
                    prev_stcell = stcell->prev;
                    stcell->prev = prev_stcell->prev
                                 ? STraceback_Cell_share(prev_stcell->prev)
                                 : NULL;
                    stcell->length += prev_stcell->length;
                    STraceback_Cell_destroy(prev_stcell,
                                            spair->straceback);
                    }
                STraceback_Cell_destroy(cell->traceback[i],
                                        spair->straceback);
                }
            }
        }
    if((!spair->scheduler->is_forward) && spair->boundary){
        /* Record new boundary */
        boundary_row = Boundary_get_last_row(spair->boundary);
        g_assert(boundary_row);
        state_id = model->start_state->state->id;
        if(cell->score[state_id][0] >= 0){
            Boundary_Row_prepend(boundary_row, -cell->query_pos,
                                 cell->score[state_id][2]);
        } else {
            for(i = 0; i < model->span_list->len; i++){
                span = model->span_list->pdata[i];
                state_id = span->span_state->id;
                if(cell->score[state_id][0] > 0){
                    Boundary_Row_prepend(boundary_row, -cell->query_pos,
                                         cell->score[state_id][2]);
                    break;
                    }
                }
            }
        }
    /**/
    Scheduler_Cell_destroy(cell, spair);
    return;
    }

static void Scheduler_Row_destroy(Scheduler_Row *row,
                                  Scheduler_Pair *spair){
    register Boundary_Row *boundary_row = NULL;
    Scheduler_TraverseData std;
    std.target_pos = row->target_pos;
    std.spair = spair;
    if((!spair->scheduler->is_forward) && spair->boundary)
        boundary_row = Boundary_add_row(spair->boundary,
                                        -row->target_pos);
    Lookahead_destroy(row->cell_index);
    SList_destroy_with_data(row->used_list,
                            Scheduler_Row_traverse_cell_destroy, &std);
    SList_destroy_with_data(row->unused_list,
                            Scheduler_Row_traverse_cell_destroy, &std);
    if(boundary_row)
        Boundary_remove_empty_last_row(spair->boundary);
    g_free(row);
    return;
    }

/* lookahead->pos uses -ve coordinates
 * terminal coordinates are always +ve
 *     cell->query_pos
 * and row->target_pos use reversed coordinates.
 */

static void Scheduler_Row_add_seed(Scheduler_Row *row,
                                   Scheduler_Pair *spair,
                                   Scheduler_Seed *seed){
    register Scheduler_Cell *cell
           = Scheduler_Cell_create(seed->query_pos,
                                  (spair->scheduler->is_forward
                                && spair->scheduler->use_boundary), spair);
    register Scheduler_Cell *first_cell;
    register gint advance;
    if(row->cell_index->mask){ /* Already have cells */
        /* Add to lookahead when within advance range */
        first_cell = Lookahead_get(row->cell_index, 0);
        g_assert(first_cell);
        advance = cell->query_pos - first_cell->query_pos;
        g_assert(advance >= 0);
        if(advance < row->cell_index->max_advance){
            /* Row should not exist already */
            g_assert(!Lookahead_get(row->cell_index, advance));
            Lookahead_set(row->cell_index, advance, cell);
        } else { /* Otherwise add to the unused list */
            SList_queue(row->unused_list, cell);
            }
    } else { /* Row is empty */
        Lookahead_move(row->cell_index, cell->query_pos);
        Lookahead_set(row->cell_index, 0, cell);
        }
    /* FIXME: Optimisation: (FOR EFFICIENT SUBOPT DP)
     *        Should not seed here unless the seed is part of the
     *        update boundary (or the seed is an update seed)
     */
    Scheduler_Cell_seed(cell, seed, spair);
    return;
    }

static void Scheduler_Row_reset(Scheduler_Row *row){
    register SList *joined_list;
    Lookahead_reset(row->cell_index);
    /* Prepend anything in used onto unused */
    joined_list = SList_join(row->used_list, row->unused_list);
    row->used_list = row->unused_list;
    row->unused_list = joined_list;
    return;
    }

static void Scheduler_Row_move(Scheduler_Row *row, gint query_pos){
    register Scheduler_Cell *cell;
    register gint advance;
    Lookahead_move(row->cell_index, query_pos);
    /* Move anything in unused list before index onto used list */
    do {
        cell = SList_first(row->unused_list);
        if(!cell)
            break;
        if(cell->query_pos >= query_pos)
            break;
        SList_queue(row->used_list, SList_pop(row->unused_list));
    } while(TRUE);
    /* Move anything in unused list in index into index */
    do {
        cell = SList_first(row->unused_list);
        if(!cell)
            break;
        advance = cell->query_pos - query_pos;
        g_assert(advance >= 0);
        if(advance > row->cell_index->max_advance)
            break;
        Lookahead_set(row->cell_index, advance,
                      SList_pop(row->unused_list));
    } while(TRUE);
    return;
    }

/* FIXME: tidy */

static gboolean Scheduler_Row_print_traverse(gpointer data,
                                       gpointer user_data){
    register Scheduler_Cell *cell = data;
    register gchar *message = user_data;
    g_message("  >[%s] CELL [%d](%p)",
              message, cell->query_pos, cell);
    return FALSE;
    }

static void Scheduler_Row_print(Scheduler_Row *row){
    register gint i;
    register Scheduler_Cell *cell;
    g_print("\n");
    g_message("--[");
    g_message("ROW target_pos [%d] (mask %d) (ma=%d)", row->target_pos,
            row->cell_index->mask, row->cell_index->max_advance);
    g_message("index pos [%d]", row->cell_index->pos);
    for(i = 0; i <= row->cell_index->max_advance; i++){
        cell = Lookahead_get(row->cell_index, i);
        if(cell){
            g_message("  >lookahead (%d) CELL [%d]",
                    i, cell->query_pos);
            }
        }
    SList_traverse(row->unused_list,
                   Scheduler_Row_print_traverse, "unused");
    SList_traverse(row->used_list,
                   Scheduler_Row_print_traverse, "used");
    g_message("]--\n");
    return;
    }

/**/

typedef struct {
        gint pos;
    gboolean is_valid;
} Scheduler_Row_Validate;

static gboolean Scheduler_Row_is_valid_traverse(gpointer data,
                                          gpointer user_data){
    register Scheduler_Row_Validate *v = user_data;
    register Scheduler_Cell *cell = data;
    if(cell->query_pos < v->pos){
        v->is_valid = FALSE;
        }
    v->pos = cell->query_pos;
    return FALSE;
    }

static gboolean Scheduler_Row_is_valid(Scheduler_Row *row){
    Scheduler_Row_Validate v;
    register gint i;
    register Scheduler_Cell *cell;
    v.pos = -987654321;
    v.is_valid = TRUE;
    SList_traverse(row->used_list, Scheduler_Row_is_valid_traverse, &v);
    for(i = 0; i <= row->cell_index->max_advance; i++){
        cell = Lookahead_get(row->cell_index, i);
        if(cell){
            if(cell->query_pos < v.pos)
                v.is_valid = FALSE;
            }
        }
    SList_traverse(row->unused_list,
                   Scheduler_Row_is_valid_traverse, &v);
    if(!v.is_valid){
        g_message("Found invalid row");
        Scheduler_Row_print(row);
        g_error("bad row");
        }
    return v.is_valid;
    }

/**/

static void Scheduler_Pair_init(Scheduler_Pair *spair,
                                Scheduler_Seed *seed){
    /* Load first seed at start of index */
    register Scheduler_Row *first_row;
    /* Make first row */
    first_row = Scheduler_Row_create(seed->target_pos, spair);
    Scheduler_Row_add_seed(first_row, spair, seed);
    /* Set row index position */
    Lookahead_move(spair->row_index, first_row->target_pos);
    Lookahead_set(spair->row_index, 0, first_row);
    return;
    }

static void Scheduler_Pair_add_seed(Scheduler_Pair *spair,
                                    Scheduler_Seed *seed){
    register Scheduler_Row *row, *first_row;
    register gint advance;
    first_row = Lookahead_get(spair->row_index, 0);
    g_assert(first_row);
    advance = seed->target_pos - first_row->target_pos;
    g_assert(advance >= 0);
    g_assert(advance <= spair->row_index->max_advance);
    row = Lookahead_get(spair->row_index, advance);
    if(!row){
        row = Scheduler_Row_create(seed->target_pos, spair);
        Lookahead_set(spair->row_index, advance, row);
        }
    g_assert(row);
    g_assert(row->target_pos == seed->target_pos);
    Scheduler_Row_add_seed(row, spair, seed);
    return;
    }

static void Scheduler_Pair_reset_rows(Scheduler_Pair *spair){
    register gint i;
    register Scheduler_Row *row;
    for(i = 0; i <= spair->row_index->max_advance; i++){
        row = Lookahead_get(spair->row_index, i);
        if(row){
            g_assert(Scheduler_Row_is_valid(row));
            Scheduler_Row_reset(row);
            }
        }
    return;
    }

void Scheduler_Pair_calculate(Scheduler_Pair *spair){
    Scheduler_Seed seed;
    register Scheduler_Row *row;
    register gint i;
    register C4_Calc *calc;
    register C4_Model *model = spair->scheduler->model;
    /* Call model init func */
    if(model->init_func)
        model->init_func(spair->region, spair->user_data);
    /* Call model calc init funcs */
    for(i = 0; i < model->calc_list->len; i++){
        calc = model->calc_list->pdata[i];
        if(calc && calc->init_func)
            calc->init_func(spair->region, spair->user_data);
        }
    /* Call spair init func */
    spair->scheduler->init_func(spair->seed_data);
    /* Compute seeded DP rows while rows and seeds remain */
    do {
        /* Add first point if row_index is empty */
        if(!spair->row_index->mask){
            if(!spair->scheduler->get_func(spair->seed_data, &seed))
                break;
            Scheduler_Pair_init(spair, &seed);
            spair->scheduler->next_func(spair->seed_data);
            }
        /* Populate initial DP rows with seeds
         * (no more than max_target_advance rows)
         */
        row = Lookahead_get(spair->row_index, 0);
        g_assert(row);
        while(spair->scheduler->get_func(spair->seed_data, &seed)){
            if((seed.target_pos - row->target_pos)
              > model->max_target_advance)
                break;
            Scheduler_Pair_add_seed(spair, &seed);
            spair->scheduler->next_func(spair->seed_data);
            }
        Scheduler_Pair_reset_rows(spair);
        Scheduler_Row_process(row, spair);
        Lookahead_next(spair->row_index);
    } while(TRUE);
    /* Call model calc exit funcs */
    for(i = 0; i < model->calc_list->len; i++){
        calc = model->calc_list->pdata[i];
        if(calc && calc->exit_func)
            calc->exit_func(spair->region,
                            spair->user_data);
        }
    /* Call model exit funcs */
    if(model->exit_func)
        model->exit_func(spair->region, spair->user_data);
    return;
    }

/**/

static void Scheduler_Lookahead_free_func(gpointer data,
                                          gpointer user_data){
    register Scheduler_Row *row = data;
    register Scheduler_Pair *spair = user_data;
    Scheduler_Row_destroy(row, spair);
    return;
    }

Scheduler_Pair *Scheduler_Pair_create(
                           Scheduler *scheduler,
                           STraceback *straceback,
                           gint query_len, gint target_len,
                           SubOpt *subopt,
                           Boundary *boundary, gint best_seed_id,
                           gpointer seed_data, gpointer user_data){
    register Scheduler_Pair *spair = g_new(Scheduler_Pair, 1);
    spair->scheduler = scheduler;
    spair->slist_set = SListSet_create();
    spair->region = Region_create(0, 0, query_len, target_len);
    spair->seed_data = seed_data;
    spair->user_data = user_data;
    spair->best_seed_id = best_seed_id;
    if(scheduler->has_traceback){
        g_assert(straceback);
        spair->straceback = STraceback_share(straceback);
    } else {
        spair->straceback = NULL;
        }
    spair->row_index = Lookahead_create(
               scheduler->is_forward?0:-target_len,
               scheduler->model->max_target_advance,
               Scheduler_Lookahead_free_func, spair);
    spair->cell_recycle = RecycleBin_create("Scheduler_Cell",
                                            scheduler->cell_size, 1024);
    if(scheduler->use_boundary){
        g_assert(boundary);
        spair->boundary = Boundary_share(boundary);
        spair->span_seed_cache = SparseCache_create(query_len+1,
                                                    Scheduler_Cache_fill_func,
                                                    Scheduler_Cache_empty_func,
                                                    NULL, spair);
        spair->span_history = Scheduler_SpanHistory_create(spair->slist_set,
                                                           scheduler->model);
    } else {
        g_assert(!boundary);
        spair->boundary = NULL;
        spair->span_history = NULL;
        spair->span_seed_cache = NULL;
        }
    spair->subopt_index = SubOpt_Index_create(subopt, spair->region);
    return spair;
    }

void Scheduler_Pair_destroy(Scheduler_Pair *spair){
    if(spair->straceback)
        STraceback_destroy(spair->straceback);
    if(spair->span_history)
        Scheduler_SpanHistory_destroy(spair->span_history);
    if(spair->boundary)
        Boundary_destroy(spair->boundary);
    Region_destroy(spair->region);
    Lookahead_destroy(spair->row_index);
    RecycleBin_destroy(spair->cell_recycle);
    SListSet_destroy(spair->slist_set);
    if(spair->span_seed_cache)
        SparseCache_destroy(spair->span_seed_cache);
    if(spair->subopt_index)
        SubOpt_Index_destroy(spair->subopt_index);
    g_free(spair);
    return;
    }

static void Scheduler_Pair_align_rows(Scheduler_Pair *spair){
    register gint i;
    /* Set other rows in row_index to same position as first row */
    register Scheduler_Row *first_row, *row;
    register Scheduler_Cell *first_cell;
    first_row = Lookahead_get(spair->row_index, 0);
    g_assert(first_row);
    first_cell = Lookahead_get(first_row->cell_index, 0);
    g_assert(first_cell);
    for(i = 1; i <= spair->row_index->max_advance; i++){
        row = Lookahead_get(spair->row_index, i);
        if(row)
            Scheduler_Row_move(row, first_cell->query_pos);
        }
    /* FIXME: optimisation
     *        can do this faster with Lookahead_next()
     */
    return;
    }

static gboolean Scheduler_Pair_is_valid(Scheduler_Pair *spair){
    register gint i;
    register Scheduler_Row *row;
    for(i = 0; i <= spair->row_index->max_advance; i++){
        row = Lookahead_get(spair->row_index, i);
        g_assert((!row) || Scheduler_Row_is_valid(row));
        }
    return TRUE;
    }

/**/

