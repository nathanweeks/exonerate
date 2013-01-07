/****************************************************************\
*                                                                *
*  C4 dynamic programming library - code for models              *
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

#include "c4.h"
#include "matrix.h"

#include <ctype.h>  /* For isalnum() */
#include <string.h> /* For strlen() */

/**/

static C4_State *C4_State_create(gchar *name){
    register C4_State *state = g_new(C4_State, 1);
    g_assert(name);
    state->name = g_strdup(name);
    state->input_transition_list = g_ptr_array_new();
    state->output_transition_list = g_ptr_array_new();
    state->src_shadow_list = g_ptr_array_new();
    return state;
    }

static void C4_State_destroy(C4_State *state){
    g_assert(state);
    g_ptr_array_free(state->input_transition_list, TRUE);
    g_ptr_array_free(state->output_transition_list, TRUE);
    g_ptr_array_free(state->src_shadow_list, TRUE);
    g_free(state->name);
    g_free(state);
    return;
    }

/**/

static C4_Calc *C4_Calc_create(gchar *name, C4_Score max_score,
                      C4_CalcFunc calc_func, gchar *calc_macro,
                      C4_PrepFunc init_func, gchar *init_macro,
                      C4_PrepFunc exit_func, gchar *exit_macro,
                      C4_Protect protect){
    register C4_Calc *calc = g_new(C4_Calc, 1);
    g_assert(name);
    g_assert(calc_macro?(calc_func?TRUE:FALSE):TRUE);
    g_assert(init_macro?(init_func?TRUE:FALSE):TRUE);
    g_assert(exit_macro?(exit_func?TRUE:FALSE):TRUE);
    if(calc_func && (!calc_macro))
        g_warning("Missing calc_macro for Calc: [%s]", name);
    if(init_func && (!init_macro))
        g_warning("Missing init_macro for Calc: [%s]", name);
    if(exit_func && (!exit_macro))
        g_warning("Missing exit_macro for Calc: [%s]", name);
    calc->name = g_strdup(name);
    calc->max_score = max_score;
    calc->calc_func = calc_func;
    calc->calc_macro = g_strdup(calc_macro);
    calc->init_func = init_func;
    calc->init_macro = g_strdup(init_macro);
    calc->exit_func = exit_func;
    calc->exit_macro = g_strdup(exit_macro);
    calc->protect = protect;
    return calc;
    }

static void C4_Calc_destroy(C4_Calc *calc){
    g_assert(calc);
    if(calc->calc_macro)
        g_free(calc->calc_macro);
    if(calc->init_macro)
        g_free(calc->init_macro);
    if(calc->exit_macro)
        g_free(calc->exit_macro);
    g_free(calc->name);
    g_free(calc);
    return;
    }

static gboolean C4_Calc_diff(C4_Calc *a, C4_Calc *b){
    if((a->max_score == b->max_score)
    && (a->calc_func == b->calc_func)
    && (a->init_func == b->init_func)
    && (a->exit_func == b->exit_func)
    && (a->protect == b->protect))
        return FALSE;
    return TRUE;
    }

/**/

static C4_StartState* C4_StartState_create(C4_State *state,
       C4_Scope scope,
       C4_CellStartFunc cell_start_func, gchar *cell_start_macro){
    register C4_StartState *start_state = g_new(C4_StartState, 1);
    g_assert(state);
    start_state->state = state;
    start_state->scope = scope;
    start_state->cell_start_func = cell_start_func;
    if(cell_start_macro)
        start_state->cell_start_macro = g_strdup(cell_start_macro);
    else
        start_state->cell_start_macro = NULL;
    return start_state;
    }

static void C4_StartState_destroy(C4_StartState *start_state){
    g_assert(start_state);
    if(start_state->cell_start_macro)
        g_free(start_state->cell_start_macro);
    g_free(start_state);
    return;
    }

/**/

static C4_EndState* C4_EndState_create(C4_State *state, C4_Scope scope,
        C4_CellEndFunc cell_end_func, gchar *cell_end_macro){
    register C4_EndState *end_state = g_new(C4_EndState, 1);
    g_assert(state);
    end_state->state = state;
    end_state->scope = scope;
    end_state->cell_end_func = cell_end_func;
    if(cell_end_macro)
        end_state->cell_end_macro = g_strdup(cell_end_macro);
    else
        end_state->cell_end_macro = NULL;
    return end_state;
    }

static void C4_EndState_destroy(C4_EndState *end_state){
    g_assert(end_state);
    if(end_state->cell_end_macro)
        g_free(end_state->cell_end_macro);
    g_free(end_state);
    return;
    }

/**/

static C4_Transition *C4_Transition_create(gchar *name,
                         C4_State *input, C4_State *output,
                         gint advance_query, gint advance_target,
                         C4_Calc *calc, C4_Label label,
                         gpointer label_data){
    register C4_Transition *transition = g_new(C4_Transition, 1);
    g_assert(name);
    transition->name = g_strdup(name);
    transition->input = input;
    transition->output = output;
    transition->advance_query = advance_query;
    transition->advance_target = advance_target;
    transition->calc = calc;
    transition->label = label;
    transition->label_data = label_data;
    transition->dst_shadow_list = g_ptr_array_new();
    return transition;
    }

static void C4_Transition_destroy(C4_Transition *transition){
    g_assert(transition);
    g_ptr_array_free(transition->dst_shadow_list, TRUE);
    g_free(transition->name);
    g_free(transition);
    return;
    }

/**/

static C4_Shadow *C4_Shadow_create(gchar *name,
                    C4_State *src, C4_Transition *dst,
                    C4_StartFunc start_func, gchar *start_macro,
                    C4_EndFunc end_func, gchar *end_macro){
    register C4_Shadow *shadow = g_new0(C4_Shadow, 1);
    g_assert(name);
    g_assert(start_macro?(start_func?TRUE:FALSE):TRUE);
    g_assert(end_macro?(end_func?TRUE:FALSE):TRUE);
    if(start_func && (!start_macro))
        g_warning("Missing start_macro for Shadow: [%s]", name);
    if(end_func && (!end_macro))
        g_warning("Missing end_macro for Shadow: [%s]", name);
    shadow->name = g_strdup(name);
    g_assert(src);
    g_assert(dst);
    shadow->src_state_list = g_ptr_array_new();
    shadow->dst_transition_list = g_ptr_array_new();
    g_ptr_array_add(shadow->src_state_list, src);
    g_ptr_array_add(shadow->dst_transition_list, dst);
    shadow->start_func = start_func;
    if(start_macro)
        shadow->start_macro = g_strdup(start_macro);
    shadow->end_func = end_func;
    if(end_macro)
        shadow->end_macro = g_strdup(end_macro);
    g_ptr_array_add(src->src_shadow_list, shadow);
    g_ptr_array_add(dst->dst_shadow_list, shadow);
    return shadow;
    }

static void C4_Shadow_destroy(C4_Shadow *shadow){
    g_ptr_array_free(shadow->src_state_list, TRUE);
    g_ptr_array_free(shadow->dst_transition_list, TRUE);
    g_free(shadow->name);
    if(shadow->start_macro)
        g_free(shadow->start_macro);
    if(shadow->end_macro)
        g_free(shadow->end_macro);
    g_free(shadow);
    return;
    }

void C4_Shadow_add_src_state(C4_Shadow *shadow, C4_State *src){
    g_assert(shadow);
    g_assert(src);
    g_ptr_array_add(shadow->src_state_list, src);
    g_ptr_array_add(src->src_shadow_list, shadow);
    return;
    }

void C4_Shadow_add_dst_transition(C4_Shadow *shadow,
                                  C4_Transition *dst){
    g_assert(shadow);
    g_assert(dst);
    g_ptr_array_add(shadow->dst_transition_list, dst);
    g_ptr_array_add(dst->dst_shadow_list, shadow);
    return;
    }

static gboolean C4_Shadow_is_valid(C4_Shadow *shadow, C4_Model *model){
    register gint i, j;
    register C4_State *state;
    register C4_Transition *transition;
    g_assert(shadow);
    g_assert(shadow->src_state_list->len);
    g_assert(shadow->dst_transition_list->len);
    for(i = 0; i < shadow->src_state_list->len; i++){
        state = shadow->src_state_list->pdata[i];
        g_assert(state);
        for(j = 0; j < shadow->dst_transition_list->len; j++){
            transition = shadow->dst_transition_list->pdata[j];
            g_assert(transition);
            /* FIXME: this is not necessarily true for derived models */
            /*
            g_assert(C4_Model_path_is_possible(model, state,
                                               transition->input));
            */
            }
        }
    return TRUE;
    }

/**/

static C4_Portal *C4_Portal_create(gchar *name, C4_Calc *calc,
           gint advance_query, gint advance_target){
    register C4_Portal *portal = g_new(C4_Portal, 1);
    g_assert(name);
    g_assert(calc);
    portal->name = g_strdup(name);
    portal->calc = calc;
    portal->advance_query = advance_query;
    portal->advance_target = advance_target;
    portal->transition_list = g_ptr_array_new();
    return portal;
    }

static void C4_Portal_destroy(C4_Portal *portal){
    g_free(portal->name);
    g_ptr_array_free(portal->transition_list, TRUE);
    g_free(portal);
    return;
    }

/**/

static void C4_Span_find_loop_transitions(C4_Span *span){
    register gint i;
    register C4_Transition *transition;
    span->query_loop = span->target_loop = NULL;
    for(i = 0; i < span->span_state->output_transition_list->len; i++){
        transition = span->span_state->output_transition_list->pdata[i];
        if(transition->output == span->span_state){ /* looping */
            /* Find a C4_Calc which is sequence position independent */
            if((!transition->calc) || (!transition->calc->calc_func)){
                g_assert((!transition->calc)
                      || (!transition->calc->max_score));
                g_assert(transition->advance_query
                      || transition->advance_target);
                g_assert(!(transition->advance_query
                        && transition->advance_target));
                if(transition->advance_query){
                    g_assert(!span->query_loop);
                    g_assert(!transition->advance_target);
                    span->query_loop = transition;
                } else {
                    g_assert(!span->target_loop);
                    g_assert(transition->advance_target);
                    g_assert(!transition->advance_query);
                    span->target_loop = transition;
                    }
                }
            }
        }
    g_assert(span->query_loop || span->target_loop);
    return;
    }

static C4_Span *C4_Span_create(gchar *name, C4_State *span_state,
                               gint min_query, gint max_query,
                               gint min_target, gint max_target){
    register C4_Span *span = g_new(C4_Span, 1);
    g_assert(name);
    g_assert(span_state);
    g_assert(min_query >= 0);
    g_assert(min_query <= max_query);
    g_assert(min_target >= 0);
    g_assert(min_target <= max_target);
    span->name = g_strdup(name);
    span->span_state = span_state;
    span->min_query = min_query;
    span->max_query = max_query;
    span->min_target = min_target;
    span->max_target = max_target;
    C4_Span_find_loop_transitions(span);
    return span;
    }

static void C4_Span_destroy(C4_Span *span){
    g_assert(span);
    g_free(span->name);
    g_free(span);
    return;
    }

/**/

C4_Calc *C4_Model_add_calc(C4_Model *model, gchar *name,
                      C4_Score max_score,
                      C4_CalcFunc calc_func, gchar *calc_macro,
                      C4_PrepFunc init_func, gchar *init_macro,
                      C4_PrepFunc exit_func, gchar *exit_macro,
                      C4_Protect protect){
    register C4_Calc *calc = C4_Calc_create(name, max_score,
                                            calc_func, calc_macro,
                                            init_func, init_macro,
                                            exit_func, exit_macro,
                                            protect);
    g_assert(model);
    g_assert(name);
    g_assert(model->is_open);
    g_ptr_array_add(model->calc_list, calc);
    return calc;
    }

static C4_Calc *C4_Model_match_calc(C4_Model *model, C4_Calc *calc){
    register C4_Calc *mcalc = NULL;
    register gint i;
    g_assert(model);
    g_assert(calc);
    for(i = 0; i < model->calc_list->len; i++){
        mcalc = model->calc_list->pdata[i];
        if(!C4_Calc_diff(mcalc, calc))
            return mcalc;
        }
    return NULL;
    }

C4_State *C4_Model_add_state(C4_Model *model, gchar *name){
    register C4_State *state = C4_State_create(name);
    g_assert(model);
    g_assert(name);
    g_assert(model->is_open);
    g_ptr_array_add(model->state_list, state);
    return state;
    }

GPtrArray *C4_Model_select_transitions(C4_Model *model, C4_Label label){
    register gint i;
    register GPtrArray *list = g_ptr_array_new();
    register C4_Transition *transition;
    for(i = 0; i < model->transition_list->len; i++){
        transition = model->transition_list->pdata[i];
        if(transition->label == label)
            g_ptr_array_add(list, transition);
        }
    if(!list->len){
        g_ptr_array_free(list, TRUE);
        return NULL;
        }
    return list;
    }

C4_Transition *C4_Model_select_single_transition(C4_Model *model,
                                                 C4_Label label){
    register C4_Transition *transition = NULL;
    register GPtrArray *transition_list
        = C4_Model_select_transitions(model, label);
    g_assert(transition_list);
    g_assert(transition_list->len == 1);
    if(transition_list){
        if(transition_list->len == 1)
            transition = transition_list->pdata[0];
        g_ptr_array_free(transition_list, TRUE);
        }
    return transition;
    }


/**/

C4_Transition *C4_Model_add_transition(C4_Model *model, gchar *name,
                       C4_State *input, C4_State *output,
                       gint advance_query, gint advance_target,
                       C4_Calc *calc, C4_Label label,
                       gpointer label_data){
    register C4_Transition *transition;
    g_assert(model);
    g_assert(name);
    g_assert(model->is_open);
    g_assert(advance_query >= 0);
    g_assert(advance_target >= 0);
    g_assert(  (label == C4_Label_NONE)
            || (advance_query || advance_target));
    g_assert(  (label != C4_Label_MATCH)
            || (advance_query && advance_target));
    g_assert(  (label != C4_Label_GAP)
            || (advance_query || advance_target));
    g_assert(  (label != C4_Label_FRAMESHIFT)
            || (advance_query || advance_target));
    if(!input)
        input = model->start_state->state;
    if(!output)
        output = model->end_state->state;
    transition = C4_Transition_create(name, input, output,
                    advance_query, advance_target, calc, label,
                    label_data);
    g_ptr_array_add(input->output_transition_list, transition);
    g_ptr_array_add(output->input_transition_list, transition);
    g_ptr_array_add(model->transition_list, transition);
    return transition;
    }

C4_Shadow *C4_Model_add_shadow(C4_Model *model, gchar *name,
                    C4_State *src, C4_Transition *dst,
                    C4_StartFunc start_func, gchar *start_macro,
                    C4_EndFunc end_func, gchar *end_macro){
    register C4_Shadow *shadow;
    register gint i;
    g_assert(model);
    g_assert(name);
    g_assert(model->is_open);
    g_assert(start_func);
    g_assert(end_func);
    if(!src)
        src = model->start_state->state;
    if(dst){
        shadow = C4_Shadow_create(name, src, dst,
                                  start_func, start_macro,
                                  end_func, end_macro);
    } else { /* Use all transitions to END as dst */
        g_assert(model->end_state->state->input_transition_list->len);
        dst = model->end_state->state->input_transition_list->pdata[0];
        shadow = C4_Shadow_create(name, src, dst,
                                  start_func, start_macro,
                                  end_func, end_macro);
        for(i = 1; i < model->end_state->state
                     ->input_transition_list->len; i++){
            dst = model->end_state->state
                ->input_transition_list->pdata[i];
            C4_Shadow_add_dst_transition(shadow, dst);
            }
        }
    g_assert(shadow);
    g_ptr_array_add(model->shadow_list, shadow);
    return shadow;
    }

C4_Portal *C4_Model_add_portal(C4_Model *model, gchar *name,
           C4_Calc *calc, gint advance_query, gint advance_target){
    register C4_Portal *portal;
    g_assert(model->is_open);
    g_assert(calc);
    g_assert(calc->calc_func); /* Must have position specific calc */
    portal = C4_Portal_create(name, calc,
                              advance_query, advance_target);
    /* Add the portal to the model */
    g_ptr_array_add(model->portal_list, portal);
    return portal;
    }

C4_Span *C4_Model_add_span(C4_Model *model, gchar *name,
                           C4_State *span_state,
                           gint min_query, gint max_query,
                           gint min_target, gint max_target){
    register C4_Span *span;
    g_assert(model->is_open);
    span = C4_Span_create(name, span_state, min_query, max_query,
                                            min_target, max_target);
    g_ptr_array_add(model->span_list, span);
    return span;
    }

C4_Model *C4_Model_create(gchar *name){
    register C4_Model *model = g_new(C4_Model, 1);
    g_assert(name);
    model->name = g_strdup(name);
    /**/
    model->global_code_list = g_ptr_array_new();
    model->local_code_list = g_ptr_array_new();
    model->cflags_add_list = g_ptr_array_new();
    model->init_func = NULL;
    model->init_macro = NULL;
    model->exit_func = NULL;
    model->exit_macro = NULL;
    /**/
    model->thread_ref = ThreadRef_create();
    model->is_open = TRUE;
    /**/
    model->state_list      = g_ptr_array_new();
    model->transition_list = g_ptr_array_new();
    model->shadow_list     = g_ptr_array_new();
    model->calc_list       = g_ptr_array_new();
    model->portal_list     = g_ptr_array_new();
    model->span_list       = g_ptr_array_new();
    /**/
    model->start_state = C4_StartState_create(
            C4_Model_add_state(model, "START"),
            C4_Scope_ANYWHERE, NULL, NULL);
    model->end_state = C4_EndState_create(
            C4_Model_add_state(model, "END"),
            C4_Scope_ANYWHERE, NULL, NULL);
    return model;
    }

void C4_Model_rename(C4_Model *model, gchar *name){
    g_assert(model);
    g_assert(name);
    g_free(model->name);
    model->name = g_strdup(name);
    return;
    }

/**/

static void C4_Model_destroy_string_list(GPtrArray *string_list){
    register gint i;
    for(i = 0; i < string_list->len; i++)
        g_free(string_list->pdata[i]);
    g_ptr_array_free(string_list, TRUE);
    return;
    }

void C4_Model_destroy(C4_Model *model){
    register gint i;
    g_assert(model);
    if(ThreadRef_destroy(model->thread_ref))
        return;
    g_free(model->name);
    C4_Model_destroy_string_list(model->global_code_list);
    C4_Model_destroy_string_list(model->local_code_list);
    C4_Model_destroy_string_list(model->cflags_add_list);
    if(model->init_macro)
        g_free(model->init_macro);
    if(model->exit_macro)
        g_free(model->exit_macro);
    C4_StartState_destroy(model->start_state);
    C4_EndState_destroy(model->end_state);
    for(i = 0; i < model->state_list->len; i++)
        C4_State_destroy(model->state_list->pdata[i]);
    for(i = 0; i < model->transition_list->len; i++)
        C4_Transition_destroy(model->transition_list->pdata[i]);
    for(i = 0; i < model->shadow_list->len; i++)
        C4_Shadow_destroy(model->shadow_list->pdata[i]);
    for(i = 0; i < model->calc_list->len; i++)
        C4_Calc_destroy(model->calc_list->pdata[i]);
    for(i = 0; i < model->portal_list->len; i++)
        C4_Portal_destroy(model->portal_list->pdata[i]);
    for(i = 0; i < model->span_list->len; i++)
        C4_Span_destroy(model->span_list->pdata[i]);
    g_ptr_array_free(model->state_list, TRUE);
    g_ptr_array_free(model->transition_list, TRUE);
    g_ptr_array_free(model->shadow_list, TRUE);
    g_ptr_array_free(model->calc_list, TRUE);
    g_ptr_array_free(model->portal_list, TRUE);
    g_ptr_array_free(model->span_list, TRUE);
    g_free(model);
    return;
    }

C4_Model *C4_Model_share(C4_Model *model){
    g_assert(model);
    g_assert(!model->is_open);
    ThreadRef_share(model->thread_ref);
    return model;
    }

static gint C4_State_compare_by_name(gconstpointer a,
                                     gconstpointer b){
    register C4_State *state_a = (C4_State*)a,
                      *state_b = (C4_State*)b;
    return strcmp(state_a->name, state_b->name);
    }

static gint C4_Transition_compare_by_name(gconstpointer a,
                                          gconstpointer b){
    register C4_Transition *tran_a = (C4_Transition*)a,
                           *tran_b = (C4_Transition*)b;
    return strcmp(tran_a->name, tran_b->name);
    }

/* FIXME: rewrite without using state names */
C4_State **C4_Model_build_state_map(C4_Model *src, C4_Model *dst){
    register C4_State **state_map = g_new0(C4_State*,
                                           src->state_list->len);
    register gint i;
    register C4_State *src_state, *dst_state;
    register GTree *state_map_tree
        = g_tree_new(C4_State_compare_by_name);
    g_assert(src);
    g_assert(dst);
    g_assert(src->state_list->len == dst->state_list->len);
    for(i = 0; i < dst->state_list->len; i++){
        dst_state = dst->state_list->pdata[i];
        if(g_tree_lookup(state_map_tree, dst_state))
            g_error("Duplicate state name [%s] in model [%s]",
                    dst_state->name, src->name);
        g_assert(!g_tree_lookup(state_map_tree, dst_state));
        g_tree_insert(state_map_tree, dst_state, dst_state);
        }
    for(i = 0; i < src->state_list->len; i++){
        src_state = src->state_list->pdata[i];
        dst_state = g_tree_lookup(state_map_tree, src_state);
        g_assert(dst_state);
        state_map[src_state->id] = dst_state;
        }
    g_tree_destroy(state_map_tree);
    return state_map;
    }

/* FIXME: rewrite without using transition names */
C4_Transition **C4_Model_build_transition_map(C4_Model *src,
                                              C4_Model *dst){
    register C4_Transition **transition_map = g_new0(C4_Transition*,
                                      src->transition_list->len);
    register gint i;
    register C4_Transition *src_transition, *dst_transition;
    register GTree *transition_map_tree
        = g_tree_new(C4_Transition_compare_by_name);
    g_assert(src);
    g_assert(dst);
    g_assert(src->transition_list->len == dst->transition_list->len);
    for(i = 0; i < dst->transition_list->len; i++){
        dst_transition = dst->transition_list->pdata[i];
        if(g_tree_lookup(transition_map_tree, dst_transition))
            g_error("Duplicate transition name [%s] in model [%s]",
                    dst_transition->name, src->name);
        g_assert(!g_tree_lookup(transition_map_tree, dst_transition));
        g_tree_insert(transition_map_tree, dst_transition,
                                           dst_transition);
        }
    for(i = 0; i < src->transition_list->len; i++){
        src_transition = src->transition_list->pdata[i];
        dst_transition = g_tree_lookup(transition_map_tree,
                                       src_transition);
        g_assert(dst_transition);
        transition_map[src_transition->id] = dst_transition;
        }
    g_tree_destroy(transition_map_tree);
    return transition_map;
    }

/**/

void C4_Model_make_stereo(C4_Model *model, gchar *suffix_a,
                                           gchar *suffix_b){
    register gint i, j,
                  prev_state_count = model->state_list->len,
                  prev_transition_count = model->transition_list->len,
                  prev_shadow_count = model->shadow_list->len;
    register C4_State *state,
                     **state_map = g_new0(C4_State*, prev_state_count);
    register C4_Transition *transition,
                     **transition_map = g_new0(C4_Transition*,
                                               prev_transition_count);
    register C4_Shadow *shadow, *new_shadow;
    register gchar *name;
    g_assert(model->is_open);
    /* Copy each state -> state.suffix_b */
    for(i = 0; i < prev_state_count; i++){
        state = model->state_list->pdata[i];
        if((state != model->start_state->state)
        && (state != model->end_state->state)){
            name = g_strconcat(state->name, " ", suffix_b, NULL);
            state_map[state->id] = C4_Model_add_state(model, name);
            g_free(name);
            }
        }
    /* Copy each transition -> transition.suffix_b */
    for(i = 0; i < prev_transition_count; i++){
        transition = model->transition_list->pdata[i];
        name = g_strconcat(transition->name, " ", suffix_b, NULL);
        transition_map[transition->id] = C4_Model_add_transition(
                model, name,
                state_map[transition->input->id],
                state_map[transition->output->id],
                transition->advance_query, transition->advance_target,
                transition->calc, transition->label,
                transition->label_data);
        g_free(name);
        }
    /* Copy each shadow -> shadow.suffix_b */
    for(i = 0; i < prev_shadow_count; i++){
        shadow = model->shadow_list->pdata[i];
        name = g_strconcat(shadow->name, " ", suffix_b, NULL);
        g_assert(shadow->src_state_list->len);
        g_assert(shadow->dst_transition_list->len);
        state = shadow->src_state_list->pdata[0];
        transition = shadow->dst_transition_list->pdata[0];
        new_shadow = C4_Model_add_shadow(model, name,
                state_map[state->id],
                transition_map[transition->id],
                shadow->start_func, shadow->start_macro,
                shadow->end_func, shadow->end_macro);
        g_free(name);
        /* Copy remainding src states */
        for(j = 1; j < shadow->src_state_list->len; j++){
            state = shadow->src_state_list->pdata[j];
            C4_Shadow_add_src_state(new_shadow, state);
            }
        /* Copy remainding dst transitions */
        for(j = 1; j < shadow->dst_transition_list->len; j++){
            transition = shadow->dst_transition_list->pdata[j];
            C4_Shadow_add_dst_transition(new_shadow, transition);
            }
        }
    /* Rename each original state -> state.suffix_a */
    for(i = 0; i < prev_state_count; i++){
        state = model->state_list->pdata[i];
        if((state != model->start_state->state)
        && (state != model->end_state->state)){
            name = state->name;
            state->name = g_strconcat(name, " ", suffix_a, NULL);
            g_free(name);
            }
        }
    /* Rename each transition -> transition.suffix_a */
    for(i = 0; i < prev_transition_count; i++){
        transition = model->transition_list->pdata[i];
        name = transition->name;
        transition->name = g_strconcat(name, " ", suffix_a, NULL);
        g_free(name);
        }
    /* Rename each shadow -> shadow.suffix_a */
    for(i = 0; i < prev_shadow_count; i++){
        shadow = model->shadow_list->pdata[i];
        name = shadow->name;
        shadow->name = g_strconcat(name, " ", suffix_a, NULL);
        g_free(name);
        }
    g_free(state_map);
    g_free(transition_map);
    return;
    }

static void C4_Model_insert_calcs(C4_Model *target, C4_Model *insert,
                                  C4_Calc **calc_map){
    register gint i;
    register C4_Calc *insert_calc, *target_calc;
    for(i = 0; i < insert->calc_list->len; i++){
        insert_calc = insert->calc_list->pdata[i];
        target_calc = C4_Model_match_calc(target, insert_calc);
        if(!target_calc)
            target_calc = C4_Model_add_calc(target, insert_calc->name,
                insert_calc->max_score,
                insert_calc->calc_func, insert_calc->calc_macro,
                insert_calc->init_func, insert_calc->init_macro,
                insert_calc->exit_func, insert_calc->exit_macro,
                insert_calc->protect);
        calc_map[insert_calc->id] = target_calc;
        }
    return;
    }

static void C4_Model_insert_states(C4_Model *target, C4_Model *insert,
                                   C4_State **state_map){
    register gint i;
    register C4_State *insert_state;
    for(i = 0; i < insert->state_list->len; i++){
        insert_state = insert->state_list->pdata[i];
        if((insert_state != insert->start_state->state)
        && (insert_state != insert->end_state->state))
            state_map[insert_state->id]
                = C4_Model_add_state(target, insert_state->name);
        }
    return;
    }

static void C4_Model_insert_transitions(C4_Model *target,
                                        C4_Model *insert,
                C4_State *src, C4_State *dst,
                C4_Calc **calc_map,
                C4_State **state_map,
                C4_Transition **transition_map){
    register gint i;
    register C4_Transition *transition;
    register C4_Calc *calc;
    for(i = 0; i < insert->transition_list->len; i++){
        transition = insert->transition_list->pdata[i];
        g_assert(transition);
        if(transition->calc){
            calc = calc_map[transition->calc->id];
        } else {
            calc = NULL;
            }
        transition_map[transition->id] = C4_Model_add_transition(
                  target, transition->name,
                  state_map[transition->input->id],
                  state_map[transition->output->id],
                  transition->advance_query, transition->advance_target,
                  calc, transition->label, transition->label_data);
        }
    return;
    }

static void C4_Model_insert_shadows(C4_Model *target,
                                    C4_Model *insert,
                                    C4_State *src, C4_State *dst,
                                    C4_State **state_map,
                                    C4_Transition **transition_map){
    register gint i, j;
    register C4_Shadow *shadow, *new_shadow;
    register C4_State *state, *src_state;
    register C4_Transition *transition, *dst_transition;
    /* For each insert shadow */
    for(i = 0; i < insert->shadow_list->len; i++){
        shadow = insert->shadow_list->pdata[i];
        g_assert(shadow);
        g_assert(shadow->src_state_list->len);
        g_assert(shadow->dst_transition_list->len);
        state = shadow->src_state_list->pdata[0];
        src_state = state_map[state->id];
        transition = shadow->dst_transition_list->pdata[0];
        dst_transition = transition_map[transition->id];
        new_shadow = C4_Model_add_shadow(target, shadow->name,
              src_state, dst_transition,
              shadow->start_func, shadow->start_macro,
              shadow->end_func, shadow->end_macro);
        for(j = 1; j < shadow->src_state_list->len; j++){
            state = shadow->src_state_list->pdata[j];
            src_state = state_map[state->id];
            C4_Shadow_add_src_state(new_shadow, src_state);
            }
        for(j = 1; j < shadow->dst_transition_list->len; j++){
            transition = shadow->dst_transition_list->pdata[j];
            dst_transition = transition_map[transition->id];
            C4_Shadow_add_dst_transition(new_shadow, dst_transition);
            }
        }
    return;
    }

static gboolean C4_Calc_is_same(C4_Calc *calc_a, C4_Calc *calc_b){
    if((calc_a->max_score == calc_b->max_score)
    && (calc_a->calc_func == calc_b->calc_func)
    && (calc_a->init_func == calc_b->init_func)
    && (calc_a->exit_func == calc_b->exit_func)
    && (calc_a->protect == calc_b->protect))
        return TRUE;
    return FALSE;
    }

static gboolean C4_Portal_is_same(C4_Portal *portal_a,
                                  C4_Portal *portal_b){
    if((portal_a->advance_query == portal_b->advance_query)
    && (portal_a->advance_target == portal_b->advance_target)
    && C4_Calc_is_same(portal_a->calc, portal_b->calc))
        return TRUE;
    return FALSE;
    }

static void C4_Model_insert_portals(C4_Model *target,
            C4_Model *insert, C4_Calc **calc_map){
    register gint i, j;
    register C4_Portal *insert_portal, *target_portal;
    register C4_Transition *transition;
    register gboolean found_same;
    for(i = 0; i < insert->portal_list->len; i++){
        insert_portal = insert->portal_list->pdata[i];
        found_same = FALSE;
        for(j = 0; j < target->portal_list->len; j++){
            target_portal = target->portal_list->pdata[j];
            if(C4_Portal_is_same(insert_portal, target_portal)){
                for(j = 0; j < insert_portal->transition_list->len;
                    j++){
                    transition
                        = insert_portal->transition_list->pdata[j];
                    g_ptr_array_add(target_portal->transition_list,
                                    transition);
                    }
                found_same = TRUE;
                break;
                }
            }
        if(!found_same){
            C4_Model_add_portal(target, insert_portal->name,
                calc_map[insert_portal->calc->id],
                insert_portal->advance_query,
                insert_portal->advance_target);
            }
        }
    return;
    }
/* FIXME: merge duplicate portals */

static void C4_Model_insert_spans(C4_Model *target,
            C4_Model *insert, C4_State **state_map){
    register gint i;
    register C4_Span *span;
    for(i = 0; i < insert->span_list->len; i++){
        span = insert->span_list->pdata[i];
        C4_Model_add_span(target, span->name,
                          state_map[span->span_state->id],
                          span->min_query, span->max_query,
                          span->min_target, span->max_target);
        }
    return;
    }

static gboolean C4_Model_string_list_contains(GPtrArray *list,
                                              gchar *str){
    register gint i;
    for(i = 0; i < list->len; i++)
        if(!strcmp(list->pdata[i], str))
            return TRUE;
    return FALSE;
    }

static void C4_Model_merge_codegen_list(GPtrArray *dst, GPtrArray *src){
    register gint i;
    for(i = 0; i < src->len; i++)
        if(!C4_Model_string_list_contains(dst, src->pdata[i]))
            g_ptr_array_add(dst, g_strdup(src->pdata[i]));
    return;
    }

static void C4_Model_merge_codegen(C4_Model *dst, C4_Model *src){
    C4_Model_merge_codegen_list(dst->global_code_list,
                                src->global_code_list);
    C4_Model_merge_codegen_list(dst->local_code_list,
                                src->local_code_list);
    C4_Model_merge_codegen_list(dst->cflags_add_list,
                                src->cflags_add_list);
    return;
    }

void C4_Model_insert(C4_Model *target, C4_Model *insert,
                     C4_State *src, C4_State *dst){
    register C4_Calc **calc_map = g_new0(C4_Calc*,
                                         insert->calc_list->len);
    register C4_State **state_map = g_new0(C4_State*,
                                         insert->state_list->len);
    register C4_Transition **transition_map = g_new0(C4_Transition*,
                                         insert->transition_list->len);
    g_assert(target->is_open);
    g_assert(!insert->is_open);
    if(!src)
        src = target->start_state->state;
    if(!dst)
        dst = target->end_state->state;
    C4_Model_insert_calcs(target, insert, calc_map);
    C4_Model_insert_states(target, insert, state_map);
    /* Add src and dst states to the state map */
    state_map[insert->start_state->state->id] = src;
    state_map[insert->end_state->state->id] = dst;
    C4_Model_insert_transitions(target, insert, src, dst,
                                calc_map, state_map, transition_map);
    C4_Model_insert_shadows(target, insert, src, dst,
                            state_map, transition_map);
    C4_Model_insert_portals(target, insert, calc_map);
    C4_Model_insert_spans(target, insert, state_map);
    C4_Model_merge_codegen(target, insert);
    C4_Model_configure_extra(target,
                             insert->init_func, insert->init_macro,
                             insert->exit_func, insert->exit_macro);
    g_free(state_map);
    g_free(transition_map);
    g_free(calc_map);
    return;
    }

/**/

void C4_Model_clear_codegen(C4_Model *model){
    C4_Model_destroy_string_list(model->global_code_list);
    model->global_code_list = g_ptr_array_new();
    C4_Model_destroy_string_list(model->local_code_list);
    model->local_code_list = g_ptr_array_new();
    C4_Model_destroy_string_list(model->cflags_add_list);
    model->cflags_add_list = g_ptr_array_new();
    return;
    }

void C4_Model_append_codegen(C4_Model *model,
                             gchar *global_code, gchar *local_code,
                             gchar *cflags_add){
    if((global_code) && (!C4_Model_string_list_contains(
                model->global_code_list, global_code)))
        g_ptr_array_add(model->global_code_list, g_strdup(global_code));
    if((local_code) && (!C4_Model_string_list_contains(
                model->local_code_list, local_code)))
        g_ptr_array_add(model->local_code_list, g_strdup(local_code));
    if((cflags_add) && (!C4_Model_string_list_contains(
                model->cflags_add_list, cflags_add)))
        g_ptr_array_add(model->cflags_add_list, g_strdup(cflags_add));
    return;
    }

void C4_Model_configure_extra(C4_Model *model,
                              C4_PrepFunc init_func, gchar *init_macro,
                              C4_PrepFunc exit_func, gchar *exit_macro){
    g_assert(model);
    g_assert(model->is_open);
    g_assert(init_macro?(init_func?TRUE:FALSE):TRUE);
    g_assert(exit_macro?(exit_func?TRUE:FALSE):TRUE);
    /* FIXME: all instances should be on a calc ... */
    /*
    g_assert(!init_func);
    g_assert(!exit_func);
    */
    /*
    if(init_func && (!init_macro))
        g_warning("Missing init_macro for Model: [%s]", model->name);
    if(exit_func && (!exit_macro))
        g_warning("Missing exit_macro for Model: [%s]", model->name);
    */
    model->init_func = init_func;
    if(model->init_macro){
        if(model->init_macro)
            g_free(model->init_macro);
        model->init_macro = g_strdup(init_macro);
        }
    model->exit_func = exit_func;
    if(model->exit_macro){
        if(model->exit_macro)
            g_free(model->exit_macro);
        model->exit_macro = g_strdup(exit_macro);
        }
    return;
    }
/* FIXME: need to allow joining of multiple funcs and macros ... */

void C4_Model_configure_start_state(C4_Model *model, C4_Scope scope,
        C4_CellStartFunc cell_start_func, gchar *cell_start_macro){
    g_assert(model);
    /*
    if(cell_start_func && (!cell_start_macro))
        g_warning("Missing cell_start_macro for Model: [%s]",
                  model->name);
    */
    model->start_state->scope = scope;
    model->start_state->cell_start_func = cell_start_func;
    if(model->start_state->cell_start_macro != cell_start_macro){
        if(model->start_state->cell_start_macro)
            g_free(model->start_state->cell_start_macro);
        if(cell_start_macro)
            model->start_state->cell_start_macro
                                     = g_strdup(cell_start_macro);
        else
            model->start_state->cell_start_macro = NULL;
    } else {
        model->start_state->cell_start_macro = NULL;
        }
    return;
    }

void C4_Model_configure_end_state(C4_Model *model, C4_Scope scope,
        C4_CellEndFunc cell_end_func, gchar *cell_end_macro){
    g_assert(model);
    /*
    if(cell_end_func && (!cell_end_macro))
        g_warning("Missing cell_end_macro for Model: [%s]",
                  model->name);
    */
    model->end_state->scope = scope;
    model->end_state->cell_end_func = cell_end_func;
    if(model->end_state->cell_end_macro != cell_end_macro){
        if(model->end_state->cell_end_macro)
            g_free(model->end_state->cell_end_macro);
        if(cell_end_macro)
            model->end_state->cell_end_macro = g_strdup(cell_end_macro);
        else
            model->end_state->cell_end_macro = NULL;
    } else {
        model->end_state->cell_end_macro = NULL;
        }
    return;
    }

/**/

gchar *C4_Scope_get_name(C4_Scope scope){
    register gchar *name = NULL;
    switch(scope){
        case C4_Scope_ANYWHERE:
            name = "anywhere";
            break;
        case C4_Scope_EDGE:
            name = "edge";
            break;
        case C4_Scope_QUERY:
            name = "query";
            break;
        case C4_Scope_TARGET:
            name = "target";
            break;
        case C4_Scope_CORNER:
            name = "corner";
            break;
        default:
            g_error("Unknown C4_Scope type [%d]", scope);
            break;
        }
    return name;
    }

gchar *C4_Label_get_name(C4_Label label){
    register gchar *name = NULL;
    switch(label){
        case C4_Label_NONE:
            name = "none";
            break;
        case C4_Label_MATCH:
            name = "match";
            break;
        case C4_Label_GAP:
            name = "gap";
            break;
        case C4_Label_NER:
            name = "ner";
            break;
        case C4_Label_5SS:
            name = "5'ss";
            break;
        case C4_Label_3SS:
            name = "3'ss";
            break;
        case C4_Label_INTRON:
            name = "intron";
            break;
        case C4_Label_SPLIT_CODON:
            name = "split codon";
            break;
        case C4_Label_FRAMESHIFT:
            name = "frameshift";
            break;
        default:
            g_error("Unknown label [%d]", label);
            break;
        }
    return name;
    }

void C4_Model_print(C4_Model *model){
    register gint i;
    register C4_State *state;
    register C4_Transition *transition;
    register C4_Calc *calc;
    register C4_Portal *portal;
    register C4_Span *span;
    g_assert(model);
    g_print("--\n"
            "    Info for model [%s]\n"
            "            status [%s]\n"
            "       start scope [%s]\n"
            "         end scope [%s]\n",
            model->name, model->is_open?"OPEN":"CLOSED",
            C4_Scope_get_name(model->start_state->scope),
            C4_Scope_get_name(model->end_state->scope));
    for(i = 0; i < model->state_list->len; i++){
        state = model->state_list->pdata[i];
        g_print("State [%s]\n", state->name);
        }
    for(i = 0; i < model->transition_list->len; i++){
        transition = model->transition_list->pdata[i];
        g_print("Transition [%s] ( [%s]->[%s] ) [%d,%d]\n",
                transition->name,
                transition->input->name,
                transition->output->name,
                transition->advance_query,
                transition->advance_query);
        }
    for(i = 0; i < model->calc_list->len; i++){
        calc = model->calc_list->pdata[i];
        g_print("Calc [%s]\n", calc->name);
        }
    for(i = 0; i < model->portal_list->len; i++){
        portal = model->portal_list->pdata[i];
        g_print("Portal [%s]\n", portal->name);
        }
    for(i = 0; i < model->span_list->len; i++){
        span = model->span_list->pdata[i];
        g_print("Span [%s]\n", span->name);
        }
    g_print("--\n");
    return;
    }

static gchar *clean_string(gchar *str){
    register gchar *new_string = g_strdup(str);
    register gint i;
    for(i = 0; new_string[i]; i++){
        if(!isalnum(new_string[i]))
            new_string[i] = '_';
        }
    return new_string;
    }

static gchar *break_string(gchar *str){
    register GString *s = g_string_sized_new(strlen(str)+10);
    register gchar *new_string;
    register gint i;
    for(i = 0; str[i]; i++){
        if(isspace(str[i]))
            g_string_append(s, "\\n");
        else
            g_string_append_c(s, str[i]);
        }
    new_string = s->str;
    g_string_free(s, FALSE);
    return new_string;
    }

void C4_Model_dump_graphviz(C4_Model *model){
    register gint i, j, k;
    register C4_State *state;
    register C4_Transition *transition;
    register C4_Shadow *shadow;
    register gchar *identifier;
    g_assert(model);
    g_assert(!model->is_open);
    identifier = clean_string(model->name);
    g_print("// --- START OF GRAPHVIZ DUMP ---\n");
    g_print("digraph %s {\n", identifier);
    g_free(identifier);
    g_print("// states\n");
    for(i = 0; i < model->state_list->len; i++){
        state = model->state_list->pdata[i];
        identifier = break_string(state->name);
        g_print("state_%d [label=\"%s\"", state->id, identifier);
        g_free(identifier);
        if((state == model->start_state->state)
        || (state == model->end_state->state))
            g_print(",shape=box");
        g_print("]\n");
        }
    g_print("// transitions\n");
    for(i = 0; i < model->transition_list->len; i++){
        transition = model->transition_list->pdata[i];
        g_print("state_%d -> state_%d",
                transition->input->id,
                transition->output->id);
        if(transition->advance_query || transition->advance_target)
            g_print(" [label=\"%d:%d\"]",
              transition->advance_query, transition->advance_target);
        g_print(" // %s", transition->name);
        g_print("\n");
        }
    if(model->shadow_list->len)
        g_print("// shadows\n");
    /* Shadows are just drawn from src to dst->input */
    for(i = 0; i < model->shadow_list->len; i++){
        shadow = model->shadow_list->pdata[i];
        for(j = 0; j < shadow->src_state_list->len; j++){
            state = shadow->src_state_list->pdata[j];
            for(k = 0; k < shadow->dst_transition_list->len; k++){
                transition = shadow->dst_transition_list->pdata[k];
                g_print("state_%d -> state_%d"
                        " [style=dotted,label=\"%s\"]\n",
                        state->id,
                        transition->input->id,
                        transition->name);
                }
            }
        }
    g_print("}\n");
    g_print("// --- END OF GRAPHVIZ DUMP ---\n");
    return;
    }
/* FIXME: include display of portals and spans
 *        with peripheries=2 etc
 */

void C4_Model_open(C4_Model *model){
    g_assert(model);
    g_assert(!model->is_open);
    model->is_open = TRUE;
    return;
    }

static gboolean C4_Model_path_is_possible_recur(C4_Model *model,
                C4_State *src, C4_State *dst, gboolean *visited){
    register C4_Transition *transition;
    register C4_State *next_state;
    register gint i;
    visited[src->id] = TRUE;
    for(i = 0; i < src->output_transition_list->len; i++){
        transition = src->output_transition_list->pdata[i];
        next_state = transition->output;
        if(next_state == dst)
            return TRUE;
        if(!visited[next_state->id]){
            if(C4_Model_path_is_possible_recur(model,
                                    next_state, dst, visited)){
                return TRUE;
                }
            }
        }
    return FALSE;
    }
/* Assumes that state->id is set */

gboolean C4_Model_path_is_possible(C4_Model *model,
                                   C4_State *src, C4_State *dst){
    register gboolean *visited;
    register gboolean is_valid;
    g_assert(model);
    g_assert(src);
    g_assert(dst);
    visited = g_new0(gboolean, model->state_list->len);
    is_valid = C4_Model_path_is_possible_recur(model, src, dst,
                                               visited);
    g_free(visited);
    return is_valid;
    }

static gboolean C4_Model_has_start_to_end_path(C4_Model *model){
    return C4_Model_path_is_possible(model,
                  model->start_state->state,
                  model->end_state->state);
    }

static void C4_Model_set_ids(C4_Model *model){
    register gint i;
    register C4_State *state;
    register C4_Transition *transition;
    register C4_Shadow *shadow;
    register C4_Calc *calc;
    register C4_Portal *portal;
    register C4_Span *span;
    g_assert(model);
    for(i = 0; i < model->state_list->len; i++){
        state = model->state_list->pdata[i];
        state->id = i;
        }
    for(i = 0; i < model->transition_list->len; i++){
        transition = model->transition_list->pdata[i];
        transition->id = i;
        }
    for(i = 0; i < model->shadow_list->len; i++){
        shadow = model->shadow_list->pdata[i];
        shadow->id = i;
        }
    for(i = 0; i < model->calc_list->len; i++){
        calc = model->calc_list->pdata[i];
        calc->id = i;
        }
    for(i = 0; i < model->portal_list->len; i++){
        portal = model->portal_list->pdata[i];
        portal->id = i;
        }
    for(i = 0; i < model->span_list->len; i++){
        span = model->span_list->pdata[i];
        span->id = i;
        }
    return;
    }

static gboolean C4_Model_is_valid(C4_Model *model){
    register gint i;
    register C4_State *state;
    g_assert(model);
    for(i = 0; i < model->state_list->len; i++){
        state = model->state_list->pdata[i];
        /* All states except start must have input transitions */
        if(state == model->start_state->state)
            g_assert(!state->input_transition_list->len);
        else
            g_assert(state->input_transition_list->len);
        /* All states except end must have output transitions */
        if(state == model->end_state->state)
            g_assert(!state->output_transition_list->len);
        else
            g_assert(state->output_transition_list->len);
        }
    return TRUE;
    }

static void c4_g_ptr_array_reverse(GPtrArray *ptr_array){
    register gint a, z;
    register gpointer swap;
    for(a = 0, z = ptr_array->len-1; a < z; a++, z--){
        swap = ptr_array->pdata[a];
        ptr_array->pdata[a] = ptr_array->pdata[z];
        ptr_array->pdata[z] = swap;
        }
    return;
    }

/* Calculate dependency ordering and transition precedence:
 */
static void C4_Model_topological_sort(C4_Model *model){
    register gint *dependent
             = g_new0(gint, model->transition_list->len);
    register GPtrArray *ordered = g_ptr_array_new();
    register gint i, j;
    register C4_Transition *transition, *input_transition;
    register C4_State *state;
    register gboolean removed_transition;
    /**/
    /* Set initial dependencies */
    for(i = 0; i < model->transition_list->len; i++){
        transition = model->transition_list->pdata[i];
        if((transition->advance_query == 0)
        && (transition->advance_target == 0)){
            state = transition->input;
            for(j = 0; j < state->input_transition_list->len; j++){
                input_transition
                    = state->input_transition_list->pdata[j];
                if((input_transition->advance_query == 0)
                && (input_transition->advance_target == 0)){
                    dependent[input_transition->id]++;
                    }
                }
            }
        }
    /* Perform topological sort of static transitions */
    do {
        removed_transition = FALSE;
        for(i = 0; i < model->transition_list->len; i++){
            if(dependent[i] != 0)
                continue;
            transition = model->transition_list->pdata[i];
            if((transition->advance_query != 0)
            || (transition->advance_target != 0))
                continue;
            removed_transition = TRUE;
            dependent[transition->id] = -1;
            g_ptr_array_add(ordered, transition);
            state = transition->input;
            for(j = 0; j < state->input_transition_list->len; j++){
                input_transition
                    = state->input_transition_list->pdata[j];
                if((transition->advance_query == 0)
                && (transition->advance_target == 0))
                    dependent[input_transition->id]--;
                }
            }
    } while(removed_transition);
    /* Take the precedence zero transitions */
    for(i = 0; i < model->transition_list->len; i++){
        transition = model->transition_list->pdata[i];
        if(transition->advance_query || transition->advance_target)
            g_ptr_array_add(ordered, transition);
        }
    /**/
    c4_g_ptr_array_reverse(ordered);
    /* Assert that there are no cycles */
    g_assert(ordered->len == model->transition_list->len);
    /* Reorder the transition list */
    for(i = 0; i < ordered->len; i++){
        transition = ordered->pdata[i];
        transition->id = i;
        model->transition_list->pdata[i] = transition;
        }
    /**/
    g_ptr_array_free(ordered, TRUE);
    g_free(dependent);
    return;
    }
/* FIXME: optimisation:
 *        transitions are only ordered;
 *        may wish to store precedence info instead
 *        to allow some optimisations.
 */

static void C4_Portal_find_transitions(C4_Portal *portal,
                                       C4_Model *model){
    register gint i;
    register C4_Transition *transition;
    /* Find applicable transitions for portal->calc */
    for(i = 0; i < model->transition_list->len; i++){
        transition = model->transition_list->pdata[i];
        if((transition->calc == portal->calc)
        && (transition->input == transition->output)){
            g_assert(transition->advance_query
                  == portal->advance_query);
            g_assert(transition->advance_target
                  == portal->advance_target);
            g_ptr_array_add(portal->transition_list, transition);
            }
        }
    g_assert(portal->transition_list->len);
    return;
    }

static void C4_Model_finalise(C4_Model *model){
    register gint i;
    register C4_Portal *portal;
    register C4_Transition *transition;
    g_assert(model);
    for(i = 0; i < model->portal_list->len; i++){
        portal = model->portal_list->pdata[i];
        if(portal->transition_list->len)
            g_ptr_array_set_size(portal->transition_list, 0);
        C4_Portal_find_transitions(portal, model);
        }
    model->max_query_advance = 0;
    model->max_target_advance = 0;
    for(i = 0; i < model->transition_list->len; i++){
        transition = model->transition_list->pdata[i];
        if(model->max_query_advance < transition->advance_query)
            model->max_query_advance = transition->advance_query;
        if(model->max_target_advance < transition->advance_target)
            model->max_target_advance = transition->advance_target;
        }
    g_assert(model->max_query_advance || model->max_target_advance);
    return;
    }

/**/

static void C4_Shadow_designate_recur(C4_Shadow *shadow,
                                      C4_Transition *transition,
                                      gboolean *des,
                                      gboolean *states_visited){
    register gint i;
    register C4_State *state = transition->input;
    if(!states_visited[state->id]){
        states_visited[state->id] = TRUE;
        /* If transition is an dst transition stop */
        for(i = 0; i < transition->dst_shadow_list->len; i++){
            if(shadow == transition->dst_shadow_list->pdata[i])
                return;
            }
        /* Visit input transitions of input state */
        for(i = 0; i < state->input_transition_list->len; i++){
            transition = state->input_transition_list->pdata[i];
            g_assert(!des[transition->id]);
            des[transition->id] = TRUE;
            C4_Shadow_designate_recur(shadow, transition,
                                      des, states_visited);
            }
        }
    return;
    }

static gboolean *C4_Shadow_get_designation(C4_Model *model,
                                           C4_Shadow *shadow){
    register gboolean *des = g_new0(gboolean,
                                    model->transition_list->len);
    register gboolean *states_visited = g_new0(gboolean,
                                               model->state_list->len);
    register gint i;
    register C4_Transition *transition;
    g_assert(C4_Shadow_is_valid(shadow, model));
    for(i = 0; i < shadow->dst_transition_list->len; i++){
        transition = shadow->dst_transition_list->pdata[i];
        des[transition->id] = TRUE;
        C4_Shadow_designate_recur(shadow, transition,
                                  des, states_visited);
        }
    g_free(states_visited);
    return des;
    }

static gboolean C4_Shadow_designation_fits(C4_Model *model,
                          gboolean *des_a, gboolean *des_b){
    register gint i;
    register gboolean *state_is_used;
    register C4_Transition *transition;
    for(i = 0; i < model->transition_list->len; i++)
        if(des_a[i] && des_b[i])
            return FALSE;
    state_is_used = g_new0(gboolean, model->state_list->len);
    /* Fail if any des_a output states are des_b inputs */
    for(i = 0; i < model->transition_list->len; i++)
        if(des_a[i]){
            transition = model->transition_list->pdata[i];
            state_is_used[transition->output->id] = TRUE;
            }
    for(i = 0; i < model->transition_list->len; i++)
        if(des_b[i]){
            transition = model->transition_list->pdata[i];
            if(state_is_used[transition->input->id]){
                g_free(state_is_used);
                return FALSE;
                }
            }
    for(i = 0; i < model->state_list->len; i++)
        state_is_used[i] = FALSE;
    /* Fail if any des_b output states are des_a inputs */
    for(i = 0; i < model->transition_list->len; i++)
        if(des_b[i]){
            transition = model->transition_list->pdata[i];
            state_is_used[transition->output->id] = TRUE;
            }
    for(i = 0; i < model->transition_list->len; i++)
        if(des_a[i]){
            transition = model->transition_list->pdata[i];
            if(state_is_used[transition->input->id]){
                g_free(state_is_used);
                return FALSE;
                }
            }
    g_free(state_is_used);
    return TRUE;
    }

static void C4_Shadow_designation_join(C4_Model *model,
                          gboolean *master, gboolean *copy){
    register gint i;
    for(i = 0; i < model->transition_list->len; i++){
        if(copy[i]){
            g_assert(!master[i]);
            master[i] = TRUE;
            }
        }
    return;
    }

static void C4_Model_designate_shadows(C4_Model *model){
    register gboolean *des, *curr_des;
    register GPtrArray *designation_list = g_ptr_array_new();
    register gint i, j;
    register C4_Shadow *shadow;
    for(i = 0; i < model->shadow_list->len; i++){
        shadow = model->shadow_list->pdata[i];
        curr_des = C4_Shadow_get_designation(model, shadow);
        shadow->designation = -1;
        for(j = 0; j < designation_list->len; j++){
            des = designation_list->pdata[j];
            if(C4_Shadow_designation_fits(model, des, curr_des)){
                C4_Shadow_designation_join(model, des, curr_des);
                shadow->designation = j;
                break;
                }
            }
        if(shadow->designation == -1){ /* not designated */
            shadow->designation = designation_list->len;
            g_ptr_array_add(designation_list, curr_des);
        } else {
            g_free(curr_des);
            }
        }
    model->total_shadow_designations = designation_list->len;
    for(i = 0; i < designation_list->len; i++)
        g_free(designation_list->pdata[i]);
    g_ptr_array_free(designation_list, TRUE);
    return;
    }

void C4_Model_close(C4_Model *model){
    g_assert(model);
    g_assert(model->is_open);
    g_assert(C4_Model_is_valid(model));
    C4_Model_set_ids(model);
    g_assert(C4_Model_has_start_to_end_path(model));
    C4_Model_topological_sort(model);
    C4_Model_designate_shadows(model);
    C4_Model_finalise(model);
    model->is_open = FALSE;
    return;
    }

/**/

void C4_Calc_init(C4_Calc *calc, Region *region, gpointer user_data){
    if(!calc)
        return;
    if(calc->init_func)
        calc->init_func(region, user_data);
    return;
    }

void C4_Calc_exit(C4_Calc *calc, Region *region, gpointer user_data){
    if(!calc)
        return;
    if(calc->exit_func)
        calc->exit_func(region, user_data);
    return;
    }

C4_Score C4_Calc_score(C4_Calc *calc, gint query_pos, gint target_pos,
                       gpointer user_data){
    register C4_Score score;
    if(!calc)
        return 0;
    if(calc->calc_func){
        score = calc->calc_func(query_pos, target_pos, user_data);
        g_assert(score <= calc->max_score);
        return score;
        }
    return calc->max_score;
    }

C4_Model *C4_Model_copy(C4_Model *old_model){
    register C4_Model *new_model;
    register gchar *name;
    register gint i, j;
    register C4_Calc *calc, **calc_map = g_new0(C4_Calc*,
                                        old_model->calc_list->len);
    register C4_State *state, **state_map = g_new0(C4_State*,
                                        old_model->state_list->len);
    register C4_State *src_state, *dst_state;
    register C4_Transition *transition,
                           **transition_map = g_new0(C4_Transition*,
                                     old_model->transition_list->len);
    register C4_Shadow *shadow, *new_shadow;
    register C4_Portal *portal;
    register C4_Span *span;
    g_assert(old_model);
    g_assert(!old_model->is_open);
    name = g_strdup_printf("Copy:%s", old_model->name);
    new_model = C4_Model_create(name);
    g_free(name);
    /* Copy calcs */
    for(i = 0; i < old_model->calc_list->len; i++){
        calc = old_model->calc_list->pdata[i];
        calc_map[i] = C4_Model_add_calc(new_model, calc->name,
                           calc->max_score,
                           calc->calc_func, calc->calc_macro,
                           calc->init_func, calc->init_macro,
                           calc->exit_func, calc->exit_macro,
                           calc->protect);
        }
    /* Copy states */
    for(i = 0; i < old_model->state_list->len; i++){
        state = old_model->state_list->pdata[i];
        if((state == old_model->start_state->state)
        || (state == old_model->end_state->state))
            continue;
        state_map[i] = C4_Model_add_state(new_model, state->name);
        }
    /* Copy transitions */
    for(i = 0; i < old_model->transition_list->len; i++){
        transition = old_model->transition_list->pdata[i];
        if(transition->input == old_model->start_state->state)
            src_state = NULL;
        else
            src_state = state_map[transition->input->id];
        if(transition->output == old_model->end_state->state)
            dst_state = NULL;
        else
            dst_state = state_map[transition->output->id];
        transition_map[i]
            = C4_Model_add_transition(new_model, transition->name,
                src_state, dst_state,
                transition->advance_query, transition->advance_target,
                transition->calc?calc_map[transition->calc->id]:NULL,
                transition->label, transition->label_data);
        }
    /* Copy shadows */
    for(i = 0; i < old_model->shadow_list->len; i++){
        shadow = old_model->shadow_list->pdata[i];
        state = shadow->src_state_list->pdata[0];
        transition = shadow->dst_transition_list->pdata[0];
        new_shadow = C4_Model_add_shadow(new_model, shadow->name,
              state_map[state->id], transition_map[transition->id],
              shadow->start_func, shadow->start_macro,
              shadow->end_func, shadow->end_macro);
        for(j = 1; j < shadow->src_state_list->len; j++){
            state = shadow->src_state_list->pdata[j];
            C4_Shadow_add_src_state(new_shadow, state_map[state->id]);
            }
        for(j = 1; j < shadow->dst_transition_list->len; j++){
            transition = shadow->dst_transition_list->pdata[j];
            C4_Shadow_add_dst_transition(new_shadow,
                                      transition_map[transition->id]);
            }
        }
    /* Copy portals */
    for(i = 0; i < old_model->portal_list->len; i++){
        portal = old_model->portal_list->pdata[i];
        C4_Model_add_portal(new_model, portal->name,
            calc_map[portal->calc->id],
            portal->advance_query, portal->advance_target);
        }
    /* Copy spans */
    for(i = 0; i < old_model->span_list->len; i++){
        span = old_model->span_list->pdata[i];
        C4_Model_add_span(new_model, span->name,
                          state_map[span->span_state->id],
                          span->min_query, span->max_query,
                          span->min_target, span->max_target);
        }
    /* Configure extras */
    C4_Model_merge_codegen(new_model, old_model);
    C4_Model_configure_extra(new_model,
             old_model->init_func, old_model->init_macro,
             old_model->exit_func, old_model->exit_macro);
    /* Configure start and end states */
    C4_Model_configure_start_state(new_model,
                     old_model->start_state->scope,
                     old_model->start_state->cell_start_func,
                     old_model->start_state->cell_start_macro);
    C4_Model_configure_end_state(new_model,
                     old_model->end_state->scope,
                     old_model->end_state->cell_end_func,
                     old_model->end_state->cell_end_macro);
    g_free(calc_map);
    g_free(state_map);
    g_free(transition_map);
    /* Close the model without a topological sort */
    C4_Model_set_ids(new_model);
    C4_Model_designate_shadows(new_model);
    C4_Model_finalise(new_model);
    new_model->is_open = FALSE;
    return new_model;
    }

/**/

void C4_Model_remove_transition(C4_Model *model,
                                C4_Transition *transition){
    register gboolean removed_transition;
    /* Remove from the input state */
    removed_transition = g_ptr_array_remove_fast(
            transition->input->output_transition_list, transition);
    g_assert(removed_transition);
    /* Remove from the output state */
    removed_transition = g_ptr_array_remove_fast(
            transition->output->input_transition_list, transition);
    g_assert(removed_transition);
    /* Remove from the Model */
    removed_transition = g_ptr_array_remove_fast(model->transition_list,
                                      transition);
    g_assert(removed_transition);
    C4_Transition_destroy(transition);
    return;
    }

static void C4_Model_remove_shadow(C4_Model *model,
                                   C4_Shadow *shadow){
    register gint i;
    register gboolean removed_shadow;
    register C4_State *state;
    register C4_Transition *transition;
    /* Remove from src states */
    for(i = 0; i < shadow->src_state_list->len; i++){
        state = shadow->src_state_list->pdata[i];
        removed_shadow = g_ptr_array_remove_fast(
            state->src_shadow_list, shadow);
        g_assert(removed_shadow);
        }
    /* Remove from dst transitions */
    for(i = 0; i < shadow->dst_transition_list->len; i++){
        transition = shadow->dst_transition_list->pdata[i];
        removed_shadow = g_ptr_array_remove_fast(
            transition->dst_shadow_list, shadow);
        g_assert(removed_shadow);
        }
    /* Remove from the Model */
    removed_shadow = g_ptr_array_remove_fast(model->shadow_list,
                                             shadow);
    g_assert(removed_shadow);
    C4_Shadow_destroy(shadow);
    return;
    }

void C4_Model_remove_state(C4_Model *model, C4_State *state){
    register gboolean was_open, removed_state;
    register gint i;
    register C4_Transition *transition;
    register C4_Shadow *shadow;
    register C4_Calc *calc;
    register C4_Portal *portal;
    register C4_Span *span;
    register gboolean *calc_used;
    /* Check not start or end */
    g_assert(model);
    g_assert(state);
    g_assert(state != model->start_state->state);
    g_assert(state != model->end_state->state);
    /* Open model if not open */
    was_open = model->is_open;
    if(!was_open)
        C4_Model_open(model);
    /* Remove from state */
    removed_state = g_ptr_array_remove_fast(model->state_list, state);
    g_assert(removed_state);
    /* Remove unused transitions */
    for(i = 0; i < state->input_transition_list->len; i++){
        transition = state->input_transition_list->pdata[i];
        C4_Model_remove_transition(model, transition);
        }
    for(i = 0; i < state->output_transition_list->len; i++){
        transition = state->input_transition_list->pdata[i];
        C4_Model_remove_transition(model, transition);
        }
    /* Remove state from shadows */
    for(i = 0; i < state->src_shadow_list->len; i++){
        shadow = state->src_shadow_list->pdata[i];
        if(shadow->src_state_list->len == 1){ /* Remove entire shadow */
            C4_Model_remove_shadow(model, shadow);
        } else { /* Remove state from shadow state list */
            removed_state = g_ptr_array_remove_fast(
                    shadow->src_state_list, state);
            g_assert(removed_state);
            }
        }
    /* Remove unused calcs */
    calc_used = g_new0(gboolean, model->calc_list->len);
    for(i = 0; i < model->calc_list->len; i++){
        calc = model->calc_list->pdata[i];
        calc->id = i;
        }
    for(i = 0; i < model->transition_list->len; i++){
        transition = model->transition_list->pdata[i];
        calc_used[transition->calc->id] = TRUE;
        }
    for(i = model->calc_list->len-1; i >= 0; i--){
        if(!calc_used[i]){
            calc = model->calc_list->pdata[i];
            g_ptr_array_remove_fast(model->calc_list, calc);
            C4_Calc_destroy(calc);
            }
        }
    /* Remove unused portals */
    for(i = model->portal_list->len; i >= 0; i--){
        portal = model->portal_list->pdata[i];
        if(!calc_used[portal->calc->id]){
            g_ptr_array_remove_fast(model->portal_list, portal);
            C4_Portal_destroy(portal);
            }
        }
    g_free(calc_used);
    /* Remove related spans */
    for(i = model->span_list->len-1; i >= 0; i--){
        span = model->span_list->pdata[i];
        if(span->span_state == state){
            g_ptr_array_remove_fast(model->span_list, span);
            C4_Span_destroy(span);
            }
        }
    /* Close model */
    if(!was_open)
        C4_Model_close(model);
    C4_State_destroy(state);
    return;
    }

gboolean C4_Model_is_global(C4_Model *model){
    if(model->start_state->scope != C4_Scope_CORNER)
        return FALSE;
    if(model->end_state->scope != C4_Scope_CORNER)
        return FALSE;
    return TRUE;
    }

gboolean C4_Model_is_local(C4_Model *model){
    if(model->start_state->scope != C4_Scope_ANYWHERE)
        return FALSE;
    if(model->end_state->scope != C4_Scope_ANYWHERE)
        return FALSE;
    return TRUE;
    }

void C4_Model_remove_all_shadows(C4_Model *model){
    register C4_Shadow *shadow;
    g_assert(model->is_open);
    while(model->shadow_list->len){
        shadow = model->shadow_list->pdata[0];
        g_assert(shadow);
        C4_Model_remove_shadow(model, shadow);
        }
    return;
    }

/**/

typedef struct {
    GPtrArray *new_src_state_list;
    GPtrArray *new_dst_transition_list;
} C4_ProtoShadow;

static C4_ProtoShadow *C4_ProtoShadow_create(void){
    register C4_ProtoShadow *cps = g_new(C4_ProtoShadow, 1);
    cps->new_src_state_list = g_ptr_array_new();
    cps->new_dst_transition_list = g_ptr_array_new();
    return cps;
    }

static void C4_ProtoShadow_destroy(C4_ProtoShadow *cps){
    g_ptr_array_free(cps->new_src_state_list, TRUE);
    g_ptr_array_free(cps->new_dst_transition_list, TRUE);
    g_free(cps);
    return;
    }

static void C4_ProtoShadow_add_state(C4_ProtoShadow *cps,
                                     C4_State *state){
    g_ptr_array_add(cps->new_src_state_list, state);
    return;
    }

static void C4_ProtoShadow_add_transition(C4_ProtoShadow *cps,
                                          C4_Transition *transition){
    g_ptr_array_add(cps->new_dst_transition_list, transition);
    return;
    }

static void C4_ProtoShadow_generate(C4_ProtoShadow *cps,
                                    C4_Model *old_model,
                                    C4_Model *new_model,
                                    C4_Shadow *old_shadow){
    register gint i;
    register C4_Shadow *new_shadow;
    g_assert(cps);
    g_assert(cps->new_src_state_list->len);
    g_assert(cps->new_dst_transition_list->len);
    new_shadow = C4_Model_add_shadow(new_model,
                old_shadow->name,
                cps->new_src_state_list->pdata[0],
                cps->new_dst_transition_list->pdata[0],
                old_shadow->start_func, old_shadow->start_macro,
                old_shadow->end_func, old_shadow->end_macro);
    for(i = 1; i < cps->new_src_state_list->len; i++)
        C4_Shadow_add_src_state(new_shadow,
                cps->new_src_state_list->pdata[i]);
    for(i = 1; i < cps->new_dst_transition_list->len; i++)
        C4_Shadow_add_dst_transition(new_shadow,
                cps->new_dst_transition_list->pdata[i]);
    return;
    }

/**/

static void C4_Model_segment_reuse_state(C4_Model *old_model,
                                         C4_Model *new_model,
                             C4_State *old_state,
                             C4_State **state_map,
                             C4_ProtoShadow **proto_shadow_map){
    register C4_State *new_state;
    register gint i;
    register C4_Shadow *shadow;
    if(old_state == old_model->start_state->state)
        return;
    if(old_state == old_model->end_state->state)
        return;
    if(state_map[old_state->id]) /* Already reused */
        return;
    new_state = C4_Model_add_state(new_model, old_state->name);
    state_map[old_state->id] = new_state;
    for(i = 0; i < old_state->src_shadow_list->len; i++){
        shadow = old_state->src_shadow_list->pdata[i];
        if(!proto_shadow_map[shadow->id])
            proto_shadow_map[shadow->id] = C4_ProtoShadow_create();
        C4_ProtoShadow_add_state(proto_shadow_map[shadow->id],
                                 new_state);
        }
    return;
    }

static void C4_Model_segment_add_transition(C4_Model *old_model,
                 C4_Model *new_model,
                 C4_Transition *transition, C4_State **state_map,
                 C4_Calc **calc_map, GTree *transition_map_tree,
                 C4_ProtoShadow **proto_shadow_map,
                 gboolean from_start, gboolean to_end){
    register C4_Calc *calc = NULL;
    register C4_Transition *new_transition;
    register gint i;
    register C4_Shadow *shadow;
    if(!from_start)
        C4_Model_segment_reuse_state(old_model, new_model,
                                     transition->input, state_map,
                                     proto_shadow_map);
    if(!to_end)
        C4_Model_segment_reuse_state(old_model, new_model,
                                     transition->output, state_map,
                                     proto_shadow_map);
    if(transition->calc){
        if(!calc_map[transition->calc->id]){
            calc_map[transition->calc->id]
                = C4_Model_add_calc(new_model,
                     transition->calc->name,
                     transition->calc->max_score,
                     transition->calc->calc_func,
                     transition->calc->calc_macro,
                     transition->calc->init_func,
                     transition->calc->init_macro,
                     transition->calc->exit_func,
                     transition->calc->exit_macro,
                     transition->calc->protect);
            }
        calc = calc_map[transition->calc->id];
        }
    new_transition = C4_Model_add_transition(new_model,
          transition->name,
          from_start?NULL:state_map[transition->input->id],
          to_end?NULL:state_map[transition->output->id],
          transition->advance_query, transition->advance_target,
          calc, transition->label, transition->label_data);
    g_tree_insert(transition_map_tree, new_transition, transition);
    for(i = 0; i < transition->dst_shadow_list->len; i++){
        shadow = transition->dst_shadow_list->pdata[i];
        if(!proto_shadow_map[shadow->id])
            proto_shadow_map[shadow->id] = C4_ProtoShadow_create();
        C4_ProtoShadow_add_transition(proto_shadow_map[shadow->id],
                                      new_transition);
        }
    return;
    }

static void C4_Model_segment_recur(C4_Model *old_model,
                                   C4_Model *new_model,
                                   C4_State *state,
                                   gboolean *visited_state,
                                   C4_State **state_map,
                                   C4_Calc  **calc_map,
                                   GTree *transition_map_tree,
                                   C4_ProtoShadow **proto_shadow_map){
    register gint i;
    register C4_Transition *transition;
    if(!state_map[state->id]) /* Check that it is on the state map */
        return;
    if(visited_state[state->id]) /* Check is it not already visited */
        return;
    if(state == old_model->start_state->state)
        return;
    if(state == old_model->end_state->state)
        return;
    visited_state[state->id] = TRUE;
    for(i = 0; i < state->output_transition_list->len; i++){
        transition = state->output_transition_list->pdata[i];
        if(transition->output == old_model->end_state->state)
            continue;
        C4_Model_segment_add_transition(old_model, new_model,
           transition, state_map, calc_map,
           transition_map_tree, proto_shadow_map, FALSE, FALSE);
        C4_Model_segment_recur(old_model, new_model,
                               transition->output,
                               visited_state, state_map, calc_map,
                               transition_map_tree, proto_shadow_map);
        }
    return;
    }

#if 0
/* FIXME: not currently used */

static void C4_State_remove_duplicate_transitions(C4_State *state,
                               C4_Model *model, GPtrArray *to_remove){
    register gint i, j;
    register C4_Score score_a, score_b;
    register C4_Transition *transition_a, *transition_b;
    for(i = 0; i < state->output_transition_list->len; i++){
        transition_a = state->output_transition_list->pdata[i];
        for(j = 0; j < i; j++){
            transition_b = state->output_transition_list->pdata[j];
            if(transition_a->output != transition_b->output)
                continue;
            if(transition_a->advance_query
            != transition_b->advance_query)
                continue;
            if(transition_a->advance_target
            != transition_b->advance_target)
                continue;
            if(transition_a->calc == transition_b->calc)
                C4_Model_remove_transition(model, transition_b);
            if(transition_a->calc && transition_a->calc->calc_func)
                continue;
            if(transition_b->calc && transition_b->calc->calc_func)
                continue;
            score_a = transition_a->calc ? transition_a->calc->max_score
                                         : 0;
            score_b = transition_b->calc ? transition_b->calc->max_score
                                         : 0;
            if(score_a > score_b){
                g_ptr_array_add(to_remove, transition_b);
            } else {
                g_ptr_array_add(to_remove, transition_a);
                }
            }
        }
    return;
    }

static void C4_Model_remove_duplicate_transitions(C4_Model *model){
    register gint i;
    register C4_State *state;
    register C4_Transition *transition;
    register GPtrArray *to_remove = g_ptr_array_new();
    /* Find transitions to remove */
    for(i = 0; i < model->state_list->len; i++){
        state = model->state_list->pdata[i];
        C4_State_remove_duplicate_transitions(state, model, to_remove);
        }
    /* Remove transitions from the model */
    for(i = 0; i < to_remove->len; i++){
        transition = to_remove->pdata[i];
        C4_Model_remove_transition(model, transition);
        }
    g_ptr_array_free(to_remove, TRUE);
    return;
    }

#endif /* 0 */

static C4_Model *C4_Model_select(C4_Model *model,
                                 C4_State *src, C4_State *dst,
                                 C4_State **state_map,
                                 GTree *transition_map_tree){
    register gchar *name = g_strdup_printf(
                                   "Segment(\"%s\"->\"%s\"):[%s]",
                                   src->name, dst->name, model->name);
    register C4_Model *segment_model = C4_Model_create(name);
    register C4_Transition *transition;
    register gboolean *visited_state = g_new0(gboolean,
                                       model->state_list->len);
    register C4_State *state;
    register C4_Calc **calc_map = g_new0(C4_Calc*,
                                         model->calc_list->len);
    register gint i;
    register C4_Shadow *shadow;
    register C4_ProtoShadow **proto_shadow_map = g_new0(C4_ProtoShadow*,
                                  model->shadow_list->len);
    /* Propagate model extras */
    C4_Model_configure_extra(segment_model,
             model->init_func, model->init_macro,
             model->exit_func, model->exit_macro);
    C4_Model_merge_codegen(segment_model, model);
    /* Propagate any shadows from start */
    for(i = 0; i < src->src_shadow_list->len; i++){
        shadow = src->src_shadow_list->pdata[i];
        if(!proto_shadow_map[shadow->id])
            proto_shadow_map[shadow->id] = C4_ProtoShadow_create();
        C4_ProtoShadow_add_state(proto_shadow_map[shadow->id],
                                 segment_model->start_state->state);
        }
    /* Transitions from start */
    for(i = 0; i < src->output_transition_list->len; i++){
        transition = src->output_transition_list->pdata[i];
        if(!C4_Model_path_is_possible(model, transition->output, dst))
            continue;
        C4_Model_segment_add_transition(model, segment_model,
           transition, state_map, calc_map,
           transition_map_tree, proto_shadow_map, TRUE, FALSE);
        }
    /* Transitions to end */
    for(i = 0; i < dst->input_transition_list->len; i++){
        transition = dst->input_transition_list->pdata[i];
        if(!C4_Model_path_is_possible(model, src, transition->input))
            continue;
        C4_Model_segment_add_transition(model, segment_model,
           transition, state_map, calc_map,
           transition_map_tree, proto_shadow_map, FALSE, TRUE);
        }
    /* Other transitions */
    for(i = 0; i < model->state_list->len; i++){
        state = model->state_list->pdata[i];
        C4_Model_segment_recur(model, segment_model,
                               state, visited_state,
                               state_map, calc_map,
                               transition_map_tree,
                               proto_shadow_map);
        }
    /* Remove duplicate equivalent transitions */
    /* C4_Model_remove_duplicate_transitions(segment_model); */
    for(i = 0; i < model->shadow_list->len; i++){
        if(proto_shadow_map[i]){
            C4_ProtoShadow_generate(proto_shadow_map[i], model,
                                    segment_model,
                                    model->shadow_list->pdata[i]);
            C4_ProtoShadow_destroy(proto_shadow_map[i]);
            }
        }
    g_free(proto_shadow_map);
    g_free(calc_map);
    g_free(visited_state);
    g_free(name);
    return segment_model;
    }

static gint C4_Transition_compare_by_address(gconstpointer a,
                                             gconstpointer b){
    return GPOINTER_TO_INT(a)-GPOINTER_TO_INT(b);
    }

C4_DerivedModel *C4_DerivedModel_create(C4_Model *original_model,
                C4_State *src, C4_State *dst,
                C4_Scope start_scope, C4_CellStartFunc cell_start_func,
                gchar *cell_start_macro,
                C4_Scope end_scope, C4_CellEndFunc cell_end_func,
                gchar *cell_end_macro){
    register C4_DerivedModel *derived_model = g_new(C4_DerivedModel, 1);
    register GTree *transition_map_tree
          = g_tree_new(C4_Transition_compare_by_address);
    register gint i;
    register C4_Transition *original_transition, *derived_transition;
    register C4_State **state_map = g_new0(C4_State*,
                                    original_model->state_list->len);
    g_assert(!original_model->is_open);
    g_assert(src);
    g_assert(dst);
    derived_model->original = C4_Model_share(original_model);
    derived_model->derived = C4_Model_select(original_model, src, dst,
                                             state_map,
                                             transition_map_tree);
    g_free(state_map);
    C4_Model_close(derived_model->derived);
    C4_Model_configure_start_state(derived_model->derived,
            start_scope, cell_start_func, cell_start_macro);
    C4_Model_configure_end_state(derived_model->derived,
            end_scope, cell_end_func, cell_end_macro);
    derived_model->transition_map = g_new(C4_Transition*,
                   derived_model->derived->transition_list->len);
    for(i = 0; i < derived_model->derived->transition_list->len; i++){
        derived_transition
            = derived_model->derived->transition_list->pdata[i];
        g_assert(derived_transition);
        original_transition = g_tree_lookup(transition_map_tree,
                                            derived_transition);
        g_assert(original_transition);
        derived_model->transition_map[derived_transition->id]
                                    = original_transition;
        }
    g_tree_destroy(transition_map_tree);
    return derived_model;
    }
/* FIXME: optimisation: when start and end are the same as original,
 *                      just do a C4_Model_copy(), (with transition map)
 */

void C4_DerivedModel_destroy(C4_DerivedModel *derived_model){
    C4_Model_destroy(derived_model->original);
    C4_Model_destroy(derived_model->derived);
    g_free(derived_model->transition_map);
    g_free(derived_model);
    return;
    }

/**/

