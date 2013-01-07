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

#ifndef INCLUDED_C4_H
#define INCLUDED_C4_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <glib.h>

#include "region.h"
#include "threadref.h"

typedef gint C4_Score;
#define C4_IMPOSSIBLY_LOW_SCORE -987654321
#define C4_IMPOSSIBLY_HIGH_SCORE -(C4_IMPOSSIBLY_LOW_SCORE)

/**/

typedef C4_Score (*C4_CalcFunc)(gint query_pos, gint target_pos,
                                gpointer user_data);
/* Receives the position at the transition source */

typedef void (*C4_PrepFunc)(Region *region, gpointer user_data);
/* To prepare subject region before starting DP */

typedef C4_Score (*C4_StartFunc)(gint query_pos, gint target_pos,
                                 gpointer user_data);
/* Receives the position at the transition source */

typedef void (*C4_EndFunc)(C4_Score score,
                           gint query_pos, gint target_pos,
                           gpointer user_data);
/* Receives the position at the transition end */

typedef C4_Score *(*C4_CellStartFunc)(gint query_pos, gint target_pos,
                                      gpointer user_data);
/* Receives the position at the transition source */

typedef void (*C4_CellEndFunc)(C4_Score *cell, gint cell_size,
                               gint query_pos, gint target_pos,
                               gpointer user_data);
/* Receives the position at the transition end */

/**/

typedef struct {
         gchar *name;
          gint  id;            /* Only set once the model is closed */
     GPtrArray *input_transition_list;
     GPtrArray *output_transition_list;
     GPtrArray *src_shadow_list; /* Shadows starting at this state */
} C4_State;

typedef enum {
    C4_Protect_NONE      = (0<<0),
    C4_Protect_OVERFLOW  = (1<<0),
    C4_Protect_UNDERFLOW = (1<<1)
} C4_Protect;

typedef struct {
          gchar *name;
           gint  id;         /* Only set once the model is closed */
       C4_Score  max_score;
    C4_CalcFunc  calc_func;  /* If NULL, use max_score     */
          gchar *calc_macro; /* If NULL, use calc_func     */
    C4_PrepFunc  init_func;  /* If NULL, no initialisation */
          gchar *init_macro; /* If NULL, use init_func     */
    C4_PrepFunc  exit_func;  /* If NULL, no initialisation */
          gchar *exit_macro; /* If NULL, use exit_func     */
     C4_Protect  protect;
} C4_Calc;
/* When {calc,init,exit}_macro is supplied,
 * {calc,init,exit}_func must be set.
 */

typedef enum {
    C4_Scope_ANYWHERE, /* Anywhere        */
    C4_Scope_EDGE,     /* Query || Target */
    C4_Scope_QUERY,    /* Query only      */
    C4_Scope_TARGET,   /* Target only     */
    C4_Scope_CORNER    /* Query && Target */
} C4_Scope;

typedef struct {
             C4_State *state;
             C4_Scope  scope;
     C4_CellStartFunc  cell_start_func;
                gchar *cell_start_macro;
} C4_StartState;
/* Start cell will be set to zero if cell_start_func is NULL */

typedef struct {
          C4_State *state;
          C4_Scope  scope;
    C4_CellEndFunc  cell_end_func;
             gchar *cell_end_macro;
} C4_EndState;

typedef enum {
    C4_Label_NONE,
    C4_Label_MATCH,
    C4_Label_GAP,
    C4_Label_NER,
    C4_Label_5SS,
    C4_Label_3SS,
    C4_Label_INTRON,
    C4_Label_SPLIT_CODON,
    C4_Label_FRAMESHIFT
} C4_Label;

typedef struct {
            gchar *name;
             gint  id;          /* Only set once the model is closed */
         C4_State *input;
         C4_State *output;
          C4_Calc *calc; /* Zero scores are emitted for a NULL calc */
             gint  advance_query;
             gint  advance_target;
         C4_Label  label;
         gpointer  label_data;
        GPtrArray *dst_shadow_list; /* Shadows ending on transition */
} C4_Transition;

typedef struct {
         gint  id;          /* Only set once the model is closed */
        gchar *name;
    GPtrArray *src_state_list;
    GPtrArray *dst_transition_list;
 C4_StartFunc  start_func;
        gchar *start_macro;
   C4_EndFunc  end_func;
        gchar *end_macro;
         gint  designation;  /* Only set once the model is closed */
} C4_Shadow;

typedef struct {
       gchar *name;
        gint  id;
        gint  advance_query;
        gint  advance_target;
     C4_Calc *calc;
   GPtrArray *transition_list; /* Transitions using this portal */
} C4_Portal;

typedef struct {
            gchar *name;
             gint  id;
         C4_State *span_state;
             gint  min_query;
             gint  max_query;
             gint  min_target;
             gint  max_target;
    C4_Transition *query_loop;  /* May be NULL */
    C4_Transition *target_loop; /* May be NULL */
} C4_Span;

typedef struct {
            gchar  *name;
        GPtrArray  *global_code_list;
        GPtrArray  *local_code_list;
        GPtrArray  *cflags_add_list;
      C4_PrepFunc  init_func;  /* If NULL, no initialisation */
            gchar *init_macro; /* If NULL, use init_func     */
      C4_PrepFunc  exit_func;  /* If NULL, no initialisation */
            gchar *exit_macro; /* If NULL, use exit_func     */
        ThreadRef *thread_ref;
         gboolean   is_open;
        GPtrArray  *state_list;
        GPtrArray  *transition_list;
        GPtrArray  *shadow_list;
        GPtrArray  *calc_list;
        GPtrArray  *portal_list;
        GPtrArray  *span_list;
    C4_StartState  *start_state;
      C4_EndState  *end_state;
             gint   max_query_advance;  /* For convenience */
             gint   max_target_advance; /* For convenience */
             gint   total_shadow_designations;
} C4_Model;

/**/

C4_Model *C4_Model_create(gchar *name);
          /* Model is created open */
    void  C4_Model_destroy(C4_Model *model);
C4_Model *C4_Model_share(C4_Model *model);

C4_State **C4_Model_build_state_map(C4_Model *src, C4_Model *dst);
C4_Transition **C4_Model_build_transition_map(C4_Model *src,
                                              C4_Model *dst);

/**/

void C4_Model_make_stereo(C4_Model *model, gchar *suffix_a,
                                           gchar *suffix_b);
void C4_Model_insert(C4_Model *target, C4_Model *insert,
                     C4_State *src, C4_State *dst);
/* Inserts <insert> into <target> between <src> and <dst> */

/**/

void C4_Model_rename(C4_Model *model, gchar *name);

C4_State *C4_Model_add_state(C4_Model *model, gchar *name);

 C4_Calc *C4_Model_add_calc(C4_Model *model, gchar *name,
                       C4_Score max_score,
                       C4_CalcFunc calc_func, gchar *calc_macro,
                       C4_PrepFunc init_func, gchar *init_macro,
                       C4_PrepFunc exit_func, gchar *exit_macro,
                       C4_Protect protect);
/* calc_macro may contain %[QT][SPRE]
 * which are expanded to {query,target}_{start,position,region_pos,end}
 * in the generated code.
 */

C4_Transition *C4_Model_add_transition(C4_Model *model, gchar *name,
                       C4_State *input, C4_State *output,
                       gint advance_query, gint advance_target,
                       C4_Calc *calc, C4_Label label,
                       gpointer label_data);
/* NULL input implies START state, NULL output implies END state */

GPtrArray *C4_Model_select_transitions(C4_Model *model, C4_Label label);
C4_Transition *C4_Model_select_single_transition(C4_Model *model,
                                                 C4_Label label);

#define C4_Transition_is_match(transition)            \
          (transition->label == C4_Label_MATCH)

#define C4_Transition_is_span(transition)             \
          ((transition->input == transition->output)  \
           &&  (!transition->calc))
/* Transition advance or target must be set, otherwise would be cyclic
 */

C4_Shadow *C4_Model_add_shadow(C4_Model *model, gchar *name,
                    C4_State *src, C4_Transition *dst,
                    C4_StartFunc start_func, gchar *start_macro,
                    C4_EndFunc end_func, gchar *end_macro);
/* NULL src state implies START state,
 * NULL output transition implies all transitions to the END state
 */
void C4_Shadow_add_src_state(C4_Shadow *shadow, C4_State *src);
void C4_Shadow_add_dst_transition(C4_Shadow *shadow,
                                  C4_Transition *dst);

C4_Portal *C4_Model_add_portal(C4_Model *model, gchar *name,
           C4_Calc *calc, gint advance_query, gint advance_target);
/* If terminal_range or join_range are NULL, a default is used.
 */

C4_Span *C4_Model_add_span(C4_Model *model, gchar *name,
                           C4_State *span_state,
                           gint min_query, gint max_query,
                           gint min_target, gint max_target);
/* If range is NULL, a default is used.
 */

void C4_Model_append_codegen(C4_Model *model,
                             gchar *global_code, gchar *local_code,
                             gchar *cflags_add);
void C4_Model_clear_codegen(C4_Model *model);

void C4_Model_configure_extra(C4_Model *model,
                              C4_PrepFunc init_func, gchar *init_macro,
                              C4_PrepFunc exit_func, gchar *exit_macro);
/* init_{func,macro} is called/included at start of DP functions,
 * (after variable declarations, before anything else)
 * exit_{func,macro} is called/included at the end of DP functions.
 *
 * {init,exit}_{func,macro} will be stripped from derived models,
 * but {global,local}_code will not.
 */

void C4_Model_configure_start_state(C4_Model *model, C4_Scope scope,
        C4_CellStartFunc cell_start_func, gchar *cell_start_macro);
void C4_Model_configure_end_state(C4_Model *model, C4_Scope scope,
        C4_CellEndFunc cell_end_func, gchar *cell_end_macro);
/* C4_Model_configure_{start,end}_state()
 *     will work on open or closed models.
 */

void C4_Model_print(C4_Model *model);
void C4_Model_dump_graphviz(C4_Model *model);

void C4_Model_open(C4_Model *model);
void C4_Model_close(C4_Model *model);

/**/

/* Utility functions */

gchar *C4_Scope_get_name(C4_Scope scope);
gchar *C4_Label_get_name(C4_Label label);

gboolean C4_Model_path_is_possible(C4_Model *model, C4_State *src,
                                                    C4_State *dst);

void C4_Calc_init(C4_Calc *calc, Region *region,
                  gpointer user_data);

void C4_Calc_exit(C4_Calc *calc, Region *region,
                  gpointer user_data);

C4_Score C4_Calc_score(C4_Calc *calc, gint query_pos, gint target_pos,
                       gpointer user_data);

C4_Model *C4_Model_copy(C4_Model *model);
/* Returns a closed copy of the model */

void C4_Model_remove_state(C4_Model *model, C4_State *state);
void C4_Model_remove_transition(C4_Model *model,
                                C4_Transition *transition);

gboolean C4_Model_is_global(C4_Model *model);
gboolean C4_Model_is_local(C4_Model *model);
void C4_Model_remove_all_shadows(C4_Model *model);

/**/

typedef struct {
         C4_Model  *original;
         C4_Model  *derived;
    C4_Transition **transition_map;
} C4_DerivedModel;
/* The transition_map maps derived_model transitions
 * back to original_model transitions.
 */

C4_DerivedModel *C4_DerivedModel_create(C4_Model *original_model,
                C4_State *src, C4_State *dst,
                C4_Scope start_scope, C4_CellStartFunc cell_start_func,
                gchar *cell_start_macro,
                C4_Scope end_scope, C4_CellEndFunc cell_end_func,
                gchar *cell_end_macro);
/* A C4_DerivedModel will be created closed
 */

void  C4_DerivedModel_destroy(C4_DerivedModel *derived_model);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* INCLUDED_C4_H */

