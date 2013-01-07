/****************************************************************\
*                                                                *
*  C4 dynamic programming library - DP layout code               *
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

#include "layout.h"
#include "matrix.h"

/**/

static gboolean Layout_model_has_state_active(C4_Model *model,
                                 C4_State *state,
                                 gint query_pos, gint target_pos,
                                 gint query_length, gint target_length){
    if(query_pos < 0)
        return FALSE;
    if(target_pos < 0)
        return FALSE;
    if(query_pos > query_length)
        return FALSE;
    if(target_pos > target_length)
        return FALSE;
    /**/
    if(state == model->start_state->state){
        switch(model->start_state->scope){
            case C4_Scope_ANYWHERE:
                break;
            case C4_Scope_EDGE:
                if((query_pos != 0) && (target_pos != 0))
                    return FALSE;
                break;
            case C4_Scope_QUERY:
                if(query_pos != 0)
                    return FALSE;
                break;
            case C4_Scope_TARGET:
                if(target_pos != 0)
                    return FALSE;
                break;
            case C4_Scope_CORNER:
                if((query_pos != 0) || (target_pos != 0))
                    return FALSE;
                break;
            default:
                g_assert(!"Bad start scope type");
                break;
            }
        }
    /**/
    if(state == model->end_state->state){
        switch(model->end_state->scope){
            case C4_Scope_ANYWHERE:
                break;
            case C4_Scope_EDGE:
                if((query_pos != query_length)
                && (target_pos != target_length))
                    return FALSE;
                break;
            case C4_Scope_QUERY:
                if(query_pos != query_length)
                    return FALSE;
                break;
            case C4_Scope_TARGET:
                if(target_pos != target_length)
                    return FALSE;
                break;
            case C4_Scope_CORNER:
                if((query_pos != query_length)
                || (target_pos != target_length))
                    return FALSE;
                break;
            default:
                g_assert(!"Bad end scope type");
                break;
            }
        }
    return TRUE;
    }

static gboolean Layout_model_state_is_valid_input(C4_Model *model,
                             C4_State *state,
                             gint query_pos, gint target_pos,
                             gint query_length, gint target_length){
    register gint i;
    register C4_Transition *transition;
    if(state == model->start_state->state)
        return TRUE;
    for(i = 0; i < state->input_transition_list->len; i++){
        transition = state->input_transition_list->pdata[i];
        if((transition->advance_query == 0)
         && (transition->advance_target == 0)){ /* If silent */
            if(Layout_model_state_is_valid_input(model,
                      transition->input,
                      query_pos, target_pos,
                      query_length, target_length))
                return TRUE;
        } else { /* Emits */
            if(Layout_model_has_state_active(model, transition->input,
                   query_pos-transition->advance_query,
                   target_pos-transition->advance_target,
                   query_length, target_length))
                return TRUE;
            }
        }
    return FALSE;
    }
/* The input state must be in scope,
 * and be traceable to a non-silent transition
 * from an state which is also in scope.
 */

static gboolean Layout_transition_is_valid(C4_Model *model,
                           C4_Transition *transition,
                           gint dst_query_pos, gint dst_target_pos,
                           gint query_length, gint target_length){
    /* Check transition input is valid */
    if(!Layout_model_has_state_active(model, transition->input,
        dst_query_pos-transition->advance_query,
        dst_target_pos-transition->advance_target,
        query_length, target_length)){
        return FALSE;
        }
    /* Check transition output is valid */
    if(!Layout_model_has_state_active(model, transition->output,
        dst_query_pos, dst_target_pos,
        query_length, target_length)){
        return FALSE;
        }
    /* Check input score will be present.
     *
     * This must be disabled for continuation DP to work properly.
     * also generated code is smaller without it,
     * and the DP seems to be just as fast.
     */
    if(FALSE){ /* DISABLED */
        if(!Layout_model_state_is_valid_input(model, transition->input,
            dst_query_pos-transition->advance_query,
            dst_target_pos-transition->advance_target,
            query_length, target_length)){
            return FALSE;
            }
        }
    return TRUE;
    }

gboolean Layout_is_transition_valid(Layout *layout, C4_Model *model,
                              C4_Transition *transition,
                              gint dst_query_pos, gint dst_target_pos,
                              gint query_length, gint target_length){
    register Layout_Row *row;
    register Layout_Cell *cell;
    register Layout_Mask *mask;
    g_assert(layout);
    row = layout->row_list->pdata[MIN(dst_target_pos,
                                     (layout->row_list->len-1))];
    g_assert(row);
    cell = row->cell_list->pdata[MIN(dst_query_pos,
                                    (row->cell_list->len-1))];
    g_assert(cell);
    if(dst_query_pos == query_length){
        if(dst_target_pos == target_length){ /* QT */
            mask = cell->corner;
            if(!mask)
                mask = cell->end_target;
            if(!mask)
                mask = cell->normal;
        } else { /* Q */
            mask = cell->end_query;
            if(!mask)
                mask = cell->normal;
            }
    } else {
        if(dst_target_pos == target_length){ /* T */
            mask = cell->end_target;
            if(!mask)
                mask = cell->normal;
        } else { /* - */
            mask = cell->normal;
            }
        }
    g_assert(mask);
    g_assert(Layout_transition_is_valid(model, transition,
                     dst_query_pos, dst_target_pos,
                     query_length, target_length)
          == g_array_index(mask->transition_mask, gboolean,
                           transition->id));
    return g_array_index(mask->transition_mask, gboolean,
                         transition->id);
    }
/*
We could just call Layout_transition_is_valid(),
but this way is slightly faster.

FIXME: should change Layout to allow for faster lookup here.
       (currently makes code generation fast, but interp slow ...)
*/

/**/

static Layout_Mask *Layout_Mask_create(C4_Model *model,
                            gint query_pos, gint target_pos,
                            gint query_length, gint target_length){
    register Layout_Mask *mask = g_new(Layout_Mask, 1);
    register gint i;
    register C4_Transition *transition;
    gboolean is_valid;
    mask->transition_mask = g_array_new(FALSE, TRUE, sizeof(gboolean));
    for(i = 0; i < model->transition_list->len; i++){
        transition = model->transition_list->pdata[i];
        is_valid = Layout_transition_is_valid(model, transition,
                    query_pos, target_pos, query_length, target_length);
        g_array_append_val(mask->transition_mask, is_valid);
        }
    return mask;
    }

static void Layout_Mask_destroy(Layout_Mask *mask){
    g_array_free(mask->transition_mask, TRUE);
    g_free(mask);
    return;
    }

static gboolean Layout_Mask_is_same(Layout_Mask *a, Layout_Mask *b){
    register gint i;
    register gboolean value_a, value_b;
    if(!(a || b))
        return TRUE;
    if(!(a && b))
        return FALSE;
    if(a->transition_mask->len != b->transition_mask->len)
        return FALSE;
    for(i = 0; i < a->transition_mask->len; i++){
        value_a = g_array_index(a->transition_mask, gboolean, i);
        value_b = g_array_index(b->transition_mask, gboolean, i);
        if(value_a != value_b)
            return FALSE;
        }
    return TRUE;
    }

static void Layout_Mask_info(Layout_Mask *mask, C4_Model *model,
                             gchar *name){
    register gint i, count = 0;
    register gboolean is_valid;
    if(!mask)
        return;
    g_assert(mask->transition_mask->len == model->transition_list->len);
    g_print("[");
    for(i = 0; i < mask->transition_mask->len; i++){
        is_valid = g_array_index(mask->transition_mask, gboolean, i);
        if(is_valid){
            count++;
            g_print("#");
        } else {
            g_print(" ");
            }
        }
    g_print("] [%s] [%d]\n", name, count);
    return;
    }

/**/

static Layout_Cell *Layout_Cell_create(C4_Model *model,
                            gint query_pos, gint target_pos,
                            gint query_length, gint target_length){
    register Layout_Cell *cell = g_new(Layout_Cell, 1);
    cell->normal = Layout_Mask_create(model,
                        query_pos, target_pos,
                        query_length, target_length);
    cell->end_query = Layout_Mask_create(model,
                        query_pos, target_pos,
                        query_pos, target_length);
    cell->end_target = Layout_Mask_create(model,
                        query_pos, target_pos,
                        query_length, target_pos);
    cell->corner = Layout_Mask_create(model,
                        query_pos, target_pos,
                        query_pos, target_pos);
    /* Skip EndQuery if same as Normal */
    if(Layout_Mask_is_same(cell->normal, cell->end_query)){
        Layout_Mask_destroy(cell->end_query);
        cell->end_query = NULL;
        }
    /* Skip Corner if same as EndTarget */
    if(Layout_Mask_is_same(cell->corner, cell->end_target)){
        Layout_Mask_destroy(cell->corner);
        cell->corner = NULL;
        /* Also skip EndTarget if same as Normal */
        if(Layout_Mask_is_same(cell->normal, cell->end_target)){
            Layout_Mask_destroy(cell->end_target);
            cell->end_target = NULL;
            }
        }
/* FIXME: disabled as this was incorrect for some non-local models */
#if 0
    /* Skip EndTarget (and Corner if present) if same as Normal
     * and Corner?Corner:Normal is same as EndQuery?EndQuery:Normal
     */
    if(Layout_Mask_is_same(cell->normal, cell->end_target)
    && Layout_Mask_is_same(cell->corner ? cell->corner
                                        : cell->end_target,
                           cell->end_query ? cell->end_query
                                           : cell->normal)){
        Layout_Mask_destroy(cell->end_target);
        cell->end_target = NULL;
        if(cell->corner){
            Layout_Mask_destroy(cell->corner);
            cell->corner = NULL;
            }
        }
#endif /* 0 */
    return cell;
    }

static void Layout_Cell_destroy(Layout_Cell *cell){
    g_assert(cell);
    g_assert(cell->normal);
    Layout_Mask_destroy(cell->normal);
    if(cell->end_query)
        Layout_Mask_destroy(cell->end_query);
    if(cell->end_target)
        Layout_Mask_destroy(cell->end_target);
    if(cell->corner)
        Layout_Mask_destroy(cell->corner);
    g_free(cell);
    return;
    }

static gboolean Layout_Cell_is_same(Layout_Cell *a, Layout_Cell *b){
    if(!Layout_Mask_is_same(a->normal, b->normal))
        return FALSE;
    if(!Layout_Mask_is_same(a->end_query, b->end_query))
        return FALSE;
    if(!Layout_Mask_is_same(a->end_target, b->end_target))
        return FALSE;
    if(!Layout_Mask_is_same(a->corner, b->corner))
        return FALSE;
    return TRUE;
    }

static void Layout_Cell_info(Layout_Cell *cell, C4_Model *model){
    Layout_Mask_info(cell->normal,     model, "normal");
    Layout_Mask_info(cell->end_query,  model, "end_query");
    Layout_Mask_info(cell->end_target, model, "end_target");
    Layout_Mask_info(cell->corner,     model, "corner");
    return;
    }

/**/

static Layout_Row *Layout_Row_create(C4_Model *model,
               gint row_number, gint query_length, gint target_length){
    register Layout_Row *row = g_new(Layout_Row, 1);
    register Layout_Cell *cell;
    row->cell_list = g_ptr_array_new();
    do {
        cell = Layout_Cell_create(model, row->cell_list->len,
                           row_number, query_length, target_length);
        if(row->cell_list->len >= model->max_query_advance)
            if(Layout_Cell_is_same(cell,
                        row->cell_list->pdata[row->cell_list->len-1])){
                Layout_Cell_destroy(cell);
                break;
                }
        g_ptr_array_add(row->cell_list, cell);
        g_assert(row->cell_list->len < 1024); /* safety check */
    } while(TRUE);
    return row;
    }

static void Layout_Row_destroy(Layout_Row *row){
    register gint i;
    for(i = 0; i < row->cell_list->len; i++)
        Layout_Cell_destroy(row->cell_list->pdata[i]);
    g_ptr_array_free(row->cell_list, TRUE);
    g_free(row);
    return;
    }

static gboolean Layout_Row_is_same(Layout_Row *a, Layout_Row *b){
    register gint i;
    register Layout_Cell *cell_a, *cell_b;
    g_assert(a);
    g_assert(b);
    if(a->cell_list->len != b->cell_list->len)
        return FALSE;
    for(i = 0; i < a->cell_list->len; i++){
        cell_a = a->cell_list->pdata[i];
        cell_b = b->cell_list->pdata[i];
        if(!Layout_Cell_is_same(cell_a, cell_b))
            return FALSE;
        }
    return TRUE;
    }

static void Layout_Row_info(Layout_Row *row, C4_Model *model){
    register gint i;
    g_print("Layout_Cell count: [%d]\n", row->cell_list->len);
    for(i = 0; i < row->cell_list->len; i++){
        g_print("Layout_Cell [%d]:\n", i);
        Layout_Cell_info(row->cell_list->pdata[i], model);
        }
    return;
    }

/**/

Layout *Layout_create(C4_Model *model){
    register Layout *layout = g_new(Layout, 1);
    register Layout_Row *row;
    g_assert(model);
    g_assert(!model->is_open);
    layout->row_list = g_ptr_array_new();
    do {
        row = Layout_Row_create(model, layout->row_list->len,
                                1024, 1024);
        if(layout->row_list->len >= model->max_target_advance)
            if(Layout_Row_is_same(row,
                   layout->row_list->pdata[layout->row_list->len-1])){
                Layout_Row_destroy(row);
                break;
                }
        g_ptr_array_add(layout->row_list, row);
        g_assert(layout->row_list->len < 1024); /* safety check */
    } while(TRUE);
    return layout;
    }
/* FIXME: optimisation
 *        may be possible to trim the pattern for some models
 *        where dimensions are <= layout->max_{query,target}_advance
 */

void Layout_destroy(Layout *layout){
    register gint i;
    for(i = 0; i < layout->row_list->len; i++)
        Layout_Row_destroy(layout->row_list->pdata[i]);
    g_ptr_array_free(layout->row_list, TRUE);
    g_free(layout);
    return;
    }

void Layout_info(Layout *layout, C4_Model *model){
    register gint i;
    register C4_Transition *transition;
    g_print("Layout for model [%s]\n", model->name);
    for(i = 0; i < model->transition_list->len; i++){
        transition = model->transition_list->pdata[i];
        g_message("transition [%d] [%s]", i, transition->name);
        }
    g_print("Layout_Row count: [%d]\n", layout->row_list->len);
    for(i = 0; i < layout->row_list->len; i++){
        g_print("Layout_Row [%d]:\n", i);
        Layout_Row_info(layout->row_list->pdata[i], model);
        }
    return;
    }

/**/

