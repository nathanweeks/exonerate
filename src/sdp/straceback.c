/****************************************************************\
*                                                                *
*  C4 dynamic programming library - Scheduler Traceback          *
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

#include "straceback.h"

/**/

STraceback *STraceback_create(C4_Model *model, gboolean is_forward){
    register STraceback *straceback = g_new(STraceback, 1);
    straceback->ref_count = 1;
    straceback->is_forward = is_forward;
    straceback->cell_recycle = RecycleBin_create(
            "STraceback_Cell", sizeof(STraceback_Cell), 1024);
    straceback->model = C4_Model_share(model);
    return straceback;
    }

STraceback *STraceback_share(STraceback *straceback){
    straceback->ref_count++;
    return straceback;
    }

void STraceback_destroy(STraceback *straceback){
    if(--straceback->ref_count)
        return;
    RecycleBin_destroy(straceback->cell_recycle);
    C4_Model_destroy(straceback->model);
    g_free(straceback);
    return;
    }

STraceback_Cell *STraceback_add(STraceback *straceback,
                                C4_Transition *transition, gint length,
                                STraceback_Cell *prev){
    register STraceback_Cell *cell
           = RecycleBin_alloc(straceback->cell_recycle);
    g_assert(transition);
    cell->ref_count = 1;
    cell->transition = transition;
    cell->length = length;
    if(prev){
        if(straceback->is_forward){
            g_assert(prev->transition->output == transition->input);
        } else {
            g_assert(prev->transition->input == transition->output);
            }
        cell->prev = STraceback_Cell_share(prev);
    } else {
        if(straceback->is_forward){
             g_assert(transition->input
                   == straceback->model->start_state->state);
        } else {
             g_assert(transition->output
                   == straceback->model->end_state->state);
            }
        cell->prev = NULL;
        }
    return cell;
    }

STraceback_Cell *STraceback_Cell_share(STraceback_Cell *cell){
    cell->ref_count++;
    return cell;
    }

void STraceback_Cell_destroy(STraceback_Cell *cell,
                             STraceback *straceback){
    if(--cell->ref_count)
        return;
    if(cell->prev)
        STraceback_Cell_destroy(cell->prev, straceback);
    RecycleBin_recycle(straceback->cell_recycle, cell);
    return;
    }

static void straceback_g_ptr_array_reverse(GPtrArray *ptr_array){
    register gint a, z;
    register gpointer swap;
    for(a = 0, z = ptr_array->len-1; a < z; a++, z--){
        swap = ptr_array->pdata[a];
        ptr_array->pdata[a] = ptr_array->pdata[z];
        ptr_array->pdata[z] = swap;
        }
    return;
    }

STraceback_List *STraceback_List_create(STraceback *straceback,
                                        STraceback_Cell *prev){
    register STraceback_List *stlist = g_new(STraceback_List, 1);
    register STraceback_Cell *cell;
    register STraceback_Operation *operation;
    g_assert(prev);
    if(straceback->is_forward){
        g_assert(prev->transition->output
              == straceback->model->end_state->state);
    } else {
        g_assert(prev->transition->input
              == straceback->model->start_state->state);
        }
    stlist->operation_list = g_ptr_array_new();
    for(cell = prev; cell; cell = cell->prev){
        operation = g_new(STraceback_Operation, 1);
        operation->transition = cell->transition;
        operation->length = cell->length;
        g_ptr_array_add(stlist->operation_list, operation);
        }
    straceback_g_ptr_array_reverse(stlist->operation_list);
    return stlist;
    }

void STraceback_List_destroy(STraceback_List *stlist){
    register gint i;
    register STraceback_Operation *operation;
    for(i = 0; i < stlist->operation_list->len; i++){
        operation = stlist->operation_list->pdata[i];
        g_free(operation);
        }
    g_ptr_array_free(stlist->operation_list, TRUE);
    g_free(stlist);
    return;
    }

/**/


