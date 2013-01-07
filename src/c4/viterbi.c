/****************************************************************\
*                                                                *
*  C4 dynamic programming library - the Viterbi implementation   *
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

#include <string.h> /* For strlen() */
#include "viterbi.h"
#include "matrix.h"
#include "cgutil.h"

#ifdef USE_COMPILED_MODELS
#include "c4_model_archive.h"
#endif /* USE_COMPILED_MODELS */

/**/

Viterbi_ArgumentSet *Viterbi_ArgumentSet_create(Argument *arg){
    register ArgumentSet *as;
    static Viterbi_ArgumentSet vas = {32};
    if(arg){
        as = ArgumentSet_create("Viterbi algorithm options");
        ArgumentSet_add_option(as, 'D', "dpmemory", "Mb",
           "Maximum memory to use for DP tracebacks (Mb)", "32",
           Argument_parse_int, &vas.traceback_memory_limit);
        Argument_absorb_ArgumentSet(arg, as);
        }
    return &vas;
    }

/**/

static gsize Viterbi_get_cell_size(Viterbi *viterbi){
    register gint cell_size = 1
                            + viterbi->model->total_shadow_designations;
    if(viterbi->mode == Viterbi_Mode_FIND_REGION){
        if(viterbi->model->start_state->scope != C4_Scope_CORNER){
            if(viterbi->model->start_state->scope != C4_Scope_QUERY)
                cell_size++; /* For query_start shadow */
            if(viterbi->model->start_state->scope != C4_Scope_TARGET)
                cell_size++; /* For target_start shadow */
            }
        }
    if(viterbi->mode == Viterbi_Mode_FIND_CHECKPOINTS)
        cell_size++; /*  For checkpoint shadow */
    return cell_size;
    }

Viterbi *Viterbi_create(C4_Model *model, gchar *name,
                        Viterbi_Mode mode, gboolean use_continuation,
                        gboolean use_codegen){
    register Viterbi *viterbi = g_new(Viterbi, 1);
    register Codegen_ArgumentSet *cas
           = Codegen_ArgumentSet_create(NULL);
    g_assert(model);
    g_assert(name);
    viterbi->vas = Viterbi_ArgumentSet_create(NULL);
    viterbi->name = Codegen_clean_path_component(name);
    if(use_continuation){
        viterbi->model = C4_Model_copy(model);
        /* Make model global for continuation */
        C4_Model_configure_start_state(viterbi->model, C4_Scope_CORNER,
             viterbi->model->start_state->cell_start_func,
             viterbi->model->start_state->cell_start_macro);
        C4_Model_configure_end_state(viterbi->model, C4_Scope_CORNER,
             viterbi->model->end_state->cell_end_func,
             viterbi->model->end_state->cell_end_macro);
    } else {
        viterbi->model = C4_Model_share(model);
        }
    viterbi->func = NULL;
    if(use_codegen && cas->use_compiled){
#ifdef USE_COMPILED_MODELS
        viterbi->func = Bootstrapper_lookup(viterbi->name);
        if(!viterbi->func)
            g_warning("Could not find compiled implementation of [%s]",
                      viterbi->name);
#else /* USE_COMPILED_MODELS */
        g_warning("No compiled models - will be very slow");
#endif /* USE_COMPILED_MODELS */
        }
    viterbi->mode = mode;
    viterbi->use_continuation = use_continuation;
    viterbi->cell_size = Viterbi_get_cell_size(viterbi);
    viterbi->layout = Layout_create(viterbi->model);
    return viterbi;
    }

void Viterbi_destroy(Viterbi *viterbi){
    C4_Model_destroy(viterbi->model);
    Layout_destroy(viterbi->layout);
    g_free(viterbi->name);
    g_free(viterbi);
    return;
    }

/**/

static gsize Viterbi_get_row_size(Viterbi *viterbi, Region *region){
    register gint mat_size;
    mat_size = Matrix4d_size(viterbi->model->max_target_advance+1,
                             region->query_length+1,
                             viterbi->model->state_list->len,
                             viterbi->cell_size,
                             sizeof(C4_Score));
    if(!mat_size)
        return 0;
    return sizeof(Viterbi_Row) + mat_size;
    }

static gsize Viterbi_traceback_memory_size(Viterbi *viterbi,
                                           Region *region){
    return Matrix3d_size(region->query_length+1,
                         region->target_length+1,
                         viterbi->model->state_list->len,
                         sizeof(C4_Transition*));
    }

gboolean Viterbi_use_reduced_space(Viterbi *viterbi, Region *region){
    register gsize row_memory = Viterbi_get_row_size(viterbi, region);
    register gsize traceback_memory
                 = Viterbi_traceback_memory_size(viterbi, region);
    register gsize memory_limit = viterbi->vas->traceback_memory_limit
                                << 20;
    g_assert(Region_is_valid(region));
    if(region->query_length
      <= (viterbi->model->max_query_advance * 6))
        return FALSE;
    if(region->target_length
      <= (viterbi->model->max_target_advance * 6))
        return FALSE;
    if((!row_memory) || (!traceback_memory)){ /* overflow */
        /* Region_print(region, "using redspace (overflow)"); */
        return TRUE;
        }
    if((row_memory  + traceback_memory ) > memory_limit){
        /* Region_print(region, "using redspace (memory limit)"); */
        return TRUE;
        }
    return FALSE;
    }

/**/

static Viterbi_Row *Viterbi_Row_create(Viterbi *viterbi,
                                       Region *region){
    register Viterbi_Row *vr = g_new0(Viterbi_Row, 1);
    register gint i, j, k;
    g_assert(region);
    vr->region_start_query_id = -1;
    vr->region_start_target_id = -1;
    vr->checkpoint_id = -1;
    vr->cell_size = 1 + viterbi->model->total_shadow_designations;
    if(viterbi->mode == Viterbi_Mode_FIND_REGION){
        if(viterbi->model->start_state->scope != C4_Scope_CORNER){
            if(viterbi->model->start_state->scope != C4_Scope_QUERY)
                vr->region_start_query_id = vr->cell_size++;
            if(viterbi->model->start_state->scope != C4_Scope_TARGET)
                vr->region_start_target_id = vr->cell_size++;
            }
        }
    if(viterbi->mode == Viterbi_Mode_FIND_CHECKPOINTS)
        vr->checkpoint_id = vr->cell_size++;
    g_assert(vr->cell_size == viterbi->cell_size);
    vr->score_matrix = (C4_Score****)
             Matrix4d_create(viterbi->model->max_target_advance+1,
                             region->query_length+1,
                             viterbi->model->state_list->len,
                             vr->cell_size,
                             sizeof(C4_Score));
    for(i = 0; i < viterbi->model->max_target_advance+1; i++)
        for(j = 0; j < region->query_length+1; j++)
            for(k = 0; k < viterbi->model->state_list->len; k++)
                vr->score_matrix[i][j][k][0] = C4_IMPOSSIBLY_LOW_SCORE;
    return vr;
    }

static void Viterbi_Row_destroy(Viterbi_Row *vr){
    g_free(vr->score_matrix);
    g_free(vr);
    return;
    }

/**/

static gsize Viterbi_Row_get_size(Viterbi *viterbi, Region *region){
    register gint mat_size;
    mat_size = Matrix4d_size(viterbi->model->max_target_advance+1,
                             region->query_length+1,
                             viterbi->model->state_list->len,
                             viterbi->cell_size,
                             sizeof(C4_Score));
    if(!mat_size)
        return 0;
    return sizeof(Viterbi_Row) + mat_size;
    }

static gint Viterbi_checkpoint_rows(Viterbi *viterbi, Region *region){
    register gsize row_memory = Viterbi_Row_get_size(viterbi, region);
    register gint avail_rows =
        ((viterbi->vas->traceback_memory_limit << 20)/row_memory) - 1;
    register gint max_rows = (region->target_length
                           / (viterbi->model->max_target_advance << 1))
                           - 2;
    g_assert(max_rows > 0);
    if(avail_rows < 1)
        return 1;
    return MIN(avail_rows, max_rows);
    }

static C4_Transition ****Viterbi_traceback_memory_create(
                         Viterbi *viterbi, Region *region){
    return (C4_Transition****)Matrix3d_create(
                               region->query_length+1,
                               region->target_length+1,
                               viterbi->model->state_list->len,
                               sizeof(C4_Transition*));
    }

/**/

static Viterbi_Checkpoint *Viterbi_Checkpoint_create(Viterbi *viterbi,
                                                     Region *region,
                                                     Viterbi_Row *vr){
    register Viterbi_Checkpoint *vc = g_new0(Viterbi_Checkpoint, 1);
    register gint cp_count = Viterbi_checkpoint_rows(viterbi, region);
    register gint i;
    register C4_Score ****cp;
    g_assert(cp_count > 0);
    vc->checkpoint_list = g_ptr_array_new();
    vc->cell_size = vr->cell_size;
    for(i = 0; i < cp_count; i++){
        cp = (C4_Score****)Matrix4d_create(
                             viterbi->model->max_target_advance+1,
                             region->query_length+1,
                             viterbi->model->state_list->len,
                             vc->cell_size,
                             sizeof(C4_Score));
        g_ptr_array_add(vc->checkpoint_list, cp);
        }
    vc->section_length = region->target_length
                       / (vc->checkpoint_list->len+1);
    vc->counter = 0;
    return vc;
    }

static void Viterbi_Checkpoint_destroy(Viterbi_Checkpoint *vc){
    register gint i;
    for(i = 0; i < vc->checkpoint_list->len; i++)
        g_free(vc->checkpoint_list->pdata[i]);
    g_ptr_array_free(vc->checkpoint_list, TRUE);
    g_free(vc);
    return;
    }

/**/

static Viterbi_Continuation *Viterbi_Continuation_create(
                  C4_State *first_state, C4_Score *first_cell,
                  C4_State *final_state, C4_Score *final_cell){
    register Viterbi_Continuation *vc = g_new(Viterbi_Continuation, 1);
    g_assert(first_state);
    g_assert(first_cell);
    g_assert(final_state);
    g_assert(final_cell);
    vc->first_state = first_state;
    vc->first_cell  = first_cell;
    vc->final_state = final_state;
    vc->final_cell  = final_cell;
    return vc;
    }

static void Viterbi_Continuation_destroy(Viterbi_Continuation *vc){
    g_assert(vc);
    g_free(vc);
    return;
    }

/**/

void Viterbi_Data_set_continuation(Viterbi_Data *vd,
              C4_State *first_state, C4_Score *first_cell,
              C4_State *final_state, C4_Score *final_cell){
    g_assert(!vd->continuation);
    vd->continuation
        = Viterbi_Continuation_create(first_state, first_cell,
                                      final_state, final_cell);
    return;
    }

void Viterbi_Data_clear_continuation(Viterbi_Data *vd){
    if(vd->continuation){
        Viterbi_Continuation_destroy(vd->continuation);
        vd->continuation = NULL;
        }
    return;
    }

/**/

Viterbi_Data *Viterbi_Data_create(Viterbi *viterbi, Region *region){
    register Viterbi_Data *vd = g_new0(Viterbi_Data, 1);
    g_assert(viterbi);
    g_assert(region);
    vd->vr = Viterbi_Row_create(viterbi, region);
    if(viterbi->mode == Viterbi_Mode_FIND_REGION)
        vd->alignment_region = Region_copy(region);
    if(viterbi->mode == Viterbi_Mode_FIND_PATH)
        vd->traceback = Viterbi_traceback_memory_create(viterbi,
                                                        region);
    if(viterbi->mode == Viterbi_Mode_FIND_CHECKPOINTS)
        vd->checkpoint = Viterbi_Checkpoint_create(viterbi,
                                                   region, vd->vr);
    return vd;
    }

void Viterbi_Data_destroy(Viterbi_Data *vd){
    g_assert(vd);
    if(vd->alignment_region)
        Region_destroy(vd->alignment_region);
    if(vd->traceback)
        g_free(vd->traceback);
    if(vd->checkpoint)
        Viterbi_Checkpoint_destroy(vd->checkpoint);
    Viterbi_Data_clear_continuation(vd);
    Viterbi_Row_destroy(vd->vr);
    g_free(vd);
    return;
    }

/**/

Alignment *Viterbi_Data_create_Alignment(Viterbi_Data *vd,
                  C4_Model *model, C4_Score score, Region *region){
    register Alignment *alignment;
    register GPtrArray *transition_path = g_ptr_array_new();
    register gint i = vd->curr_query_end, j = vd->curr_target_end;
    register Region *alignment_region;
    register C4_Transition *transition;
    g_assert(vd->traceback);
    if(vd->continuation)
        transition = vd->traceback[i][j]
                                  [vd->continuation->final_state->id];
    else
        transition = vd->traceback[i][j]
                                  [model->end_state->state->id];
    g_assert(transition);
    do {
        g_ptr_array_add(transition_path, transition);
        i -= transition->advance_query;
        j -= transition->advance_target;
        transition = vd->traceback[i][j][transition->input->id];
        if(!transition)
            break;
        if(transition->input == model->start_state->state){
            g_ptr_array_add(transition_path, transition);
            i -= transition->advance_query;
            j -= transition->advance_target;
            break;
            }
        if(vd->continuation){
            if(!(i|j)){
                if(transition->output->id ==
                   vd->continuation->first_state->id){
                    break; /* continuation DP exit */
                    }
                }
            }
        g_assert(transition);
    } while(TRUE);
    alignment_region = Region_create(region->query_start + i,
                                     region->target_start + j,
                                     vd->curr_query_end - i,
                                     vd->curr_target_end - j);
    alignment = Alignment_create(model, alignment_region, score);
    for(i = transition_path->len-1; i >= 0; i--){ /* Reverse order */
        transition = transition_path->pdata[i];
        Alignment_add(alignment, transition, 1);
        }
    g_ptr_array_free(transition_path, TRUE);
    Region_destroy(alignment_region);
    return alignment;
    }

/**/

static void Viterbi_Row_shadow_start(Viterbi_Row *vr,
                C4_Model *model, C4_Transition *transition,
                Region *region,
                C4_Score *src, gint query_pos, gint target_pos,
                gpointer user_data){
    register gint i;
    register C4_Shadow *shadow;
    if(transition->input == model->start_state->state){
        if(vr->region_start_query_id != -1){
            src[vr->region_start_query_id]
                = query_pos - transition->advance_query;
            }
        if(vr->region_start_target_id != -1){
            src[vr->region_start_target_id]
                = target_pos - transition->advance_target;
            }
        }
    /* Call start funcs from src */
    for(i = 0; i < transition->input->src_shadow_list->len; i++){
        shadow = transition->input->src_shadow_list->pdata[i];
        src[shadow->designation+1] = shadow->start_func(
             region->query_start
             + query_pos - transition->advance_query,
             region->target_start
             + target_pos - transition->advance_target,
         user_data);
        }
    return;
    }

static void Viterbi_Row_shadow_end(Viterbi_Row *vr,
                C4_Model *model, C4_Transition *transition,
                Region *region, C4_Score *src,
                gint query_pos, gint target_pos,
                gpointer user_data){
    register gint i;
    register C4_Shadow *shadow;
    for(i = 0; i < transition->dst_shadow_list->len; i++){
        shadow = transition->dst_shadow_list->pdata[i];
        shadow->end_func(src[shadow->designation+1],
            region->query_start
            + query_pos - transition->advance_query,
            region->target_start
            + target_pos - transition->advance_target,
            user_data);
        }
    return;
    }

static void Viterbi_Data_assign(Viterbi *viterbi, Viterbi_Data *vd,
                    C4_Score *src, C4_Score *dst, C4_Score score,
                    gint query_pos, gint target_pos,
                    C4_Transition *transition, Region *region,
                    gpointer user_data){
    register gint i;
    dst[0] = score;
    /* Call shadow start on src */
    Viterbi_Row_shadow_start(vd->vr, viterbi->model, transition,
                           region, src, query_pos, target_pos,
                           user_data);
    for(i = 1; i < vd->vr->cell_size; i++) /* Shadow transport */
        dst[i] = src[i];
    if(viterbi->mode == Viterbi_Mode_FIND_PATH)
        vd->traceback[query_pos][target_pos]
                     [transition->output->id] = transition;
    return;
    }

static void Viterbi_Data_register_end(Viterbi_Data *vd,
                                      C4_Score *cell,
                                      gint query_pos, gint target_pos){
    vd->curr_query_end = query_pos;
    vd->curr_target_end = target_pos;
    if(vd->vr->region_start_query_id != -1){
        vd->curr_query_start
            = cell[vd->vr->region_start_query_id];
        }
    if(vd->vr->region_start_target_id != -1){
        vd->curr_target_start
            = cell[vd->vr->region_start_target_id];
        }
    return;
    }

/**/

static Viterbi_SubAlignment *Viterbi_SubAlignment_create(
                         gint query_start, gint target_start,
                         gint query_length, gint target_length,
                         C4_State *first_state, C4_Score *final_cell,
                         gint cell_length){
    register Viterbi_SubAlignment *vsa = g_new(Viterbi_SubAlignment, 1);
    register gint i;
    vsa->region = Region_create(query_start, target_start,
                                query_length, target_length);
    vsa->first_state = first_state;
    vsa->final_cell = g_new(C4_Score, cell_length);
    for(i = 0; i < cell_length; i++)
        vsa->final_cell[i] = final_cell[i];
    return vsa;
    }

void Viterbi_SubAlignment_destroy(Viterbi_SubAlignment *vsa){
    Region_destroy(vsa->region);
    g_free(vsa->final_cell);
    g_free(vsa);
    return;
    }

/**/

typedef struct {
    C4_State *state;
        gint  row;
        gint  pos;
} Viterbi_Checkpoint_SRP;
/* This encoding is safe because we can address the checkpoint matrix.
 */

static gint Viterbi_Checkpoint_SRP_encode(C4_Model *model,
                                          gint state_id,
                                          gint row_id, gint query_pos){
    return ( ( (query_pos * model->state_list->len)
              + state_id)
            * model->max_target_advance)
           + row_id;
    }

static void Viterbi_Checkpoint_SRP_decode(C4_Model *model, gint srp,
            Viterbi_Checkpoint_SRP *result){
    register gint remainder;
    result->row = srp % model->max_target_advance;
    remainder = srp / model->max_target_advance;
    result->state = model->state_list->pdata
                   [remainder % model->state_list->len];
    result->pos = remainder / model->state_list->len;
    g_assert(srp == Viterbi_Checkpoint_SRP_encode(model,
                result->state->id, result->row, result->pos));
    return;
    }

GPtrArray *Viterbi_Checkpoint_traceback(Viterbi *viterbi,
        Viterbi_Data *vd, Region *region,
        C4_State *first_state, C4_Score *final_cell){
    register gint prev_row;
    register Viterbi_SubAlignment *vsa;
    register gint query_start, target_start;
    register C4_Score cp_srp;
    register C4_Score ****checkpoint, *cell = final_cell;
    register gint i;
    register GPtrArray *vsa_list = g_ptr_array_new();
    Viterbi_Checkpoint_SRP srp;
    g_assert(vd->checkpoint);
    Viterbi_Checkpoint_SRP_decode(viterbi->model,
                                  vd->checkpoint->last_srp, &srp);
    query_start = region->query_start + srp.pos;
    target_start = region->target_start
                 + (vd->checkpoint->section_length
                   * vd->checkpoint->checkpoint_list->len)
                 - srp.row;
    vsa = Viterbi_SubAlignment_create(
                   query_start,
                   target_start,
                   Region_query_end(region)-query_start,
                   Region_target_end(region)-target_start,
                   srp.state,
                   final_cell, vd->checkpoint->cell_size);
    g_assert(Region_is_within(region, vsa->region));
    g_ptr_array_add(vsa_list, vsa);
    prev_row = srp.row;
    for(i = vd->checkpoint->checkpoint_list->len-1; i >= 1; i--){
        checkpoint = vd->checkpoint->checkpoint_list->pdata[i];
        prev_row = srp.row;
        cp_srp = checkpoint[prev_row]
                           [vsa->region->query_start
                           -region->query_start]
                           [vsa->first_state->id]
                           [vd->checkpoint->cell_size-1];
        Viterbi_Checkpoint_SRP_decode(viterbi->model, cp_srp, &srp);
        query_start = region->query_start + srp.pos;
        target_start = vsa->region->target_start
                     - vd->checkpoint->section_length
                     - srp.row + prev_row;
        cell = checkpoint[prev_row]
                         [vsa->region->query_start-region->query_start]
                         [vsa->first_state->id];
        vsa = Viterbi_SubAlignment_create(query_start, target_start,
                vsa->region->query_start-query_start,
                vsa->region->target_start-target_start,
                srp.state, cell, vd->checkpoint->cell_size);
        g_assert(Region_is_within(region, vsa->region));
        g_ptr_array_add(vsa_list, vsa);
        }
    /* Set checkpoint */
    checkpoint = vd->checkpoint->checkpoint_list->pdata[0];
    /* Set final cell */
    cell = checkpoint[srp.row]
                     [vsa->region->query_start-region->query_start]
                     [vsa->first_state->id];
    vsa = Viterbi_SubAlignment_create(region->query_start,
                                      region->target_start,
            query_start-region->query_start,
            target_start-region->target_start,
            first_state, cell, vd->checkpoint->cell_size);
    g_ptr_array_add(vsa_list, vsa);
    return vsa_list;
    }

/**/

static void Viterbi_Checkpoint_process(Viterbi_Checkpoint *vc,
                C4_Model *model, Region *region, gint target_pos,
                C4_Score ****prev_row){
    register gint i, j, k, l;
    register C4_Score ****checkpoint;
    if((!(target_pos % vc->section_length)) && target_pos){
        if(vc->counter < vc->checkpoint_list->len){
            checkpoint = vc->checkpoint_list->pdata[vc->counter++];
            for(i = 0; i < model->max_target_advance; i++){
                for(j = 0; j <= region->query_length; j++){
                    for(k = 0; k < model->state_list->len; k++){
                        /* Copy cell to checkpoint */
                        for(l = 0; l < vc->cell_size; l++){
                            checkpoint[i][j][k][l]
                            = prev_row[i][j][k][l];
                            }
                        /* Load state/pos to last shadow */
                        prev_row[i][j][k][vc->cell_size-1]
                            = Viterbi_Checkpoint_SRP_encode(model,
                                                            k, i, j);
                        }
                    }
                }
            }
        }
    return;
    }

static void Viterbi_Data_finalise(Viterbi_Data *vd, Region *region){
    if(vd->alignment_region){
        if(vd->vr->region_start_query_id != -1)
            vd->alignment_region->query_start = vd->curr_query_start
                                              + region->query_start;
        if(vd->vr->region_start_target_id != -1)
            vd->alignment_region->target_start = vd->curr_target_start
                                               + region->target_start;
        vd->alignment_region->query_length
                            = (vd->curr_query_end
                             - vd->curr_query_start);
        vd->alignment_region->target_length
                            = (vd->curr_target_end
                             - vd->curr_target_start);
        /* vd->alignment_region is in the coordinates of the sequences,
         * not of the current region
         */
        g_assert(Region_is_valid(vd->alignment_region));
        }
    return;
    }

static C4_Score Viterbi_interpreted(Viterbi *viterbi, Region *region,
                                    Viterbi_Data *vd,
                                    SubOpt_Index *soi,
                                    gpointer user_data){
    register C4_Score t, score = C4_IMPOSSIBLY_LOW_SCORE;
    register C4_Score ***swap_row, *src, *dst;
    register C4_Score ****prev_row
            = g_new(C4_Score***, viterbi->model->max_target_advance+1);
    register gint i, j, k, l;
    register C4_Transition *transition;
    register C4_Calc *calc;
    register gboolean *state_is_set = g_new(gboolean,
                                            viterbi->model->state_list->len);
    register gboolean end_is_set = FALSE;
    register C4_State *final_state = NULL;
    register C4_Score *dummy_start = g_new(C4_Score, vd->vr->cell_size),
                      *tmp_start;
    /**/
    g_assert(!viterbi->model->is_open);
    if(viterbi->model->init_func)
        viterbi->model->init_func(region, user_data);
    if(vd->continuation)
        final_state = vd->continuation->final_state;
    else
        final_state = viterbi->model->end_state->state;
    for(i = 0; i < viterbi->model->calc_list->len; i++){
        calc = viterbi->model->calc_list->pdata[i];
        if(calc && calc->init_func)
            calc->init_func(region, user_data);
        }
    /**/
    for(i = 0; i <= viterbi->model->max_target_advance; i++)
        prev_row[i] = vd->vr->score_matrix[i];
    for(j = 0; j <= region->target_length; j++){
        SubOpt_Index_set_row(soi, j);
        for(i = 0; i <= region->query_length; i++){
            for(k = 0; k < viterbi->model->state_list->len; k++){
                state_is_set[k] = FALSE;
                prev_row[0][i][k][0] = C4_IMPOSSIBLY_LOW_SCORE;
                }
            for(k = 0; k < viterbi->model->transition_list->len; k++){
                transition = viterbi->model->transition_list->pdata[k];
                if(!Layout_is_transition_valid(viterbi->layout,
                           viterbi->model, transition, i, j,
                           region->query_length, region->target_length))
                    continue;
                if(C4_Transition_is_match(transition)
                && soi
                && SubOpt_Index_is_blocked_fast(soi, i))
                    continue;
                if(vd->continuation
                && (transition->input
                    == viterbi->model->start_state->state)){
                    src = vd->continuation->first_cell;
                    dst = prev_row[0][0]
                                  [vd->continuation->first_state->id];
                    for(l = 0; l < vd->vr->cell_size; l++)
                        dst[l] = src[l];
                    state_is_set[vd->continuation->first_state->id] = TRUE;
                    }
                src = prev_row[transition->advance_target]
                              [i-transition->advance_query]
                              [transition->input->id];
                dst = prev_row[0][i][transition->output->id];
                /**/
                t = 0;
                if(transition->input
                   == viterbi->model->start_state->state){
                    if(vd->continuation){
                        src = prev_row[0][0]
                              [viterbi->model->start_state->state->id];
                        t = src[0];
                    } else {
                        if(viterbi->model
                           ->start_state->cell_start_func){
                            tmp_start = viterbi->model
                                ->start_state->cell_start_func(
                                region->query_start
                                + i - transition->advance_query,
                                region->target_start
                                + j - transition->advance_target,
                                user_data);
                            for(l = 0; l < vd->vr->cell_size; l++)
                                dummy_start[l] = tmp_start[l];
                            src = dummy_start;
                            t = src[0];
                            }
                        }
                } else {
                    t = src[0];
                    }
                /* Call shadow end from src */
                Viterbi_Row_shadow_end(vd->vr, viterbi->model,
                           transition, region, src, i, j, user_data);
                t += C4_Calc_score(transition->calc,
                     region->query_start+i-transition->advance_query,
                     region->target_start+j-transition->advance_target,
                         user_data);
                /**/
                if(transition->calc){
                    if(transition->calc->protect
                    & C4_Protect_UNDERFLOW){
                        if(t < C4_IMPOSSIBLY_LOW_SCORE)
                            t = C4_IMPOSSIBLY_LOW_SCORE;
                        }
                    if(transition->calc->protect
                    & C4_Protect_OVERFLOW){
                        if(t > C4_IMPOSSIBLY_HIGH_SCORE)
                            t = C4_IMPOSSIBLY_HIGH_SCORE;
                        }
                    }
                if(state_is_set[transition->output->id]){
                    if(dst[0] < t){
                        Viterbi_Data_assign(viterbi, vd, src, dst,
                            t, i, j, transition, region, user_data);
                        }
                } else {
                    state_is_set[transition->output->id] = TRUE;
                    Viterbi_Data_assign(viterbi, vd, src, dst,
                        t, i, j, transition, region, user_data);
                    }
                }
            /* corner */
            if(state_is_set[viterbi->model->end_state->state->id]){
                t = prev_row[0][i][final_state->id][0];
                if(end_is_set){
                    if(score < t){
                        score = t;
                        Viterbi_Data_register_end(vd,
                            prev_row[0][i][final_state->id], i, j);
                        }
                } else {
                    score = t;
                    end_is_set = TRUE;
                    Viterbi_Data_register_end(vd,
                        prev_row[0][i][final_state->id], i, j);
                    }
                if(viterbi->model->end_state->cell_end_func){
                    viterbi->model->end_state->cell_end_func(
                        prev_row[0][i][viterbi->model->end_state->state->id],
                        vd->vr->cell_size,
                        region->query_start+i, region->target_start+j,
                        user_data);
                    }
                }
            }
        if(viterbi->mode == Viterbi_Mode_FIND_CHECKPOINTS)
            Viterbi_Checkpoint_process(vd->checkpoint, viterbi->model,
                                       region, j, prev_row);
        /* Rotate rows backwards */
        swap_row = prev_row[viterbi->model->max_target_advance];
        for(i = viterbi->model->max_target_advance; i > 0; i--)
            prev_row[i] = prev_row[i-1];
        prev_row[0] = swap_row;
        }
    /**/
    g_assert(end_is_set);
    Viterbi_Data_finalise(vd, region);
    if(viterbi->mode == Viterbi_Mode_FIND_CHECKPOINTS){
        vd->checkpoint->last_srp
            = prev_row[1]
                      [region->query_length]
                      [final_state->id]
                      [vd->vr->cell_size-1];
        }
    /**/
    for(i = 0; i < viterbi->model->calc_list->len; i++){
        calc = viterbi->model->calc_list->pdata[i];
        if(calc && calc->exit_func)
            calc->exit_func(region, user_data);
        }
    if(viterbi->model->exit_func)
        viterbi->model->exit_func(region, user_data);
    if(vd->continuation){ /* Record copy of final_cell */
        src = prev_row[1][region->query_length][final_state->id];
        for(i = 0; i < vd->vr->cell_size; i++)
            vd->continuation->final_cell[i] = src[i];
        }
    g_free(prev_row);
    g_free(state_is_set);
    g_free(dummy_start);
    return score;
    }
/* This is supposed to be slow, but robust.
 */
/* FIXME: tidy / split to call Viterbi_interpreted_transition()
 */
/* FIXME: optimisation: avoid doing any calc when input score is -INF
 *        (can do the same with codegen ?)
 */

C4_Score Viterbi_calculate(Viterbi *viterbi, Region *region,
                           Viterbi_Data *vd, gpointer user_data,
                           SubOpt *subopt){
    register C4_Score score;
    register SubOpt_Index *soi = NULL;
    g_assert(viterbi);
    g_assert(Region_is_valid(region));
    if(subopt)
        soi = SubOpt_Index_create(subopt, region);
    if(viterbi->func){ /* Use compiled version */
        score = viterbi->func(viterbi->model, region, vd, soi,
                              user_data);
    } else {
        score = Viterbi_interpreted(viterbi, region, vd, soi,
                                    user_data);
        }
    if(soi)
        SubOpt_Index_destroy(soi);
    return score;
    }

/**/

static gchar *Viterbi_expand_macro(gchar *macro,
              gint advance_query, gint advance_target,
              C4_Shadow *shadow, gchar *cell){
    register GString *str = g_string_sized_new(strlen(macro));
    register gchar *expanded_macro, *tmp;
    register gint i;
    for(i = 0; macro[i]; i++){
        if((macro[i] == '%')){
            i++;
            switch(macro[i]){
                case 'Q':
                    i++;
                    switch(macro[i]){
                        case 'S':
                            g_string_append(str, "region->query_start");
                            break;
                        case 'P':
                            g_string_append(str, "region->query_start+i-");
                            g_string_append_c(str, '0' + advance_query);
                            break;
                        case 'R':
                            g_string_append(str, "i-");
                            g_string_append_c(str, '0' + advance_query);
                            break;
                        case 'E':
                            g_string_append(str, "Region_query_end(region)");
                            break;
                        case 'L':
                            g_string_append(str, "region->query_length");
                            break;
                        default:
                            g_error("Bad query macro [%%Q%c] in [%s]",
                                    macro[i], macro);
                            break;
                        }
                    break;
                case 'T':
                    i++;
                    switch(macro[i]){
                        case 'S':
                            g_string_append(str, "region->target_start");
                            break;
                        case 'P':
                            g_string_append(str, "region->target_start+j-");
                            g_string_append_c(str, '0' + advance_target);
                            break;
                        case 'R':
                            g_string_append(str, "j-");
                            g_string_append_c(str, '0' + advance_target);
                            break;
                        case 'E':
                            g_string_append(str, "Region_target_end(region)");
                            break;
                        case 'L':
                            g_string_append(str, "region->target_length");
                            break;
                        default:
                            g_error("Bad target macro [%%T%c] in [%s]",
                                    macro[i], macro);
                            break;
                        }
                    break;
                case 'S':
                    i++;
                    switch(macro[i]){
                        case 'S':
                            g_assert(shadow);
                            tmp = g_strdup_printf("src[%d]",
                                                  shadow->designation+1);
                            g_string_append(str, tmp);
                            g_free(tmp);
                            break;
                        default:
                            g_error("Bad shadow macro [%%S%c] in [%s]",
                                    macro[i], macro);
                            break;
                        }
                    break;
                case 'C':
                    g_assert(cell);
                    g_string_append(str, cell);
                    break;
                case '%':
                    g_string_append(str, "%");
                    break;
                default:
                    g_error("Bad macro [%%%c] in [%s]",
                            macro[i], macro);
                    break;
                }
        } else {
            g_string_append_c(str, macro[i]);
            }
        }
    expanded_macro = str->str;
    g_string_free(str, FALSE);
    return expanded_macro;
    }
/*   Macro expansion rules:
 *   ---------------
 *   %%  -> %
 *   %QS -> query_start
 *   %QP -> query_position
 *   %QR -> query_position in region
 *   %QE -> query_end
 *   %QL -> query_length
 *   %TS -> target_start
 *   %TP -> target_position
 *   %TR -> target_position in region
 *   %TE -> target_end
 *   %TL -> target_length
 *   %SS -> shadow_score
 *   %C  -> cell
 */

/**/

static gboolean Viterbi_implement_require_shadow(Viterbi *viterbi){
    register gint i;
    register C4_Shadow *shadow;
    for(i = 0; i < viterbi->model->shadow_list->len; i++){
        shadow = viterbi->model->shadow_list->pdata[i];
        g_assert(shadow->start_func);
        if(!shadow->start_macro)
            return TRUE;
        g_assert(shadow->end_func);
        if(!shadow->end_macro)
            return TRUE;
        }
    return FALSE;
    }

static gboolean Viterbi_implement_require_transition(Viterbi *viterbi){
    register gint i;
    register C4_Calc *calc;
    if(viterbi->use_continuation
    || (viterbi->mode == Viterbi_Mode_FIND_PATH))
        return TRUE;
    for(i = 0; i < viterbi->model->calc_list->len; i++){
        calc = viterbi->model->calc_list->pdata[i];
        if(calc->calc_func && (!calc->calc_macro))
            return TRUE;
        }
    return FALSE;
    }

static gboolean Viterbi_implement_require_calc(Viterbi *viterbi){
    register gint i;
    register C4_Calc *calc;
    for(i = 0; i < viterbi->model->calc_list->len; i++){
        calc = viterbi->model->calc_list->pdata[i];
        if(calc->init_func && (!calc->init_macro))
            return TRUE;
        if(calc->exit_func && (!calc->exit_macro))
            return TRUE;
        }
    return FALSE;
    }

static void Viterbi_implement_shadow_start(Viterbi *viterbi,
                                           Codegen *codegen,
                                           C4_Transition *transition){
    register gint i;
    register C4_Shadow *shadow;
    register gchar *expanded_macro;
    if(transition->input == viterbi->model->start_state->state){
        if(viterbi->mode == Viterbi_Mode_FIND_REGION){
            if(viterbi->model->start_state->scope != C4_Scope_CORNER){
                if(viterbi->model->start_state->scope
                                                   != C4_Scope_QUERY)
                    Codegen_printf(codegen,
                        "src[vd->vr->region_start_query_id] = i-%d;\n",
                        transition->advance_query);
                if(viterbi->model->start_state->scope != C4_Scope_TARGET)
                    Codegen_printf(codegen,
                        "src[vd->vr->region_start_target_id] = j-%d;\n",
                        transition->advance_target);
                }
            }
        }
    /* Call start funcs from src */
    for(i = 0; i < transition->input->src_shadow_list->len; i++){
        shadow = transition->input->src_shadow_list->pdata[i];
        if(shadow->start_macro){
            expanded_macro = Viterbi_expand_macro(shadow->start_macro,
                                 transition->advance_query,
                                 transition->advance_target,
                                 NULL, NULL);
            Codegen_printf(codegen, "src[%d] = %s;\n",
                    shadow->designation+1, expanded_macro);
            g_free(expanded_macro);
        } else {
            g_assert(shadow->start_func);
            Codegen_printf(codegen, "shadow = "
               "transition->input->output_shadow_list->pdata[%d];\n",
               i);
            Codegen_printf(codegen,
                 "src[%d] = shadow->start_func("
                 "region->query_start+i-%d, region->target_start+j-%d,"
                 " user_data);\n",
                 shadow->designation+1,
                 transition->advance_query,
                 transition->advance_target);
            }
        }
    return;
    }

static void Viterbi_implement_shadow_end(Viterbi *viterbi,
                                         Codegen *codegen,
                                         C4_Transition *transition){
    register gint i;
    register C4_Shadow *shadow;
    register gchar *expanded_macro;
    for(i = 0; i < transition->dst_shadow_list->len; i++){
        shadow = transition->dst_shadow_list->pdata[i];
        if(shadow->end_macro){
            expanded_macro = Viterbi_expand_macro(shadow->end_macro,
                                 transition->advance_query,
                                 transition->advance_target,
                                 shadow, NULL);
            Codegen_printf(codegen, "%s;\n", expanded_macro);
            g_free(expanded_macro);
        } else {
            g_assert(shadow->end_func);
            Codegen_printf(codegen, "shadow = "
               "transition->input->input_shadow_list->pdata[%d];\n",
               i);
            Codegen_printf(codegen,
                 "shadow->end_func(src[%d],"
                 " region->query_start+i-%d, region->target_start+j-%d,"
                 " user_data);\n",
                 shadow->designation+1,
                 transition->advance_query, transition->advance_target);
            }
        }
    return;
    }

static void Viterbi_implement_dp_set_cell(Viterbi *viterbi,
                    Codegen *codegen, C4_Transition *transition,
                    gint cell_size){
    register gint i;
    /* Call shadow start on src */
    Viterbi_implement_shadow_start(viterbi, codegen, transition);
    Codegen_printf(codegen, "dst[0] = t;\n");
    for(i = 1; i < cell_size; i++)
        Codegen_printf(codegen, "dst[%d] = src[%d];\n", i, i);
    /* FIXME: optimisation:
     *        do not need to do copy for transitions
     *        when the shadow is out of scope.
     */
    if(viterbi->mode == Viterbi_Mode_FIND_PATH){
        Codegen_printf(codegen,
                "vd->traceback[i][j][%d] = transition;\n",
                transition->output->id);
        }
    return;
    }

static void Viterbi_implement_register_end(Viterbi *viterbi,
                                           Codegen *codegen, gchar *end_cell){
    if((viterbi->mode == Viterbi_Mode_FIND_REGION)
    || (viterbi->mode == Viterbi_Mode_FIND_PATH)){
        Codegen_printf(codegen, "vd->curr_query_end = i;\n");
        Codegen_printf(codegen, "vd->curr_target_end = j;\n");
        }
    if(viterbi->mode == Viterbi_Mode_FIND_REGION){
        if(viterbi->model->start_state->scope != C4_Scope_CORNER){
            if(viterbi->model->start_state->scope
                                               != C4_Scope_QUERY){
                Codegen_printf(codegen,
                        "vd->curr_query_start = "
                        "%s[vd->vr->region_start_query_id];\n", end_cell);
                }
            if(viterbi->model->start_state->scope
                                               != C4_Scope_TARGET){
                Codegen_printf(codegen,
                        "vd->curr_target_start = "
                        "%s[vd->vr->region_start_target_id];\n", end_cell);
                }
            }
        }
    return;
    }

static void Viterbi_implement_dp_cell(Viterbi *viterbi,
                Codegen *codegen, Layout_Mask *mask, gint cell_size){
    register gint i, j;
    register C4_Transition *transition;
    register gchar *end_cell, *expanded_macro;
    register gboolean is_valid, t_is_set;
    g_assert(mask);
    /**/
    /* Zero dst cell */
    for(i = 0; i < viterbi->model->state_list->len; i++)
        Codegen_printf(codegen, "prev_row_0[i][%d][0] = %d;\n",
                       i, C4_IMPOSSIBLY_LOW_SCORE);
    for(i = 0; i < viterbi->model->state_list->len; i++)
        Codegen_printf(codegen, "state_is_set[%d] = FALSE;\n", i);
    for(i = 0; i < viterbi->model->transition_list->len; i++){
        is_valid = g_array_index(mask->transition_mask, gboolean, i);
        if(!is_valid)
            continue;
        transition = viterbi->model->transition_list->pdata[i];
        g_assert(transition->id == i);
        /* Open SubOpt conditional */
        if(C4_Transition_is_match(transition)){
            Codegen_printf(codegen,
              "if(!(soi && SubOpt_Index_is_blocked_fast(soi, i))){\n");
            Codegen_indent(codegen, 1);
            }
        if(transition->input == viterbi->model->start_state->state){
            if(viterbi->use_continuation){
                Codegen_printf(codegen,
                    "src = vd->continuation->first_cell;\n");
                Codegen_printf(codegen,
                    "dst = prev_row_0[0]"
                    "[vd->continuation->first_state->id];\n");
                for(j = 0; j < cell_size; j++)
                    Codegen_printf(codegen,
                                  "dst[%d] = src[%d];\n", j, j);
                Codegen_printf(codegen,
                    "state_is_set[vd->continuation->first_state->id] = TRUE;\n");
                }
            }
        Codegen_printf(codegen, "/* transition [%s] */\n",
                      transition->name);
        Codegen_printf(codegen, "src = prev_row_%d[i-%d][%d];\n",
                      transition->advance_target,
                      transition->advance_query,
                      transition->input->id);
        Codegen_printf(codegen, "dst = prev_row_0[i][%d];\n",
                      transition->output->id);
        t_is_set = FALSE;
        if(Viterbi_implement_require_transition(viterbi))
            Codegen_printf(codegen,
                "transition = model->transition_list->pdata[%d];\n",
                transition->id);
        /* Call shadow end from dst */
        Viterbi_implement_shadow_end(viterbi, codegen, transition);
        if(transition->calc){
            if(transition->calc->calc_func){
                if(transition->calc->calc_macro){
                    expanded_macro = Viterbi_expand_macro(
                                     transition->calc->calc_macro,
                                     transition->advance_query,
                                     transition->advance_target,
                                     NULL, NULL);
                    Codegen_printf(codegen, "t %s= %s;\n",
                                  t_is_set?"+":"",expanded_macro);
                    t_is_set = TRUE;
                    g_free(expanded_macro);
                } else {
                    Codegen_printf(codegen,
                        "t %s= transition->calc->calc_func(\n"
                        "region->query_start+i-%d,\n"
                        "region->target_start+j-%d,\n"
                        "user_data);\n",
                        t_is_set?"+":"",
                        transition->advance_query,
                        transition->advance_target);
                    t_is_set = TRUE;
                    }
            } else {
                if(transition->calc->max_score){
                    Codegen_printf(codegen, "t %c= %d;\n",
                                  t_is_set?"+":"",
                                  transition->calc->max_score);
                    t_is_set = TRUE;
                    }
                }
            }
        if(transition->input == viterbi->model->start_state->state){
            if(viterbi->use_continuation){
                Codegen_printf(codegen, "t %s= prev_row_0[0][%d][0];\n",
                    t_is_set?"+":"",
                    viterbi->model->start_state->state->id);
                t_is_set = TRUE;
            } else {
                if(viterbi->model->start_state->cell_start_func){
                    if(viterbi->model->start_state->cell_start_macro){
                        expanded_macro = Viterbi_expand_macro(
                          viterbi->model->start_state->cell_start_macro,
                          transition->advance_query,
                          transition->advance_target,
                          NULL, NULL);
                        Codegen_printf(codegen, "src = %s;\n",
                                      expanded_macro);
                        g_free(expanded_macro);
                    } else {
                        Codegen_printf(codegen,
                          "src = model->start_state->cell_start_func(\n"
                          "region->query_start+i-%d,\n"
                          "region->target_start+j-%d,\n"
                          "user_data);\n",
                           transition->advance_query,
                           transition->advance_target);
                        }
                    Codegen_printf(codegen, "t %s= src[0];\n",
                        t_is_set?"+":"");
                    t_is_set = TRUE;
                    }
                }
        } else {
            Codegen_printf(codegen, "t %s= src[0];\n", t_is_set?"+":"");
            t_is_set = TRUE;
            }
        /* Check against overflow and underflow */
        if(t_is_set){
            if(transition->calc){
                if(transition->calc->protect & C4_Protect_UNDERFLOW){
                    Codegen_printf(codegen,
                                  "if(t < C4_IMPOSSIBLY_LOW_SCORE)\n");
                    Codegen_indent(codegen, 1);
                    Codegen_printf(codegen,
                                  "t = C4_IMPOSSIBLY_LOW_SCORE;\n");
                    Codegen_indent(codegen, -1);
                    }
                if(transition->calc->protect & C4_Protect_OVERFLOW){
                    Codegen_printf(codegen,
                                  "if(t > C4_IMPOSSIBLY_HIGH_SCORE)\n");
                    Codegen_indent(codegen, 1);
                    Codegen_printf(codegen,
                                  "t = C4_IMPOSSIBLY_HIGH_SCORE;\n");
                    Codegen_indent(codegen, -1);
                    }
                }
        } else {
            Codegen_printf(codegen, "t = 0;\n");
            }
        /* Challenge with t */
        Codegen_printf(codegen, "if(state_is_set[%d]){\n",
                       transition->output->id);
        Codegen_indent(codegen, 1);
        Codegen_printf(codegen, "if(dst[0] < t){\n");
        Codegen_indent(codegen, 1);
        Viterbi_implement_dp_set_cell(viterbi, codegen,
            transition, cell_size);
        Codegen_printf(codegen, "}\n");
        Codegen_indent(codegen, -2);
        Codegen_printf(codegen, "} else {\n");
        Codegen_indent(codegen, 1);
        Codegen_printf(codegen, "state_is_set[%d] = TRUE;\n",
                      transition->output->id);
        Viterbi_implement_dp_set_cell(viterbi, codegen,
            transition, cell_size);
        Codegen_printf(codegen, "}\n");
        Codegen_indent(codegen, -1);
        /* Close SubOpt conditional */
        if(C4_Transition_is_match(transition)){
            Codegen_printf(codegen, "}\n");
            Codegen_indent(codegen, -1);
            }
        }
    /* Take score if end set */
    Codegen_printf(codegen, "if(state_is_set[%d]){\n",
            viterbi->model->end_state->state->id);
    Codegen_indent(codegen, 1);
    if(viterbi->use_continuation){
        end_cell
            = g_strdup_printf("prev_row_0[i]"
                  "[vd->continuation->final_state->id]");
    } else {
        end_cell = g_strdup_printf("prev_row_0[i][%d]",
                         viterbi->model->end_state->state->id);
        }
    Codegen_printf(codegen, "if(end_is_set){\n");
    Codegen_indent(codegen, 1);
    Codegen_printf(codegen, "if(score < %s[0]){\n", end_cell);
    Codegen_indent(codegen, 1);
    Codegen_printf(codegen, "score = %s[0];\n", end_cell);
    Viterbi_implement_register_end(viterbi, codegen, end_cell);
    Codegen_printf(codegen, "}\n");
    Codegen_indent(codegen, -2);
    Codegen_printf(codegen, "} else {\n");
    Codegen_indent(codegen, 1);
    Codegen_printf(codegen, "score = %s[0];\n", end_cell);
    Codegen_printf(codegen, "end_is_set = TRUE;\n");
    Viterbi_implement_register_end(viterbi, codegen, end_cell);
    Codegen_printf(codegen, "}\n");
    Codegen_indent(codegen, -1);
    if(viterbi->model->end_state->cell_end_func){
        if(viterbi->model->end_state->cell_end_macro){
            expanded_macro = Viterbi_expand_macro(
              viterbi->model->end_state->cell_end_macro,
              0, 0, NULL, end_cell);
            Codegen_printf(codegen, "%s;\n", expanded_macro);
            g_free(expanded_macro);
        } else {
            Codegen_printf(codegen,
                "model->end_state->cell_end_func(\n"
                "%s, %d, region->query_start+i,\n"
                "region->target_start+j,\n"
                "user_data);\n",
                end_cell, cell_size);
            }
        }
    g_free(end_cell);
    Codegen_printf(codegen, "}\n");
    Codegen_indent(codegen, -1);
    return;
    }
/* FIXME: tidy / remove redundancy
 */
/* FIXME: optimisation
 *        - one-time calculation of {input,output}_{query,target}_pos
 *        - avoid unnecessary assignments
 *        - do not add onto zero scores
 *        - only protect scores after addition/subtraction
 *        - other things to simplify generated code
 */

/* FIXME: optimisation : need to find repeated rows and cells
 *        in the Layout (before max_advance_{query,target}),
 *        and remove resulting redundnacy in codegen.
 */

/* FIXME: optimisation : allow start/end transitions for global
 *        and semi-global models to be performed outside of main
 *        alignment (eg. for Find_Path_continuation() etc).
 *        Do this by removing corner start/end transitions from model,
 *        and performing them as part of the continuation DP.
 */

static void Viterbi_implement_opc(Viterbi *viterbi, Codegen *codegen,
                gint row_id, gint cell_id, gint cell_size){
    register Layout_Row *row = viterbi->layout->row_list->pdata[row_id];
    register Layout_Cell *cell = row->cell_list->pdata[cell_id];
    register gchar *next_row;
    /* The Keyword That Dare Not Speak Its Name: */
    gchar tktdnsin[5] = {0x67, 0x6F, 0x74, 0x6F, 0x00};
    Codegen_printf(codegen, "/* CELL [%d,%d] */\n", row_id, cell_id);
    Codegen_printf(codegen, "i++;\n");
    if(row_id >= (viterbi->layout->row_list->len-2))
        next_row = g_strdup("nth_row");
    else
        next_row = g_strdup_printf("row_%d", row_id+1);
    if(cell->end_target){
        Codegen_printf(codegen, "if(j == region->target_length){\n");
        Codegen_indent(codegen, 1);
        if(cell->corner){
            Codegen_printf(codegen, "if(i == region->query_length){\n");
            Codegen_indent(codegen, 1);
            Viterbi_implement_dp_cell(viterbi, codegen, cell->corner,
                cell_size);
            Codegen_printf(codegen, "%s end_of_dp;\n", tktdnsin);
            Codegen_indent(codegen, -1);
            Codegen_printf(codegen, "} else {\n");
            Codegen_indent(codegen, 1);
            Viterbi_implement_dp_cell(viterbi, codegen,
                                      cell->end_target, cell_size);
            Codegen_printf(codegen, "}\n");
            Codegen_indent(codegen, -1);
        } else {
            Viterbi_implement_dp_cell(viterbi, codegen,
                                      cell->end_target, cell_size);
            Codegen_printf(codegen, "if(i == region->query_length){\n");
            Codegen_indent(codegen, 1);
            Codegen_printf(codegen, "%s end_of_dp;\n", tktdnsin);
            Codegen_printf(codegen, "}\n");
            Codegen_indent(codegen, -1);
            }
        Codegen_indent(codegen, -1);
        Codegen_printf(codegen, "} else {\n");
        Codegen_indent(codegen, 1);
        if(cell->end_query){
            Codegen_printf(codegen, "if(i == region->query_length){\n");
            Codegen_indent(codegen, 1);
            Viterbi_implement_dp_cell(viterbi, codegen,
                                      cell->end_query, cell_size);
            Codegen_printf(codegen, "%s %s;\n", tktdnsin, next_row);
            Codegen_indent(codegen, -1);
            Codegen_printf(codegen, "} else {\n");
            Codegen_indent(codegen, 1);
            Viterbi_implement_dp_cell(viterbi, codegen, cell->normal,
                cell_size);
            Codegen_printf(codegen, "}\n");
            Codegen_indent(codegen, -1);
        } else {
            Viterbi_implement_dp_cell(viterbi, codegen,
                                      cell->normal, cell_size);
            Codegen_printf(codegen, "if(i == region->query_length){\n");
            Codegen_indent(codegen, 1);
            Codegen_printf(codegen, "%s %s;\n", tktdnsin, next_row);
            Codegen_printf(codegen, "}\n");
            Codegen_indent(codegen, -1);
            }
        Codegen_printf(codegen, "}\n");
        Codegen_indent(codegen, -1);
    } else { /* cell->end_target missing */
        if(cell->end_query){
            Codegen_printf(codegen,
                           "if(i == region->query_length){\n");
            Codegen_indent(codegen, 1);
            Viterbi_implement_dp_cell(viterbi, codegen,
                                      cell->end_query, cell_size);
            Codegen_printf(codegen,
                           "if(j == region->target_length){\n");
            Codegen_indent(codegen, 1);
            Codegen_printf(codegen, "%s end_of_dp;\n", tktdnsin);
            Codegen_indent(codegen, -1);
            Codegen_printf(codegen, "} else {\n");
            Codegen_indent(codegen, 1);
            Codegen_printf(codegen, "%s %s;\n", tktdnsin, next_row);
            Codegen_printf(codegen, "}\n");
            Codegen_indent(codegen, -1);
            Codegen_indent(codegen, -1);
            Codegen_printf(codegen, "} else {\n");
            Codegen_indent(codegen, 1);
            Viterbi_implement_dp_cell(viterbi, codegen,
                                      cell->normal, cell_size);
            Codegen_printf(codegen,
                          "if(j == region->target_length){\n");
            Codegen_indent(codegen, 1);
            Codegen_printf(codegen, "%s end_of_dp;\n", tktdnsin);
            Codegen_printf(codegen, "}\n");
            Codegen_indent(codegen, -1);
            Codegen_printf(codegen, "}\n");
            Codegen_indent(codegen, -1);
        } else {
            Viterbi_implement_dp_cell(viterbi, codegen,
                                      cell->normal, cell_size);
            Codegen_printf(codegen,
                           "if(j == region->target_length){\n");
            Codegen_indent(codegen, 1);
            Codegen_printf(codegen, "if(i == region->query_length){\n");
            Codegen_indent(codegen, 1);
            Codegen_printf(codegen, "%s end_of_dp;\n", tktdnsin);
            Codegen_printf(codegen, "}\n");
            Codegen_indent(codegen, -1);
            Codegen_indent(codegen, -1);
            Codegen_printf(codegen, "} else {\n");
            Codegen_indent(codegen, 1);
            Codegen_printf(codegen, "if(i == region->query_length){\n");
            Codegen_indent(codegen, 1);
            Codegen_printf(codegen, "%s %s;\n", tktdnsin, next_row);
            Codegen_printf(codegen, "}\n");
            Codegen_indent(codegen, -1);
            Codegen_printf(codegen, "}\n");
            Codegen_indent(codegen, -1);
            }
        }
    g_free(next_row);
    return;
    }

static void Viterbi_implement_dp_row(Viterbi *viterbi, Codegen *codegen,
               gint row_id, gint cell_size){
    register gint i, j, cp_row;
    register Layout_Row *row = viterbi->layout->row_list->pdata[row_id];
    /**/
    if(row_id){
        if(row_id == (viterbi->layout->row_list->len-1))
            Codegen_printf(codegen, "nth_row:\n");
        else
            Codegen_printf(codegen, "row_%d:\n", row_id);
        if(viterbi->mode == Viterbi_Mode_FIND_CHECKPOINTS){
            Codegen_printf(codegen,
                "if((!(j %% vd->checkpoint->section_length)) && j){\n");
            Codegen_indent(codegen, 1);
            Codegen_printf(codegen,
                "if(vd->checkpoint->counter < "
                "vd->checkpoint->checkpoint_list->len){\n");
            Codegen_indent(codegen, 1);
            Codegen_printf(codegen, "checkpoint = "
                "vd->checkpoint->checkpoint_list->pdata"
                "[vd->checkpoint->counter++];\n");
            for(cp_row = 0;
                cp_row < viterbi->model->max_target_advance; cp_row++){
                Codegen_printf(codegen,
                     "for(i = 0; i <= region->query_length; i++){\n");
                Codegen_indent(codegen, 1);
                Codegen_printf(codegen, "for(k = 0; k < %d; k++){\n",
                    viterbi->model->state_list->len);
                Codegen_indent(codegen, 1);
                /* Copy cell to checkpoint */
                for(j = 0; j < cell_size; j++)
                    Codegen_printf(codegen,
                        "checkpoint[%d][i][k][%d]"
                        " = prev_row_%d[i][k][%d];\n",
                        cp_row, j, cp_row, j);
                /* Load state/pos to last shadow */
                Codegen_printf(codegen,
                    "prev_row_%d[i][k][%d] = "
                    "(((i*%d)+k)*%d)+%d;\n",
                    cp_row, cell_size-1,
                    viterbi->model->state_list->len,
                    viterbi->model->max_target_advance,
                    cp_row);
                Codegen_printf(codegen, "}\n");
                Codegen_indent(codegen, -1);
                Codegen_printf(codegen, "}\n");
                Codegen_indent(codegen, -1);
                }
            Codegen_indent(codegen, -1);
            Codegen_printf(codegen, "}\n");
            Codegen_indent(codegen, -1);
            Codegen_printf(codegen, "}\n");
            }
        Codegen_printf(codegen, "j++;\n");
        /* Set SubOpt row */
        Codegen_printf(codegen, "SubOpt_Index_set_row(soi, j);\n");
        /* Rotate prev_row_%d pointers backwards */
        Codegen_printf(codegen, "swap_row = prev_row_%d;\n",
                viterbi->model->max_target_advance);
        for(i = viterbi->model->max_target_advance; i > 0; i--)
            Codegen_printf(codegen, "prev_row_%d = prev_row_%d;\n",
                          i, i-1);
        Codegen_printf(codegen, "prev_row_0 = swap_row;\n");
        }
    Codegen_printf(codegen, "/* start block */\n");
    Codegen_printf(codegen, "i = -1;\n");
    /* Start of cell */
    for(i = 0; i < (row->cell_list->len-1); i++){
        Codegen_printf(codegen, "/* explicit cell [%d] */\n", i);
        Viterbi_implement_opc(viterbi, codegen, row_id, i, cell_size);
        }
    /* Last cell */
    Codegen_printf(codegen, "/* last cell [%d] */\n", i);
    Codegen_printf(codegen, "do {\n");
    Codegen_indent(codegen, 1);
    Viterbi_implement_opc(viterbi, codegen, row_id, i, cell_size);
    Codegen_indent(codegen, -1);
    Codegen_printf(codegen, "} while(TRUE);\n");
    /* End of cell */
    return;
    }
/*    +--------+     +-----------+
 *    | NORMAL | <-- | QUERY END |
 *    +--------+     +-----------+
 *         ^
 *         |
 *    +------------+     +--------+
 *    | TARGET END | <-- | CORNER |
 *    +------------+     +--------+
 */

/**/

static void Viterbi_implement_finalise(Viterbi *viterbi,
                                       Codegen *codegen){
    if(viterbi->mode == Viterbi_Mode_FIND_REGION){
        if(viterbi->model->start_state->scope != C4_Scope_CORNER){
            if(viterbi->model->start_state->scope
                                               != C4_Scope_QUERY){
                Codegen_printf(codegen,
                  "vd->alignment_region->query_start"
                  " = vd->curr_query_start + region->query_start;\n");
                }
            if(viterbi->model->start_state->scope
                                               != C4_Scope_TARGET)
                Codegen_printf(codegen,
                  "vd->alignment_region->target_start"
                  " = vd->curr_target_start + region->target_start;\n");
                }
        Codegen_printf(codegen,
            "vd->alignment_region->query_length"
            " = vd->curr_query_end - vd->curr_query_start;\n");
        Codegen_printf(codegen,
            "vd->alignment_region->target_length"
            " = vd->curr_target_end - vd->curr_target_start;\n");
        /* vd->alignment_region is in the coordinates of the sequences,
         * not of the current region
         */
        }
    return;
    }

static void Viterbi_implement_dp(Viterbi *viterbi, Codegen *codegen){
    register gint i;
    register gint cell_size = Viterbi_get_cell_size(viterbi);
    Codegen_printf(codegen, "/* Implementing [%d] explicit rows */\n",
                  viterbi->layout->row_list->len);
    Codegen_printf(codegen, "C4_Score %s%s{\n",
                  viterbi->name, Viterbi_DP_Func_ARGS_STR);
    Codegen_indent(codegen, 1);
    Codegen_printf(codegen, "register C4_Score score"
                          " = C4_IMPOSSIBLY_LOW_SCORE;\n");
    /* Preliminary code */
    Codegen_printf(codegen, "register gint i, j = 0;\n");
    for(i = 0; i <= viterbi->model->max_target_advance; i++){
        Codegen_printf(codegen,
             "register C4_Score ***prev_row_%d"
             " = vd->vr->score_matrix[%d];\n",
              i, i);
        }
    Codegen_printf(codegen, "register C4_Score ***swap_row,"
                          " *src, *dst, t;\n");
    if(Viterbi_implement_require_shadow(viterbi)){
        Codegen_printf(codegen,
                      "register C4_Transition *transition;\n");
        Codegen_printf(codegen, "register C4_Shadow *shadow;\n");
    } else {
        if(Viterbi_implement_require_transition(viterbi))
            Codegen_printf(codegen,
                     "register C4_Transition *transition;\n");
        }
    if(Viterbi_implement_require_calc(viterbi))
        Codegen_printf(codegen, "register C4_Calc *calc;\n");
    if(viterbi->mode == Viterbi_Mode_FIND_CHECKPOINTS){
        Codegen_printf(codegen, "register C4_Score ****checkpoint;\n");
        Codegen_printf(codegen, "register gint k;\n");
        }
    Codegen_printf(codegen, "register gboolean *state_is_set\n"
                            " = g_new(gboolean, %d);\n",
                            viterbi->model->state_list->len);
    Codegen_printf(codegen, "register gboolean end_is_set = FALSE;\n");
    /* Local code must go between declarations and first block */
    for(i = 0; i < viterbi->model->local_code_list->len; i++)
        Codegen_printf(codegen,
                       viterbi->model->local_code_list->pdata[i]);
    CGUtil_prep(codegen, viterbi->model, TRUE);
    /* Set initial SubOpt row */
    Codegen_printf(codegen, "SubOpt_Index_set_row(soi, 0);\n");
    /* Start of viterbi */
    for(i = 0; i < (viterbi->layout->row_list->len-1); i++){
        Codegen_printf(codegen, "/* explicit row [%d] */\n", i);
        Viterbi_implement_dp_row(viterbi, codegen, i, cell_size);
        }
    /* Last row */
    Codegen_printf(codegen, "/* last row [%d] */\n", i);
    Codegen_printf(codegen, "do {\n");
    Codegen_indent(codegen, 1);
    Viterbi_implement_dp_row(viterbi, codegen, i, cell_size);
    Codegen_indent(codegen, -1);
    Codegen_printf(codegen, "} while(TRUE);\n");
    /* End code */
    Codegen_indent(codegen, -1);
    Codegen_printf(codegen, "end_of_dp:\n");
    Codegen_indent(codegen, 1);
    /* End of viterbi */
    /**/
    CGUtil_prep(codegen, viterbi->model, FALSE);
    Viterbi_implement_finalise(viterbi, codegen);
    /* Copy final cell */
    if(viterbi->use_continuation){
        Codegen_printf(codegen,
            "src = prev_row_0[region->query_length]"
            "[vd->continuation->final_state->id];\n");
        for(i = 0; i < cell_size; i++)
            Codegen_printf(codegen,
                "vd->continuation->final_cell[%d] = src[%d];\n",
                i, i);
        }
    if(viterbi->mode == Viterbi_Mode_FIND_CHECKPOINTS){
        g_assert(viterbi->use_continuation);
        Codegen_printf(codegen,
            "vd->checkpoint->last_srp = prev_row_0"
                      "[region->query_length]"
                      "[vd->continuation->final_state->id]"
                      "[vd->vr->cell_size-1];\n");
        }
    Codegen_printf(codegen, "g_free(state_is_set);\n");
    Codegen_printf(codegen, "return score;\n");
    Codegen_printf(codegen, "}\n\n");
    Codegen_indent(codegen, -1);
    return;
    }

Codegen *Viterbi_make_Codegen(Viterbi *viterbi){
    register Codegen *codegen = Codegen_create(NULL, viterbi->name);
    CGUtil_print_header(codegen, viterbi->model);
    Codegen_printf(codegen,
        "/* Exhaustive Viterbi DP implementation. */\n"
        "\n"
        "#include \"viterbi.h\"\n"
        "\n");
    Codegen_printf(codegen, "/* Viterbi model:\n"
                            " * [%s]\n"
                            " *\n"
                            " * max query advance: %d\n"
                            " * max target advance: %d\n"
                            " *\n"
                            " * scope [%s]->[%s]\n"
                            " */\n\n",
                  viterbi->name,
                  viterbi->model->max_query_advance,
                  viterbi->model->max_target_advance,
                  C4_Scope_get_name(viterbi->model->start_state->scope),
                  C4_Scope_get_name(viterbi->model->end_state->scope));
    Viterbi_implement_dp(viterbi, codegen);
    CGUtil_print_footer(codegen);
    CGUtil_compile(codegen, viterbi->model);
    return codegen;
    }
/* FIXME: optimisation : could reuse the Layout
 *        for different DP modes (although will still need separate
 *        global Layout for continuation DP).
 */

