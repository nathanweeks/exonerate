/****************************************************************\
*                                                                *
*  C4 dynamic programming library - Boundary Object              *
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

#include "boundary.h"

/**/

static Boundary_Interval *Boundary_Interval_create(gint query_pos,
                                                   gint seed_id){
    register Boundary_Interval *interval = g_new(Boundary_Interval, 1);
    interval->query_pos = query_pos;
    interval->length = 1;
    interval->seed_id = seed_id;
    return interval;
    }

static void Boundary_Interval_destroy(Boundary_Interval *interval){
    g_free(interval);
    return;
    }

static Boundary_Interval *Boundary_Interval_copy(
                          Boundary_Interval *interval){
    register Boundary_Interval *new_interval
           = Boundary_Interval_create(interval->query_pos,
                                      interval->seed_id);
    new_interval->length = interval->length;
    return new_interval;
    }

#define Boundary_Interval_end(interval)            \
                            ((interval)->query_pos \
                           + (interval)->length)

/**/

static Boundary_Row *Boundary_Row_create(gint target_pos){
    register Boundary_Row *boundary_row = g_new(Boundary_Row, 1);
    boundary_row->target_pos = target_pos;
    boundary_row->interval_list = g_ptr_array_new();
    return boundary_row;
    }

static void Boundary_Row_destroy(Boundary_Row *boundary_row){
    register gint i;
    register Boundary_Interval *interval;
    for(i = 0; i < boundary_row->interval_list->len; i++){
        interval = boundary_row->interval_list->pdata[i];
        Boundary_Interval_destroy(interval);
        }
    g_ptr_array_free(boundary_row->interval_list, TRUE);
    g_free(boundary_row);
    return;
    }

static void Boundary_Row_print_gnuplot(Boundary_Row *boundary_row,
                                       gint colour){
    register gint i;
    register Boundary_Interval *interval;
    for(i = 0; i < boundary_row->interval_list->len; i++){
        interval = boundary_row->interval_list->pdata[i];
        g_print("set arrow from %d,%d to %d,%d heads ls %d\n",
                interval->query_pos,
                boundary_row->target_pos,
                interval->query_pos+interval->length,
                boundary_row->target_pos, colour);
        }
    return;
    }

static Boundary_Row *Boundary_Row_copy(Boundary_Row *boundary_row){
    register Boundary_Row *new_row = Boundary_Row_create(
                                           boundary_row->target_pos);
    register gint i;
    register Boundary_Interval *interval;
    g_assert(boundary_row->interval_list->len);
    for(i = 0; i < boundary_row->interval_list->len; i++){
        interval = boundary_row->interval_list->pdata[i];
        g_ptr_array_add(new_row->interval_list,
                        Boundary_Interval_copy(interval));
        }
    return new_row;
    }

static void Boundary_Row_reverse(Boundary_Row *boundary_row){
    register gint a, z;
    register Boundary_Interval *swap;
    for(a = 0, z = boundary_row->interval_list->len-1; a < z; a++, z--){
        swap = boundary_row->interval_list->pdata[a];
        boundary_row->interval_list->pdata[a]
            = boundary_row->interval_list->pdata[z];
        boundary_row->interval_list->pdata[z] = swap;
        }
    return;
    }

static void Boundary_Row_add_interval(Boundary_Row *boundary_row,
                          gint query_pos, gint length, gint seed_id){
    register Boundary_Interval *interval
           = Boundary_Interval_create(query_pos, seed_id);
    interval->length = length;
    g_ptr_array_add(boundary_row->interval_list, interval);
    return;
    }

static void Boundary_Row_select(Boundary_Row *boundary_row,
                                Boundary_Row *sub_boundary_row,
                                Region *region){
    register gint i, query_start, query_end;
    register Boundary_Interval *interval;
    for(i = 0; i < boundary_row->interval_list->len; i++){
        interval = boundary_row->interval_list->pdata[i];
        if((interval->query_pos+interval->length) < region->query_start)
            continue;
        if(interval->query_pos > Region_query_end(region))
            break;
        query_start = MAX(interval->query_pos, region->query_start);
        query_end = MIN(interval->query_pos+interval->length,
                        Region_query_end(region));
        if(query_end-query_start > 0)
            Boundary_Row_add_interval(sub_boundary_row,
                                      query_start, query_end-query_start,
                                      interval->seed_id);
        }
    return;
    }

static Boundary_Interval *Boundary_Row_get_last_interval(
                          Boundary_Row *boundary_row){
    register Boundary_Interval *interval;
    g_assert(boundary_row);
    g_assert(boundary_row->interval_list);
    if(!boundary_row->interval_list->len)
        return NULL;
    interval = boundary_row->interval_list->pdata
              [boundary_row->interval_list->len-1];
    g_assert(interval);
    return interval;
    }

void Boundary_Row_prepend(Boundary_Row *boundary_row,
                          gint query_pos, gint seed_id){
    register Boundary_Interval *interval
           = Boundary_Row_get_last_interval(boundary_row);
    if(interval
    && (interval->seed_id == seed_id)
    && ((interval->query_pos - 1) == query_pos)){
        interval->query_pos = query_pos;
        interval->length++;
        return;
        }
    Boundary_Row_add_interval(boundary_row, query_pos, 1, seed_id);
    return;
    }

static void Boundary_Row_append_interval(Boundary_Row *boundary_row,
                                         Boundary_Interval *interval){
    register Boundary_Interval *last_interval
           = Boundary_Row_get_last_interval(boundary_row);
    g_assert(interval->query_pos >= 0);
    g_assert(interval->seed_id >= 0);
    g_assert(interval->length > 0);
    if(last_interval
    && (last_interval->seed_id == interval->seed_id)
    && ((Boundary_Interval_end(last_interval) == interval->query_pos))){
        last_interval->length += interval->length;
        return;
        }
    Boundary_Row_add_interval(boundary_row, interval->query_pos,
                                            interval->length,
                                            interval->seed_id);
    return;
    }

typedef enum {
    Boundary_State_NULL     = 0,
    Boundary_State_START    = (1<<0),
    Boundary_State_END      = (1<<1),
    Boundary_State_BOUNDARY = (1<<2),
    Boundary_State_INSERT   = (1<<3)
} Boundary_State;

#define Boundary_State_pair(prev_state, curr_state) \
                           (((prev_state) << 4) | (curr_state))

static void Boundary_Row_insert_append(Boundary_Row *boundary_row,
            Boundary_Interval *interval,
            Boundary_Interval *boundary_stored,
            Boundary_Interval *insert_stored,
            Boundary_State *prev_state, Boundary_State *curr_state){
    Boundary_Interval temp_interval;
    switch(Boundary_State_pair(*prev_state, *curr_state)){
        case Boundary_State_pair(Boundary_State_START,
                                 Boundary_State_INSERT):
            /* S->I: report I, store I */
            Boundary_Row_append_interval(boundary_row, interval);
            insert_stored->query_pos = interval->query_pos;
            insert_stored->length    = interval->length;
            insert_stored->seed_id   = interval->seed_id;
            break;
        case Boundary_State_pair(Boundary_State_START,
                                 Boundary_State_BOUNDARY):
            /* S->B: store B */
            boundary_stored->query_pos = interval->query_pos;
            boundary_stored->length    = interval->length;
            boundary_stored->seed_id   = interval->seed_id;
            break;
        case Boundary_State_pair(Boundary_State_INSERT,
                                 Boundary_State_BOUNDARY):
        case Boundary_State_pair(Boundary_State_BOUNDARY,
                                 Boundary_State_BOUNDARY):
            /* ON B: report any stored_B
             *       truncate B (after last I)
             *       store B if +ve len else clear
             */
            if(boundary_stored->seed_id != -1){
                Boundary_Row_append_interval(boundary_row,
                                             boundary_stored);
                }
            boundary_stored->query_pos = interval->query_pos;
            boundary_stored->length    = interval->length;
            boundary_stored->seed_id   = interval->seed_id;
            if(insert_stored->seed_id != -1){
                if(interval->query_pos
                 < Boundary_Interval_end(insert_stored)){
                    boundary_stored->length
                                = (Boundary_Interval_end(interval)
                                - Boundary_Interval_end(insert_stored));
                    boundary_stored->query_pos
                               = Boundary_Interval_end(insert_stored);
                    boundary_stored->seed_id = interval->seed_id;
                    if(boundary_stored->length <= 0){
                        boundary_stored->query_pos = -1;
                        boundary_stored->length    = -1;
                        boundary_stored->seed_id   = -1;
                        }
                    }
                }
            break;
        case Boundary_State_pair(Boundary_State_INSERT,
                                 Boundary_State_INSERT):
        case Boundary_State_pair(Boundary_State_BOUNDARY,
                                 Boundary_State_INSERT):
            /* ON I: report any stored_B before I
             *       store any remaining stored B after I
             *       report I
             *       store I
             */
            /**/
            if(boundary_stored->seed_id != -1){
                if(interval->query_pos
                  < Boundary_Interval_end(boundary_stored)){
                        temp_interval.query_pos
                                         = boundary_stored->query_pos;
                        temp_interval.length = interval->query_pos
                                         - boundary_stored->query_pos;
                        temp_interval.seed_id
                                       = boundary_stored->seed_id;
                        if(temp_interval.length > 0)
                            Boundary_Row_append_interval(boundary_row,
                                                       &temp_interval);
                    } else {
                        Boundary_Row_append_interval(boundary_row,
                                                     boundary_stored);
                        }
                /**/
                if(Boundary_Interval_end(interval)
                 < Boundary_Interval_end(boundary_stored)){
                    boundary_stored->length
                        = Boundary_Interval_end(boundary_stored)
                        - Boundary_Interval_end(interval);
                    boundary_stored->query_pos
                        = Boundary_Interval_end(interval);
                    /* (keeps seed_id) */
                } else {
                    boundary_stored->query_pos = -1;
                    boundary_stored->length    = -1;
                    boundary_stored->seed_id   = -1;
                    }
                }
            /**/
            Boundary_Row_append_interval(boundary_row, interval);
            insert_stored->query_pos = interval->query_pos;
            insert_stored->length    = interval->length;
            insert_stored->seed_id   = interval->seed_id;
            break;
        case Boundary_State_pair(Boundary_State_INSERT,
                                 Boundary_State_END):
            /*fallthrough*/
        case Boundary_State_pair(Boundary_State_BOUNDARY,
                                 Boundary_State_END):
            /* [IB]->E: report stored B */
            if(boundary_stored->seed_id != -1)
                Boundary_Row_append_interval(boundary_row,
                                             boundary_stored);
            break;
        default:
            g_error("Bad boundary state pair [%d,%d]",
                   *prev_state, *curr_state);
            break;
        }
    return;
    }

static gboolean Boundary_Row_is_valid(Boundary_Row *boundary_row){
    register gint i;
    register Boundary_Interval *interval, *prev;
    g_assert(boundary_row->interval_list->len);
    prev = boundary_row->interval_list->pdata[0];
    g_assert(prev->length > 0);
    g_assert(prev->seed_id >= 0);
    for(i = 1; i < boundary_row->interval_list->len; i++){
        interval = boundary_row->interval_list->pdata[i];
        g_assert(interval->length > 0);
        g_assert(interval->seed_id >= 0);
        g_assert(prev->query_pos < interval->query_pos);
        g_assert(Boundary_Interval_end(prev) <= interval->query_pos);
        prev = interval;
        }
    return TRUE;
    }

static void Boundary_Row_insert(Boundary_Row *boundary_row,
                                Boundary_Row *insert_row){
    register gint i = 0, j = 0;
    register Boundary_Interval *boundary_interval, *insert_interval;
    register Boundary_Row *combined_row
           = Boundary_Row_create(boundary_row->target_pos);
    Boundary_Interval boundary_stored = {-1, -1, -1},
                      insert_stored   = {-1, -1, -1};
    Boundary_State curr_state = Boundary_State_START,
                   prev_state = Boundary_State_NULL;
    g_assert(boundary_row->interval_list->len);
    g_assert(insert_row->interval_list->len);
    g_assert(boundary_row->target_pos == insert_row->target_pos);
    do {
        if(i == insert_row->interval_list->len)
            break;
        if(j == boundary_row->interval_list->len)
            break;
        insert_interval = insert_row->interval_list->pdata[i];
        boundary_interval = boundary_row->interval_list->pdata[j];
        if(insert_interval->query_pos < boundary_interval->query_pos){
            prev_state = curr_state;
            curr_state = Boundary_State_INSERT;
            Boundary_Row_insert_append(combined_row, insert_interval,
                                &boundary_stored, &insert_stored,
                                &prev_state, &curr_state);
            i++;
        } else {
            prev_state = curr_state;
            curr_state = Boundary_State_BOUNDARY;
            Boundary_Row_insert_append(combined_row, boundary_interval,
                                &boundary_stored, &insert_stored,
                                &prev_state, &curr_state);
            j++;
            }
    } while(TRUE);
    while(i < insert_row->interval_list->len){
        prev_state = curr_state;
        curr_state = Boundary_State_INSERT;
        insert_interval = insert_row->interval_list->pdata[i++];
        Boundary_Row_insert_append(combined_row, insert_interval,
                                   &boundary_stored, &insert_stored,
                                   &prev_state, &curr_state);
        }
    while(j < boundary_row->interval_list->len){
        prev_state = curr_state;
        curr_state = Boundary_State_BOUNDARY;
        boundary_interval = boundary_row->interval_list->pdata[j++];
        Boundary_Row_insert_append(combined_row, boundary_interval,
                                   &boundary_stored, &insert_stored,
                                   &prev_state, &curr_state);
        }
    prev_state = curr_state;
    curr_state = Boundary_State_END;
    /* B -> E  report stored */
    Boundary_Row_insert_append(combined_row, NULL,
                               &boundary_stored, &insert_stored,
                               &prev_state, &curr_state);
    for(i = 0; i < boundary_row->interval_list->len; i++){
        boundary_interval = boundary_row->interval_list->pdata[i];
        Boundary_Interval_destroy(boundary_interval);
        }
    g_ptr_array_free(boundary_row->interval_list, TRUE);
    boundary_row->interval_list = combined_row->interval_list;
    g_free(combined_row);
    g_assert(boundary_row->interval_list->len);
    g_assert(Boundary_Row_is_valid(boundary_row));
    return;
    }

/**/

Boundary *Boundary_create(void){
    register Boundary *boundary = g_new(Boundary, 1);
    boundary->row_list = g_ptr_array_new();
    boundary->ref_count = 1;
    return boundary;
    }

void Boundary_destroy(Boundary *boundary){
    register gint i;
    register Boundary_Row *boundary_row;
    g_assert(boundary);
    if(--boundary->ref_count)
        return;
    for(i = 0; i < boundary->row_list->len; i++){
        boundary_row = boundary->row_list->pdata[i];
        Boundary_Row_destroy(boundary_row);
        }
    g_ptr_array_free(boundary->row_list, TRUE);
    g_free(boundary);
    return;
    }

Boundary *Boundary_share(Boundary *boundary){
    g_assert(boundary);
    boundary->ref_count++;
    return boundary;
    }

/**/

Boundary_Row *Boundary_add_row(Boundary *boundary, gint target_pos){
    register Boundary_Row *boundary_row;
    g_assert(boundary);
    g_assert(boundary->row_list);
    boundary_row = Boundary_Row_create(target_pos);
    g_ptr_array_add(boundary->row_list, boundary_row);
    return boundary_row;
    }

Boundary_Row *Boundary_get_last_row(Boundary *boundary){
    register Boundary_Row *boundary_row;
    g_assert(boundary);
    g_assert(boundary->row_list);
    if(!boundary->row_list->len)
        return NULL;
    boundary_row = boundary->row_list->pdata[boundary->row_list->len-1];
    g_assert(boundary_row);
    return boundary_row;
    }

void Boundary_remove_empty_last_row(Boundary *boundary){
    register Boundary_Row *boundary_row
           = Boundary_get_last_row(boundary);
    if(!boundary_row)
        return;
    if(boundary_row->interval_list->len)
        return; /* Don't remove if not empty */
    Boundary_Row_destroy(boundary_row);
    g_ptr_array_set_size(boundary->row_list, boundary->row_list->len-1);
    return;
    }

/* FIXME: optimisation: use RecycleBin for boundary rows and cells */


/**/

void Boundary_reverse(Boundary *boundary){
    register gint a, z;
    register Boundary_Row *swap;
    for(a = 0, z = boundary->row_list->len-1; a < z; a++, z--){
        swap = boundary->row_list->pdata[a];
        boundary->row_list->pdata[a] = boundary->row_list->pdata[z];
        boundary->row_list->pdata[z] = swap;
        Boundary_Row_reverse(boundary->row_list->pdata[a]);
        Boundary_Row_reverse(boundary->row_list->pdata[z]);
        }
    if(boundary->row_list->len & 1){ /* Reverse central row */
        a = boundary->row_list->len >> 1;
        Boundary_Row_reverse(boundary->row_list->pdata[a]);
        }
    return;
    }
/* Reverse order of the rows and intervals */

Boundary *Boundary_select(Boundary *boundary, Region *region){
    register Boundary *sub_boundary = Boundary_create();
    register Boundary_Row *boundary_row, *sub_boundary_row;
    register gint i;
    g_assert(boundary);
    for(i = 0; i < boundary->row_list->len; i++){
        boundary_row = boundary->row_list->pdata[i];
        if(boundary_row->target_pos < region->target_start)
            continue;
        if(boundary_row->target_pos > Region_target_end(region))
            break;
        sub_boundary_row = Boundary_add_row(sub_boundary,
                                            boundary_row->target_pos);
        Boundary_Row_select(boundary_row, sub_boundary_row, region);
        Boundary_remove_empty_last_row(sub_boundary);
        }
    return sub_boundary;
    }
/* Returns only the boundary within region
 */

static gboolean Boundary_is_valid(Boundary *boundary){
    register gint i;
    register Boundary_Row *prev_row, *row;
    g_assert(boundary);
    g_assert(boundary->row_list->len);
    prev_row = boundary->row_list->pdata[0];
    g_assert(prev_row);
    g_assert(Boundary_Row_is_valid(prev_row));
    for(i = 1; i < boundary->row_list->len; i++){
        row = boundary->row_list->pdata[i];
        g_assert(row);
        g_assert(prev_row->target_pos < row->target_pos);
        g_assert(Boundary_Row_is_valid(row));
        prev_row = row;
        }
    return TRUE;
    }

void Boundary_insert(Boundary *boundary, Boundary *insert){
    register gint i = 0, j = 0;
    register Boundary_Row *boundary_row, *insert_row,
                          *last_boundary_row;
    register GPtrArray *combined_row_list = g_ptr_array_new();
    do {
        if(i == insert->row_list->len)
            break;
        if(j == boundary->row_list->len)
            break;
        insert_row = insert->row_list->pdata[i];
        boundary_row = boundary->row_list->pdata[j];
        if(insert_row->target_pos < boundary_row->target_pos){
            i++;
            if(j){
                last_boundary_row = boundary->row_list->pdata[j-1];
                if(insert_row->target_pos
                == last_boundary_row->target_pos){
                    Boundary_Row_insert(last_boundary_row, insert_row);
                } else {
                    g_ptr_array_add(combined_row_list,
                                    Boundary_Row_copy(insert_row));
                    }
            } else {
                g_ptr_array_add(combined_row_list,
                                Boundary_Row_copy(insert_row));
                }
        } else {
            j++;
            g_assert(boundary_row->interval_list->len);
            g_ptr_array_add(combined_row_list, boundary_row);
            }
    } while(TRUE);
    /* Handle repeated last row */
    if((i < insert->row_list->len) && j){
        last_boundary_row = boundary->row_list->pdata[j-1];
        insert_row = insert->row_list->pdata[i];
        if(insert_row->target_pos == last_boundary_row->target_pos){
            Boundary_Row_insert(last_boundary_row, insert_row);
            i++;
            }
        }
    while(i < insert->row_list->len){
        insert_row = insert->row_list->pdata[i++];
        g_ptr_array_add(combined_row_list,
                        Boundary_Row_copy(insert_row));
        }
    while(j < boundary->row_list->len){
        boundary_row = boundary->row_list->pdata[j++];
        g_assert(boundary_row->interval_list->len);
        g_ptr_array_add(combined_row_list, boundary_row);
        }
    g_ptr_array_free(boundary->row_list, TRUE);
    boundary->row_list = combined_row_list;
    g_assert(Boundary_is_valid(boundary));
    return;
    }

/**/

static Region *Boundary_find_bounding_region(Boundary *boundary){
    register Region *region = Region_create_blank();
    register gint i, qmin, qmax;
    register Boundary_Row *boundary_row;
    register Boundary_Interval *interval;
    g_assert(boundary->row_list->len);
    /* take first row */
    boundary_row = boundary->row_list->pdata[0];
    region->target_start = boundary_row->target_pos;
    /* take 1st interval */
    interval = boundary_row->interval_list->pdata[0];
    qmin = interval->query_pos;
    /* take last interval */
    interval = boundary_row->interval_list->pdata
              [boundary_row->interval_list->len-1];
    qmax = interval->query_pos + interval->length;
    /* take last row */
    boundary_row = boundary->row_list->pdata[boundary->row_list->len-1];
    region->target_length = boundary_row->target_pos;
    for(i = 1; i < boundary->row_list->len; i++){
        boundary_row = boundary->row_list->pdata[i];
        g_assert(boundary_row->interval_list->len);
        /* take first interval */
        interval = boundary_row->interval_list->pdata[0];
        if(qmin > interval->query_pos)
            qmin = interval->query_pos;
        /* take last interval */
        interval = boundary_row->interval_list->pdata
                  [boundary_row->interval_list->len-1];
        if(qmax < (interval->query_pos + interval->length))
            qmax = interval->query_pos + interval->length;
        }
    region->query_start = qmin;
    region->query_length = qmax - qmin;
    return region;
    }

void Boundary_print_gnuplot(Boundary *boundary, gint colour){
    register gint i;
    register Boundary_Row *boundary_row;
    /* Find bounding region */
    register Region *region = Boundary_find_bounding_region(boundary);
    g_print("# Begin gnuplot commands\n");
    g_print("set xrange [%d:%d]\n",
           region->query_start, Region_query_end(region));
    g_print("set yrange [%d:%d]\n",
           region->target_start, Region_target_end(region));
    g_print("set style line %d\n", colour);
    g_print("set grid\n");
    g_print("unset key\n");
    /**/
    for(i = 0; i < boundary->row_list->len; i++){
        boundary_row = boundary->row_list->pdata[i];
        Boundary_Row_print_gnuplot(boundary_row, colour);
        }
    g_print("plot 0\n");
    g_print("# End gnuplot commands\n");
    Region_destroy(region);
    return;
    }

/**/

