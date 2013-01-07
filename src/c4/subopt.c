/****************************************************************\
*                                                                *
*  C4 dynamic programming library - suboptimal alignments        *
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

#include "subopt.h"

SubOpt *SubOpt_create(gint query_length, gint target_length){
    register SubOpt *subopt = g_new(SubOpt, 1);
    subopt->ref_count = 1;
    subopt->query_length = query_length;
    subopt->target_length = target_length;
    subopt->range_tree = RangeTree_create();
    subopt->path_count = 0;
    return subopt;
    }

SubOpt *SubOpt_share(SubOpt *subopt){
    subopt->ref_count++;
    return subopt;
    }

void SubOpt_destroy(SubOpt *subopt){
    if(--subopt->ref_count)
        return;
    RangeTree_destroy(subopt->range_tree, NULL, NULL);
    g_free(subopt);
    return;
    }

/* Greatest Common Divisor : Euclid's algorithm */
static gint SubOpt_get_gcd(gint a, gint b){
    register gint t;
    g_assert(a > 0);
    g_assert(b > 0);
    while(a > 0){
        if(a < b){
            t = a;
            a = b;
            b = t;
            }
        a -= b;
        }
    return b;
    }

static gboolean SubOpt_check_pos(SubOpt *subopt,
                                 gint query_pos, gint target_pos){
    return RangeTree_check_pos(subopt->range_tree,
                               query_pos, target_pos);
    }

static void SubOpt_add_AlignmentOperation(SubOpt *subopt,
              AlignmentOperation *ao, gint query_pos, gint target_pos){
    register gint i;
    register gint gcd = SubOpt_get_gcd(ao->transition->advance_query,
                                       ao->transition->advance_target);
    register gint q_move = ao->transition->advance_query / gcd,
                  t_move = ao->transition->advance_target / gcd;
    register gint q_limit = query_pos, t_limit = target_pos;
    register gint qp, tp;
    /* Add operation */
    for(i = 0; i < ao->length; i++){
        /* Block positions in transition */
        qp = q_limit;
        tp = t_limit;
        q_limit += ao->transition->advance_query;
        t_limit += ao->transition->advance_target;
        while(qp < q_limit){
            /* We need to ensure that the position is unused,
             * otherwise lead-out positions could clash with lead-in
             * positions from a previous alignment.
             */
            if(!SubOpt_check_pos(subopt,
                     /* qp + ao->transition->advance_query,
                        tp + ao->transition->advance_target)){ */
                        qp, tp)){
                RangeTree_add(subopt->range_tree,
                    /* qp + ao->transition->advance_query,
                    tp + ao->transition->advance_target, */
                    qp, tp,
                    GINT_TO_POINTER(subopt->path_count));
                }
            qp += q_move;
            tp += t_move;
            }
        }
    /* Block leading in positions */
    qp = query_pos - ao->transition->advance_query + q_move;
    tp = target_pos - ao->transition->advance_target + t_move;
    while(qp < query_pos){
        /* Only insert if not already present */
        if(!SubOpt_check_pos(subopt,
                    /*
                         qp + ao->transition->advance_query,
                         tp + ao->transition->advance_target)){
                     */
                             qp, tp)){
            /* FIXME: */
            if((qp >= 0) && (tp >= 0))
            RangeTree_add(subopt->range_tree,
                    /* qp + ao->transition->advance_query,
                    tp + ao->transition->advance_target, */
                    qp, tp,
                    GINT_TO_POINTER(subopt->path_count));
            }
        qp += q_move;
        tp += t_move;
        }
    return;
    }
/* FIXME: tidy */

void SubOpt_add_alignment(SubOpt *subopt, Alignment *alignment){
    register gint i;
    register gint query_pos, target_pos;
    register AlignmentOperation *ao;
    query_pos = alignment->region->query_start;
    target_pos = alignment->region->target_start;
    for(i = 0; i < alignment->operation_list->len; i++){
        ao = alignment->operation_list->pdata[i];
        if(C4_Transition_is_match(ao->transition)){
            SubOpt_add_AlignmentOperation(subopt, ao,
                                          query_pos, target_pos);
            }
        query_pos += (ao->transition->advance_query * ao->length);
        target_pos += (ao->transition->advance_target * ao->length);
        }
    subopt->path_count++;
    return;
    }

/**/

typedef struct {
    gpointer user_data;
    SubOpt_FindFunc find_func;
} SubOpt_FindData;

static gboolean SubOpt_RangeTree_find(gint x, gint y,
                                 gpointer info, gpointer user_data){
    register SubOpt_FindData *sufd = user_data;
    register gint path_id = GPOINTER_TO_INT(info);
    if(sufd->find_func(x, y, path_id, sufd->user_data))
        return TRUE;
    return FALSE;
    }

static gboolean SubOpt_overlap_find_func(gint query_pos, gint target_pos,
                                         gint path_id, gpointer user_data){
    return TRUE;
    }

gboolean SubOpt_find(SubOpt *subopt, Region *region,
                     SubOpt_FindFunc find_func, gpointer user_data){
    SubOpt_FindData sufd;
    sufd.user_data = user_data;
    sufd.find_func = find_func;
    return RangeTree_find(subopt->range_tree,
            region->query_start, region->query_length,
            region->target_start, region->target_length,
            SubOpt_RangeTree_find, &sufd);
    }

gboolean SubOpt_overlaps_alignment(SubOpt *subopt, Alignment *alignment){
    register gint i, j;
    register AlignmentOperation *ao;
    Region region;
    register gint qp = alignment->region->query_start,
                  tp = alignment->region->target_start;
    for(i = 0; i < alignment->operation_list->len; i++){
        ao = alignment->operation_list->pdata[i];
        if(C4_Transition_is_match(ao->transition)){
            for(j = 0; j < ao->length; j++){
                region.query_start = qp;
                region.target_start = tp;
                region.query_length = ao->transition->advance_query;
                region.target_length = ao->transition->advance_target;
                if(SubOpt_find(subopt, &region,
                               SubOpt_overlap_find_func, NULL))
                    return TRUE;
                qp += ao->transition->advance_query;
                tp += ao->transition->advance_target;
                }
        } else {
            qp += (ao->transition->advance_query * ao->length);
            tp += (ao->transition->advance_target * ao->length);
            }
        }
    return FALSE;
    }

/**/

static SubOpt_Index_Row *SubOpt_Index_Row_create(gint target_pos){
    register SubOpt_Index_Row *soir = g_new(SubOpt_Index_Row, 1);
    soir->target_pos = target_pos;
    soir->total = 0;
    soir->query_pos = NULL;
    return soir;
    }

static void SubOpt_Index_Row_destroy(SubOpt_Index_Row *soir){
    if(soir->query_pos)
        g_free(soir->query_pos);
    g_free(soir);
    return;
    }

/**/

typedef struct {
    gint query_pos;
    gint target_pos;
} SubOpt_Point;

static gboolean SubOpt_RangeTree_traverse(gint x, gint y, gpointer info,
                                          gpointer user_data){
    register GPtrArray *point_list = user_data;
    register SubOpt_Point *sop = g_new(SubOpt_Point, 1);
    sop->query_pos = x;
    sop->target_pos = y;
    g_ptr_array_add(point_list, sop);
    return FALSE;
    }

static int SubOpt_sort_by_target_then_query_pos(const void *a,
                                                const void *b){
    register SubOpt_Point **point_a = (SubOpt_Point**)a,
                          **point_b = (SubOpt_Point**)b;
    register gint target_diff = (*point_a)->target_pos
                              - (*point_b)->target_pos;
    if(!target_diff)
        return (*point_a)->query_pos - (*point_b)->query_pos;
    return target_diff;
    }

SubOpt_Index *SubOpt_Index_create(SubOpt *subopt, Region *region){
    register SubOpt_Index *soi;
    register GPtrArray *point_list = g_ptr_array_new();
    register gint i, j, start;
    register SubOpt_Point *point, *prev_point = NULL;
    register SubOpt_Index_Row *soir = NULL;
    /* Copy all points in region into point_list */
    RangeTree_find(subopt->range_tree,
                   region->query_start, region->query_length+1,
                   region->target_start, region->target_length+1,
                   SubOpt_RangeTree_traverse, point_list);
    if(!point_list->len){
        g_ptr_array_free(point_list, TRUE);
        return NULL; /* Found no points in region */
        }
    /* Convert all points to region coordinates */
    for(i = 0; i < point_list->len; i++){
        point = point_list->pdata[i];
        point->query_pos -= region->query_start;
        point->target_pos -= region->target_start;
        }
    /* Sort points to correct order for DP */
    qsort(point_list->pdata, point_list->len,
          sizeof(gpointer), SubOpt_sort_by_target_then_query_pos);
    soi = g_new(SubOpt_Index, 1);
    soi->region = Region_share(region);
    soi->row_list = g_ptr_array_new();
    soi->curr_row_index = 0;
    soi->curr_query_index = 0;
    /* Set the row sizes */
    for(i = 0; i < point_list->len; i++){
        point = point_list->pdata[i];
        if((!soir) || (soir->target_pos != point->target_pos)){
            soir = SubOpt_Index_Row_create(point->target_pos);
            soir->total = 1;
            g_ptr_array_add(soi->row_list, soir);
            soir->query_pos = GINT_TO_POINTER(i); /* << Hack */
        } else {
            g_assert(soir->target_pos == point->target_pos);
            soir->total++;
            }
        }
    /* Fill the rows */
    for(i = 0; i < soi->row_list->len; i++){
        soir = soi->row_list->pdata[i];
        g_assert(soir->total);
        start = GPOINTER_TO_INT(soir->query_pos); /* >> Hack */
        soir->query_pos = g_new(gint, soir->total+1);
        prev_point = NULL;
        for(j = 0; j < soir->total; j++){
            point = point_list->pdata[start+j];
            soir->query_pos[j] = point->query_pos;
            g_assert(soir->target_pos == point->target_pos);
            g_assert((!prev_point) /* Ensure points are unique */
                  || (prev_point->query_pos != point->query_pos));
            prev_point = point;
            }
        /* Tag on dummy end position (query_length+1) */
        soir->query_pos[soir->total] = subopt->query_length + 1;
        }
    /**/
    for(i = 0; i < point_list->len; i++){
        point = point_list->pdata[i];
        g_free(point);
        }
    g_ptr_array_free(point_list, TRUE);
    soi->curr_row = soi->blank_row;
    soi->blank_row = SubOpt_Index_Row_create(subopt->target_length+1);
    soi->blank_row->total = 1;
    soi->blank_row->query_pos = g_new(gint, 1);
    soi->blank_row->query_pos[0] = subopt->query_length+1;
    g_ptr_array_add(soi->row_list, soi->blank_row);
    return soi;
    }

void SubOpt_Index_destroy(SubOpt_Index *soi){
    register gint i;
    register SubOpt_Index_Row *soir;
    /* Free all rows (including final blank row) */
    for(i = 0; i < soi->row_list->len; i++){
        soir = soi->row_list->pdata[i];
        SubOpt_Index_Row_destroy(soir);
        }
    g_ptr_array_free(soi->row_list, TRUE);
    Region_destroy(soi->region);
    g_free(soi);
    return;
    }

void SubOpt_Index_set_row(SubOpt_Index *soi, gint target_pos){
    register SubOpt_Index_Row *soir = NULL;
    g_assert(target_pos >= 0);
    if(!soi)
        return;
    /**/
    soir = soi->row_list->pdata[soi->curr_row_index];
    while(soir->target_pos < target_pos){
        if(soi->curr_row_index < soi->row_list->len){
            soi->curr_row_index++;
            soir = soi->row_list->pdata[soi->curr_row_index];
        } else {
            soir = NULL;
            break;
            }
        }
    if(soir){
        while(soir->target_pos > target_pos){
            if(soi->curr_row_index > 0){
                soi->curr_row_index--;
                soir = soi->row_list->pdata[soi->curr_row_index];
            } else {
                soir = NULL;
                break;
                }
            }
        if(soir && (soir->target_pos == target_pos)){
            soi->curr_row = soir;
        } else {
            soi->curr_row = soi->blank_row;
            }
        }
    /**/
    soi->curr_query_index = 0;
    return;
    }

gboolean SubOpt_Index_is_blocked(SubOpt_Index *soi, gint query_pos){
    g_assert(query_pos >= 0);
    while(soi->curr_row->query_pos[soi->curr_query_index] <
            query_pos){
        soi->curr_query_index++;
        }
    while(soi->curr_row->query_pos[soi->curr_query_index] > query_pos){
        if(!soi->curr_query_index)
            break;
        soi->curr_query_index--;
        }
    return (soi->curr_row->query_pos[soi->curr_query_index]
            == query_pos);
    }
/* FIXME: optimisation: use fwd and rev only versions of this for sdp */


/**/

