/****************************************************************\
*                                                                *
*  BSDP : Bounded Sparse Dynamic Programming                     *
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

#include "bsdp.h"

/**/

BSDP_ArgumentSet *BSDP_ArgumentSet_create(Argument *arg){
    register ArgumentSet *as;
    static BSDP_ArgumentSet bas = {0};
    if(arg){
        as = ArgumentSet_create("BSDP algorithm options");
        ArgumentSet_add_option(as, '\0', "joinfilter", NULL,
            "BSDP join filter threshold", "0",
            Argument_parse_int, &bas.join_filter);
        Argument_absorb_ArgumentSet(arg, as);
        }
    return &bas;
    }

/**/

static BSDP_Edge *BSDP_Edge_create(BSDP *bsdp, gpointer edge_data,
             BSDP_Node *dst, C4_Score bound_score){
    register BSDP_Edge *bsdp_edge
        = RecycleBin_alloc(bsdp->edge_recycle);
    g_assert(dst);
    bsdp_edge->edge_data = edge_data;
    bsdp_edge->dst = dst;
    bsdp_edge->join_score = bound_score;
    bsdp_edge->stored_partial = 0;
    bsdp_edge->mailbox = -1;
    return bsdp_edge;
    }

static void BSDP_Edge_destroy(BSDP *bsdp, BSDP_Edge *bsdp_edge){
    if(bsdp->destroy_edge_data_func)
        bsdp->destroy_edge_data_func(bsdp_edge->edge_data);
    RecycleBin_recycle(bsdp->edge_recycle, bsdp_edge);
    return;
    }

/**/

static BSDP_Node *BSDP_Node_create(gpointer node_data,
                                   C4_Score node_score,
                                   gboolean is_valid_start,
                                   gboolean is_valid_end,
                                   C4_Score start_bound,
                                   C4_Score end_bound){
    register BSDP_Node *bsdp_node = g_new(BSDP_Node, 1);
    bsdp_node->mask = BSDP_Node_Mask_IS_NEW;
    if(is_valid_start){
        bsdp_node->mask |= BSDP_Node_Mask_IS_VALID_START;
        bsdp_node->start_score = start_bound;
        }
    if(is_valid_end){
        bsdp_node->mask |= BSDP_Node_Mask_IS_VALID_END;
        bsdp_node->end_score = end_bound;
        }
    bsdp_node->node_data = node_data;
    bsdp_node->node_score = node_score;
    bsdp_node->stored_total = node_score;
    /* BSDP_Edge list is only created when required */
    bsdp_node->edge.list = NULL;
    bsdp_node->start_mailbox = -1;
    bsdp_node->end_mailbox = -1;
    return bsdp_node;
    }

static void BSDP_Node_destroy_pqueue_func(gpointer data,
                                          gpointer user_data){
    register BSDP *bsdp = user_data;
    register BSDP_Edge *bsdp_edge = (BSDP_Edge*)data;
    BSDP_Edge_destroy(bsdp, bsdp_edge);
    return;
    }

static void BSDP_Node_destroy(BSDP *bsdp, BSDP_Node *bsdp_node){
    register gint i;
    if(bsdp_node->mask & BSDP_Node_Mask_IS_USED){
        if(bsdp_node->edge.used) /* Can be NULL if alignment end */
            BSDP_Edge_destroy(bsdp, bsdp_node->edge.used);
    } else if(bsdp_node->mask & BSDP_Node_Mask_IS_INITIALISED){
        PQueue_destroy(bsdp_node->edge.pqueue,
                       BSDP_Node_destroy_pqueue_func, bsdp);
    } else if(bsdp_node->mask & BSDP_Node_Mask_IS_NEW){
        if(bsdp_node->edge.list){
            for(i = 0; i < bsdp_node->edge.list->len; i++)
                BSDP_Edge_destroy(bsdp, bsdp_node->edge.list->pdata[i]);
            g_ptr_array_free(bsdp_node->edge.list, TRUE);
            }
        }
    if(bsdp->destroy_node_data_func)
        bsdp->destroy_node_data_func(bsdp_node->node_data);
    g_free(bsdp_node);
    return;
    }

/**/

static gboolean BSDP_Node_pqueue_compare(gpointer low,
                                         gpointer high,
                                         gpointer user_data){
    register BSDP_Node *low_node = low, *high_node = high;
    return low_node->stored_total > high_node->stored_total;
    }

BSDP *BSDP_create(BSDP_ConfirmEdgeFunc confirm_edge_func,
                  BSDP_ConfirmTerminalFunc confirm_start_func,
                  BSDP_ConfirmTerminalFunc confirm_end_func,
                  BSDP_UpdateEdgeFunc update_edge_func,
                  BSDP_UpdateTerminalFunc update_start_func,
                  BSDP_UpdateTerminalFunc update_end_func,
                  BSDP_DestroyFunc destroy_node_data_func,
                  BSDP_DestroyFunc destroy_edge_data_func,
                  gpointer user_data){
    register BSDP *bsdp = g_new(BSDP, 1);
    g_assert(confirm_edge_func);
    g_assert(confirm_start_func);
    g_assert(confirm_end_func);
    g_assert(update_edge_func);
    g_assert(update_start_func);
    g_assert(update_end_func);
    bsdp->bas = BSDP_ArgumentSet_create(NULL);
    bsdp->confirm_edge_func = confirm_edge_func;
    bsdp->confirm_start_func = confirm_start_func;
    bsdp->confirm_end_func = confirm_end_func;
    bsdp->update_edge_func = update_edge_func;
    bsdp->update_start_func = update_start_func;
    bsdp->update_end_func = update_end_func;
    bsdp->destroy_node_data_func = destroy_node_data_func;
    bsdp->destroy_edge_data_func = destroy_edge_data_func;
    bsdp->user_data = user_data;
    bsdp->edge_recycle = RecycleBin_create("BSDP_Edge",
                                            sizeof(BSDP_Edge), 1024);
    bsdp->node_list = NULL;
    bsdp->node_pqueue = NULL;
    bsdp->path_count = 0;
    bsdp->pqueue_set = PQueueSet_create();
    bsdp->edge_filter = NULL;
    if(bsdp->bas->join_filter){
        bsdp->potential_recycle = RecycleBin_create(
            "BSDP_Potential", sizeof(BSDP_Potential), 1024);
    } else {
        bsdp->potential_recycle = NULL;
        }
    return bsdp;
    }

void BSDP_destroy(BSDP *bsdp){
    register gint i;
    if(bsdp->node_list)
        for(i = 0; i < bsdp->node_list->len; i++)
            BSDP_Node_destroy(bsdp, bsdp->node_list->pdata[i]);
    if(bsdp->node_list)
        g_ptr_array_free(bsdp->node_list, TRUE);
    RecycleBin_destroy(bsdp->edge_recycle);
    /* Don't need to free the bsdp->pqueue, as we free the set */
    PQueueSet_destroy(bsdp->pqueue_set);
    if(bsdp->edge_filter)
        g_free(bsdp->edge_filter);
    if(bsdp->potential_recycle)
        RecycleBin_destroy(bsdp->potential_recycle);
    g_free(bsdp);
    return;
    }

/**/

gint BSDP_add_node(BSDP *bsdp, gpointer node_data,
                               C4_Score node_score,
                               gboolean is_valid_start,
                               gboolean is_valid_end,
                               C4_Score start_bound,
                               C4_Score end_bound){
    register BSDP_Node *bsdp_node;
    g_assert(node_data);
    bsdp_node = BSDP_Node_create(node_data, node_score,
                                 is_valid_start, is_valid_end,
                                 start_bound, end_bound);
    if(!bsdp->node_list)
        bsdp->node_list = g_ptr_array_new();
    g_ptr_array_add(bsdp->node_list, bsdp_node);
    return bsdp->node_list->len - 1;
    }

static gboolean BSDP_Potential_pqueue_compare(gpointer low,
                                              gpointer high,
                                              gpointer user_data){
    register BSDP_Potential *low_potential= low,
                            *high_potential = high;
    return low_potential->score > high_potential->score;
    }

static BSDP_Edge_Filter *BSDP_Edge_Filter_create(BSDP *bsdp){
    register BSDP_Edge_Filter *edge_filter = g_new(BSDP_Edge_Filter, 1);
    edge_filter->src_edge_pqueue = PQueue_create(bsdp->pqueue_set,
                                BSDP_Potential_pqueue_compare, NULL);
    edge_filter->dst_edge_pqueue = PQueue_create(bsdp->pqueue_set,
                                BSDP_Potential_pqueue_compare, NULL);
    return edge_filter;
    }

static void BSDP_Edge_Filter_destroy(BSDP_Edge_Filter *edge_filter){
    PQueue_destroy(edge_filter->src_edge_pqueue, NULL, NULL);
    PQueue_destroy(edge_filter->dst_edge_pqueue, NULL, NULL);
    g_free(edge_filter);
    return;
    }

/**/

static BSDP_Potential *BSDP_Potential_create(BSDP *bsdp,
                                             BSDP_Node *bsdp_node,
                                             BSDP_Edge *bsdp_edge){
    register BSDP_Potential *bsdp_potential
        = RecycleBin_alloc(bsdp->potential_recycle);
    g_assert(bsdp);
    g_assert(bsdp_node);
    g_assert(bsdp_edge);
    bsdp_potential->ref_count = 1;
    bsdp_potential->score = bsdp_node->start_score
                          + bsdp_node->node_score
                          + bsdp_edge->join_score
                          + bsdp_edge->dst->node_score
                          + bsdp_edge->dst->end_score;
    bsdp_potential->bsdp_edge = bsdp_edge;
    bsdp_potential->src = bsdp_node;
    return bsdp_potential;
    }

static BSDP_Potential *BSDP_Potential_share(
                       BSDP_Potential *bsdp_potential){
    g_assert(bsdp_potential);
    bsdp_potential->ref_count++;
    return bsdp_potential;
    }

static void BSDP_Potential_destroy(BSDP *bsdp,
                                   BSDP_Potential *bsdp_potential){
    g_assert(bsdp_potential);
    if(--bsdp_potential->ref_count)
        return;
    BSDP_Edge_destroy(bsdp, bsdp_potential->bsdp_edge);
    RecycleBin_recycle(bsdp->potential_recycle, bsdp_potential);
    return;
    }

static BSDP_Edge *BSDP_Potential_release(BSDP *bsdp,
                                   BSDP_Potential *bsdp_potential){
    register BSDP_Edge *bsdp_edge;
    g_assert(bsdp_potential);
    bsdp_edge = bsdp_potential->bsdp_edge;
    RecycleBin_recycle(bsdp->potential_recycle, bsdp_potential);
    return bsdp_edge;
    }

static void BSDP_Potential_submit_queue(BSDP *bsdp,
            BSDP_Potential *bsdp_potential, PQueue *pqueue){
    register BSDP_Potential *top, *prev;
    /* Admit extra edge for tie-breaker removal */
    if(PQueue_total(pqueue) <= bsdp->bas->join_filter){
        PQueue_push(pqueue, bsdp_potential);
    } else {
        top = PQueue_top(pqueue);
        if(top->score < bsdp_potential->score){
            prev = PQueue_pop(pqueue);
            BSDP_Potential_destroy(bsdp, prev);
            PQueue_push(pqueue, bsdp_potential);
        } else {
            BSDP_Potential_destroy(bsdp, bsdp_potential);
            }
        }
    return;
    }

/**/

static void BSDP_Edge_submit(BSDP *bsdp,
                             BSDP_Edge *bsdp_edge,
                             gint src_node_id, gint dst_node_id){
    register BSDP_Edge_Filter *src_filter, *dst_filter;
    register BSDP_Node *src = bsdp->node_list->pdata[src_node_id];
    register BSDP_Potential *bsdp_potential
           = BSDP_Potential_create(bsdp, src, bsdp_edge);
    g_assert(bsdp);
    g_assert(bsdp_edge);
    /**/
    src_filter = bsdp->edge_filter[src_node_id];
    if(!src_filter){
        bsdp->edge_filter[src_node_id] = BSDP_Edge_Filter_create(bsdp);
        src_filter = bsdp->edge_filter[src_node_id];
        }
    g_assert(src_filter);
    /**/
    dst_filter = bsdp->edge_filter[dst_node_id];
    if(!dst_filter){
        bsdp->edge_filter[dst_node_id] = BSDP_Edge_Filter_create(bsdp);
        dst_filter = bsdp->edge_filter[dst_node_id];
        }
    g_assert(dst_filter);
    /**/
    BSDP_Potential_submit_queue(bsdp,
                           BSDP_Potential_share(bsdp_potential),
                           src_filter->src_edge_pqueue);
    BSDP_Potential_submit_queue(bsdp, bsdp_potential,
                           dst_filter->dst_edge_pqueue);
    return;
    }

void BSDP_add_edge(BSDP *bsdp, gpointer edge_data,
                   gint src_node_id, gint dst_node_id,
                   C4_Score bound_score){
    register BSDP_Edge *bsdp_edge;
    register BSDP_Node *src, *dst;
    g_assert(bsdp);
    g_assert(edge_data);
    g_assert(src_node_id >= 0);
    g_assert(src_node_id < bsdp->node_list->len);
    g_assert(dst_node_id >= 0);
    g_assert(dst_node_id < bsdp->node_list->len);
    src = bsdp->node_list->pdata[src_node_id];
    dst = bsdp->node_list->pdata[dst_node_id];
    g_assert(src);
    g_assert(dst);
    g_assert(src->mask & BSDP_Node_Mask_IS_NEW);
    g_assert(dst->mask & BSDP_Node_Mask_IS_NEW);
    bsdp_edge = BSDP_Edge_create(bsdp, edge_data, dst, bound_score);
    if(bsdp->bas->join_filter){ /* Store in PQueue for src and dst */
        if(!bsdp->edge_filter)
            bsdp->edge_filter = g_new0(BSDP_Edge_Filter*,
                                       bsdp->node_list->len);
        BSDP_Edge_submit(bsdp, bsdp_edge, src_node_id, dst_node_id);
    } else { /* Just add to edge list */
        if(!src->edge.list)
            src->edge.list = g_ptr_array_new();
        g_ptr_array_add(src->edge.list, bsdp_edge);
        }
    return;
    }

/**/

static void BSDP_update(BSDP *bsdp,
                        BSDP_Node *bsdp_node, BSDP_Edge *bsdp_edge,
                        gboolean update);

static C4_Score BSDP_Node_top_partial(BSDP *bsdp, BSDP_Node *bsdp_node,
                                      gboolean update){
    register BSDP_Edge *bsdp_edge;
    register C4_Score score = C4_IMPOSSIBLY_LOW_SCORE;
    g_assert(bsdp_node);
    g_assert(bsdp_node->mask & BSDP_Node_Mask_IS_INITIALISED);
    g_assert(!(bsdp_node->mask & BSDP_Node_Mask_IS_USED));
    bsdp_node->mask &= (~BSDP_Node_Mask_SCORED_TERMINAL);
    /* Score as end if possible */
    if(bsdp_node->mask & BSDP_Node_Mask_IS_VALID_END){
        score = bsdp_node->node_score
              + bsdp_node->end_score;
        bsdp_node->mask |= BSDP_Node_Mask_SCORED_TERMINAL;
        }
    /* Get best edge */
    do {
        bsdp_edge = PQueue_top(bsdp_node->edge.pqueue);
        if(!bsdp_edge)
            break;
        if(bsdp_edge->dst->mask & BSDP_Node_Mask_IS_USED){
            bsdp_edge = PQueue_pop(bsdp_node->edge.pqueue);
            g_assert(bsdp_edge);
            BSDP_Edge_destroy(bsdp, bsdp_edge);
        } else {
            break;
            }
    } while(TRUE);
    /* If have edge, score as edge if possible */
    if(bsdp_edge){
        g_assert(!(bsdp_edge->dst->mask & BSDP_Node_Mask_IS_USED));
        if(update){
            do {
                bsdp_edge = PQueue_pop(bsdp_node->edge.pqueue);
                if(!bsdp_edge)
                    break;
                if(bsdp_edge->dst->mask & BSDP_Node_Mask_IS_USED){
                    BSDP_Edge_destroy(bsdp, bsdp_edge);
                    continue;
                    }
                BSDP_update(bsdp, bsdp_node, bsdp_edge, TRUE);
            } while(bsdp_edge != PQueue_top(bsdp_node->edge.pqueue));
            }
        if(bsdp_edge && (score < bsdp_edge->stored_partial)){
            bsdp_node->mask &= (~BSDP_Node_Mask_SCORED_TERMINAL);
            score = bsdp_edge->stored_partial;
            }
        }
    return score;
    }

static C4_Score BSDP_Node_stored_total(BSDP *bsdp,
                                       BSDP_Node *bsdp_node,
                                       gboolean update){
    if(!(bsdp_node->mask & BSDP_Node_Mask_IS_VALID_START))
        return C4_IMPOSSIBLY_LOW_SCORE;
    return bsdp_node->start_score
         + BSDP_Node_top_partial(bsdp, bsdp_node, update);
    }

static void BSDP_update(BSDP *bsdp,
                        BSDP_Node *bsdp_node, BSDP_Edge *bsdp_edge,
                        gboolean update){
    bsdp_edge->stored_partial = bsdp_node->node_score
                              + bsdp_edge->join_score
                              + BSDP_Node_top_partial(bsdp,
                                                      bsdp_edge->dst,
                                                      update);
    g_assert(bsdp_node->mask & BSDP_Node_Mask_IS_INITIALISED);
    g_assert(!(bsdp_node->mask & BSDP_Node_Mask_IS_USED));
    PQueue_push(bsdp_node->edge.pqueue, bsdp_edge);
    return;
    }

static gboolean BSDP_Edge_pqueue_compare(gpointer low,
                                         gpointer high,
                                         gpointer user_data){
    register BSDP_Edge *low_edge = low, *high_edge = high;
    return low_edge->stored_partial > high_edge->stored_partial;
    }

static void BSDP_initialise_recur(BSDP *bsdp, BSDP_Node *bsdp_node){
    register gint i;
    register BSDP_Edge *bsdp_edge;
    register GPtrArray *edge_list;
    if(bsdp_node->mask & BSDP_Node_Mask_IS_INITIALISED)
        return;
    edge_list = bsdp_node->edge.list;
    bsdp_node->edge.pqueue = PQueue_create(bsdp->pqueue_set,
                                       BSDP_Edge_pqueue_compare, NULL);
    bsdp_node->mask &= (~BSDP_Node_Mask_IS_NEW);
    bsdp_node->mask |= BSDP_Node_Mask_IS_INITIALISED;
    if(edge_list){
        for(i = 0; i < edge_list->len; i++){
            bsdp_edge = edge_list->pdata[i];
            BSDP_initialise_recur(bsdp, bsdp_edge->dst);
            BSDP_update(bsdp, bsdp_node, bsdp_edge, FALSE);
            }
        g_ptr_array_free(edge_list, TRUE);
        }
    return;
    }

static void BSDP_initialise_remove_tiebreakers(BSDP *bsdp,
                                               PQueue *pqueue){
    register BSDP_Potential *bsdp_potential;
    register C4_Score score;
    g_assert(pqueue);
    if(PQueue_total(pqueue) > bsdp->bas->join_filter){
        bsdp_potential = PQueue_pop(pqueue);
        score = bsdp_potential->score;
        BSDP_Potential_destroy(bsdp, bsdp_potential);
        while(PQueue_total(pqueue)){
            bsdp_potential = PQueue_top(pqueue);
            g_assert(bsdp_potential);
            if(score != bsdp_potential->score)
                break;
            bsdp_potential = PQueue_pop(pqueue);
            BSDP_Potential_destroy(bsdp, bsdp_potential);
            }
        }
    return;
    }

static void BSDP_initialise_filter(BSDP *bsdp, PQueue *pqueue){
    register BSDP_Potential *bsdp_potential;
    register BSDP_Node *bsdp_node;
    g_assert(pqueue);
    while((bsdp_potential = PQueue_pop(pqueue))){
        if(bsdp_potential->ref_count == 2){ /* In src + dst */
            bsdp_node = bsdp_potential->src;
            if(!bsdp_node->edge.list)
                bsdp_node->edge.list = g_ptr_array_new();
            g_ptr_array_add(bsdp_node->edge.list,
                            bsdp_potential->bsdp_edge);
            bsdp_potential->ref_count = 0;
        } else {
            if(bsdp_potential->ref_count)
                BSDP_Potential_destroy(bsdp, bsdp_potential);
            else
                BSDP_Potential_release(bsdp, bsdp_potential);
            }
        }
    return;
    }

void BSDP_initialise(BSDP *bsdp, C4_Score threshold){
    register gint i;
    register BSDP_Node *bsdp_node;
    register BSDP_Edge_Filter *edge_filter;
    if(!bsdp->node_list)
        return;
    if(bsdp->edge_filter){
        for(i = 0; i < bsdp->node_list->len; i++){
            edge_filter = bsdp->edge_filter[i];
            if(edge_filter)
                BSDP_initialise_remove_tiebreakers(bsdp,
                                   edge_filter->src_edge_pqueue);
            }
        for(i = 0; i < bsdp->node_list->len; i++){
            edge_filter = bsdp->edge_filter[i];
            if(edge_filter){
                BSDP_initialise_filter(bsdp,
                                   edge_filter->src_edge_pqueue);
                BSDP_initialise_filter(bsdp,
                                   edge_filter->dst_edge_pqueue);
                BSDP_Edge_Filter_destroy(edge_filter);
                }
            }
        g_free(bsdp->edge_filter);
        bsdp->edge_filter = NULL;
        RecycleBin_destroy(bsdp->potential_recycle);
        bsdp->potential_recycle = NULL;
        }
    for(i = 0; i < bsdp->node_list->len; i++){ /* Initialise scores */
        bsdp_node = bsdp->node_list->pdata[i];
        BSDP_initialise_recur(bsdp, bsdp_node);
        bsdp_node->stored_total = BSDP_Node_stored_total(bsdp,
                                                     bsdp_node, FALSE);
        if(bsdp_node->stored_total >= threshold){
            if(!bsdp->node_pqueue)
                bsdp->node_pqueue = PQueue_create(bsdp->pqueue_set,
                                    BSDP_Node_pqueue_compare, NULL);
            PQueue_push(bsdp->node_pqueue, bsdp_node);
            }
        }
    return;
    }

/* Score calculations:
 * ------------------
 *
 * Top_partial(node) = Max(node_score + node_end_score,
 *                         Top(node->edge.pqueue) );
 *
 * Edge: stored_partial = src_node_score
 *                      + edge_score
 *                      + Top_Partial( dst_node )
 *
 * Node: stored_total   = start_score
 *                      + Top_Partial( node )
 */

/**/

static void BSDP_path_validate_recur(BSDP *bsdp,
                                     BSDP_Node *bsdp_node){
    register BSDP_Edge *bsdp_edge;
    g_assert(bsdp_node->mask & BSDP_Node_Mask_IS_INITIALISED);
    g_assert(!(bsdp_node->mask & BSDP_Node_Mask_IS_USED));
    if(bsdp_node->mask & BSDP_Node_Mask_SCORED_TERMINAL)
        return;
    do {
        bsdp_edge = PQueue_pop(bsdp_node->edge.pqueue);
        if(bsdp_edge){
            if(bsdp_edge->dst->mask & BSDP_Node_Mask_IS_USED){
                BSDP_Edge_destroy(bsdp, bsdp_edge);
                continue;
                }
            BSDP_path_validate_recur(bsdp, bsdp_edge->dst);
            BSDP_update(bsdp, bsdp_node, bsdp_edge, FALSE);
            }
    } while(PQueue_top(bsdp_node->edge.pqueue) != bsdp_edge);
    return;
    }

static gboolean BSDP_path_is_validated(BSDP *bsdp,
                                       C4_Score threshold){
    register BSDP_Node *bsdp_node = PQueue_top(bsdp->node_pqueue);
    register BSDP_Edge *bsdp_edge;
    register C4_Score score, check_score;
    g_assert(bsdp_node);
    g_assert(bsdp_node->mask & BSDP_Node_Mask_IS_VALID_START);
    score = bsdp_node->stored_total;
    check_score = bsdp_node->start_score;
    do {
        g_assert(bsdp_node);
        g_assert(bsdp_node->mask & BSDP_Node_Mask_IS_INITIALISED);
        g_assert(!(bsdp_node->mask & BSDP_Node_Mask_IS_USED));
        check_score += bsdp_node->node_score;
        if((bsdp_node->mask & BSDP_Node_Mask_SCORED_TERMINAL)){
            check_score += bsdp_node->end_score;
            break;
            }
        bsdp_edge = PQueue_top(bsdp_node->edge.pqueue);
        g_assert(bsdp_edge);
        check_score += bsdp_edge->join_score;
        bsdp_node = bsdp_edge->dst;
    } while(bsdp_node);
    g_assert(check_score == score);
    g_assert(score >= threshold);
    return TRUE;
    }

static gboolean BSDP_path_validate(BSDP *bsdp, C4_Score threshold){
    register BSDP_Node *bsdp_node;
    if(bsdp->node_pqueue){
        do {
            bsdp_node = PQueue_pop(bsdp->node_pqueue);
            if(!bsdp_node)
                return FALSE;
            if(bsdp_node->mask & BSDP_Node_Mask_IS_USED)
                continue;
            BSDP_path_validate_recur(bsdp, bsdp_node);
            bsdp_node->stored_total = BSDP_Node_stored_total(bsdp,
                                               bsdp_node, TRUE);
            if(bsdp_node->stored_total >= threshold){
                PQueue_push(bsdp->node_pqueue, bsdp_node);
            } else {
                continue;
                }
        } while(PQueue_top(bsdp->node_pqueue) != bsdp_node);
        g_assert(BSDP_path_is_validated(bsdp, threshold));
        return TRUE;
        }
    return FALSE;
    }

static gint BSDP_path_confirm(BSDP *bsdp){
    register BSDP_Node *first_node = PQueue_top(bsdp->node_pqueue),
                       *bsdp_node = first_node;
    register BSDP_Edge *bsdp_edge;
    register C4_Score prev_score, confirmed_score;
    register gint confirm_count = 0;
    do {
        g_assert(bsdp_node->node_data);
        g_assert(bsdp_node->mask & BSDP_Node_Mask_IS_INITIALISED);
        g_assert(!(bsdp_node->mask & BSDP_Node_Mask_IS_NEW));
        if(bsdp_node->mask & BSDP_Node_Mask_SCORED_TERMINAL)
            break;
        bsdp_edge = PQueue_top(bsdp_node->edge.pqueue);
        if(!bsdp_edge)
            break;
        if(bsdp_edge->mailbox == -1){ /* Not confirmed */
            bsdp_edge->mailbox = bsdp->path_count;
            confirmed_score = bsdp->confirm_edge_func(
                                  bsdp_node->node_data,
                                  bsdp_edge->edge_data,
                                  bsdp_edge->dst->node_data,
                                  bsdp->user_data);
            g_assert(bsdp_edge->join_score >= confirmed_score);
            if(bsdp_edge->join_score != confirmed_score){
                bsdp_edge->join_score = confirmed_score;
                confirm_count++;
                }
        } else { /* Already confirmed */
            /* Check for subopt edge clashes if out of date */
            if(bsdp_edge->mailbox != bsdp->path_count){
                prev_score = bsdp_edge->join_score;
                bsdp_edge->join_score = bsdp->update_edge_func(
                                              bsdp_node->node_data,
                                              bsdp_edge->edge_data,
                                              bsdp_edge->dst->node_data,
                                              bsdp->user_data,
                                              prev_score,
                                              bsdp_edge->mailbox);
                g_assert(bsdp_edge->join_score <= prev_score);
                bsdp_edge->mailbox = bsdp->path_count;
                if(bsdp_edge->join_score != prev_score)
                    confirm_count++;
                }
            }
        bsdp_node = bsdp_edge->dst;
    } while(bsdp_node);
    /* Confirm the start */
    g_assert(bsdp_node);
    if(first_node->mask & BSDP_Node_Mask_CONFIRMED_START){
        /* Check for subopt start clashes */
        if(first_node->start_mailbox != bsdp->path_count){
            prev_score = first_node->start_score;
            first_node->start_score = bsdp->update_start_func(
                    first_node->node_data, bsdp->user_data,
                    prev_score, first_node->start_mailbox);
            g_assert(first_node->start_score <= prev_score);
            first_node->start_mailbox = bsdp->path_count;
            if(first_node->start_score != prev_score)
                confirm_count++;
            }
    } else {
        first_node->start_mailbox = bsdp->path_count;
        confirmed_score = bsdp->confirm_start_func(
                              first_node->node_data, bsdp->user_data);
        first_node->mask |= BSDP_Node_Mask_CONFIRMED_START;
        g_assert(first_node->start_score >= confirmed_score);
        if(first_node->start_score != confirmed_score){
            first_node->start_score = confirmed_score;
            confirm_count++;
            }
        }
    /* Chose the end */
    if(bsdp_node->mask & BSDP_Node_Mask_CONFIRMED_END){
        /* Check for subopt end clashes  */
        if(bsdp_node->end_mailbox != bsdp->path_count){
            prev_score = bsdp_node->end_score;
            bsdp_node->end_score = bsdp->update_end_func(
                    bsdp_node->node_data, bsdp->user_data,
                    prev_score, bsdp_node->end_mailbox);
            g_assert(bsdp_node->end_score <= prev_score);
            bsdp_node->end_mailbox = bsdp->path_count;
            if(bsdp_node->end_score != prev_score)
                confirm_count++;
            }
    } else {
        bsdp_node->end_mailbox = bsdp->path_count;
        confirmed_score = bsdp->confirm_end_func(
                              bsdp_node->node_data, bsdp->user_data);
        bsdp_node->mask |= BSDP_Node_Mask_CONFIRMED_END;
        g_assert(bsdp_node->end_score >= confirmed_score);
        if(bsdp_node->end_score != confirmed_score){
            bsdp_node->end_score = confirmed_score;
            confirm_count++;
            }
        }
    return confirm_count;
    }

static BSDP_Path *BSDP_Path_create(void){
    register BSDP_Path *bsdp_path = g_new(BSDP_Path, 1);
    bsdp_path->node_list = g_ptr_array_new();
    return bsdp_path;
    }

void BSDP_Path_destroy(BSDP_Path *bsdp_path){
    g_ptr_array_free(bsdp_path->node_list, TRUE);
    g_free(bsdp_path);
    return;
    }

static gboolean BSDP_Path_check(BSDP_Path *bsdp_path){
    register C4_Score check_score = 0;
    register BSDP_Node *start_node, *end_node, *bsdp_node;
    register gint i;
    g_assert(bsdp_path->node_list->len);
    /* Check start node */
    start_node = bsdp_path->node_list->pdata[0];
    g_assert(start_node->mask & BSDP_Node_Mask_IS_VALID_START);
    g_assert(start_node->mask & BSDP_Node_Mask_CONFIRMED_START);
    g_assert(start_node->mask & BSDP_Node_Mask_USED_AS_START);
    check_score = start_node->start_score;
    /* Check end node */
    end_node = bsdp_path->node_list->pdata[bsdp_path->node_list->len-1];
    g_assert(end_node->mask & BSDP_Node_Mask_IS_VALID_END);
    g_assert(end_node->mask & BSDP_Node_Mask_CONFIRMED_END);
    g_assert(end_node->mask & BSDP_Node_Mask_USED_AS_END);
    check_score += end_node->end_score;
    /* Check all other nodes */
    for(i = 0; i < bsdp_path->node_list->len; i++){
        bsdp_node = bsdp_path->node_list->pdata[i];
        g_assert(!(bsdp_node->mask & BSDP_Node_Mask_IS_NEW));
        g_assert(bsdp_node->mask & BSDP_Node_Mask_IS_INITIALISED);
        g_assert(bsdp_node->mask & BSDP_Node_Mask_IS_USED);
        check_score += bsdp_node->node_score;
        if(bsdp_node != end_node){
            check_score += bsdp_node->edge.used->join_score;
            }
        }
    g_assert(check_score == bsdp_path->score);
    return TRUE;
    }

static BSDP_Path *BSDP_path_extract(BSDP *bsdp){
    register BSDP_Path *bsdp_path = BSDP_Path_create();
    register BSDP_Node *bsdp_node = PQueue_top(bsdp->node_pqueue);
    register BSDP_Edge *bsdp_edge;
    bsdp_path->score = bsdp_node->stored_total;
    g_assert(bsdp_node->mask & BSDP_Node_Mask_CONFIRMED_START);
    bsdp_node->mask |= BSDP_Node_Mask_USED_AS_START;
    do {
        g_assert(bsdp_node);
        g_assert(bsdp_node->mask & BSDP_Node_Mask_IS_INITIALISED);
        g_assert(!(bsdp_node->mask & BSDP_Node_Mask_IS_USED));
        g_ptr_array_add(bsdp_path->node_list, bsdp_node);
        bsdp_node->mask |= BSDP_Node_Mask_IS_USED;
        bsdp_edge = PQueue_pop(bsdp_node->edge.pqueue);
        PQueue_destroy(bsdp_node->edge.pqueue,
                       BSDP_Node_destroy_pqueue_func, bsdp);
        bsdp_node->edge.used = bsdp_edge;
        if(bsdp_node->mask & BSDP_Node_Mask_SCORED_TERMINAL){
            g_assert(bsdp_node->mask & BSDP_Node_Mask_CONFIRMED_END);
            bsdp_node->mask |= BSDP_Node_Mask_USED_AS_END;
            break;
        } else {
            g_assert(bsdp_edge);
            }
        bsdp_node = bsdp_edge->dst;
    } while(bsdp_node);
    g_assert(BSDP_Path_check(bsdp_path));
    return bsdp_path;
    }

BSDP_Path *BSDP_next_path(BSDP *bsdp, C4_Score threshold){
    register BSDP_Path *bsdp_path;
    do {
        if(!BSDP_path_validate(bsdp, threshold))
            return NULL;
    } while(BSDP_path_confirm(bsdp));
    bsdp_path = BSDP_path_extract(bsdp);
    bsdp->path_count++;
    return bsdp_path;
    }

/**/

