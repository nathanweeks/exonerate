/****************************************************************\
*                                                                *
*  BSDP: Bounded Sparse Dynamic Programming                      *
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

#ifndef INCLUDED_BSDP_H
#define INCLUDED_BSDP_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <glib.h>

#include "c4.h"
#include "pqueue.h"
#include "recyclebin.h"
#include "argument.h"

typedef struct {
    gint join_filter;
} BSDP_ArgumentSet;

BSDP_ArgumentSet *BSDP_ArgumentSet_create(Argument *arg);

typedef struct {
            gpointer  edge_data;
    struct BSDP_Node *dst;
            C4_Score  join_score;     /* Current edge cost       */
            C4_Score  stored_partial; /* Score for edge.pqueue   */
                gint  mailbox;        /* 0 indicates confirmed   */
} BSDP_Edge;
/* The confirm cost currently used is (DP_cells * states)
 */

typedef enum {
    BSDP_Node_Mask_IS_NEW          = (1<<0),
    BSDP_Node_Mask_IS_INITIALISED  = (1<<1),
    BSDP_Node_Mask_IS_USED         = (1<<2),
    BSDP_Node_Mask_IS_VALID_START  = (1<<3),
    BSDP_Node_Mask_IS_VALID_END    = (1<<4),
    BSDP_Node_Mask_CONFIRMED_START = (1<<5),
    BSDP_Node_Mask_CONFIRMED_END   = (1<<6),
    BSDP_Node_Mask_USED_AS_START   = (1<<7),
    BSDP_Node_Mask_USED_AS_END     = (1<<8),
    BSDP_Node_Mask_SCORED_TERMINAL = (1<<9)
} BSDP_Node_Mask;

typedef struct BSDP_Node {
    BSDP_Node_Mask  mask;
          gpointer  node_data;
          C4_Score  node_score;     /* Score for HSP only       */
          C4_Score  stored_total;   /* Score for node_pqueue    */
             union {
                 GPtrArray *list;   /* Used once IS_NEW         */
                   PQueue  *pqueue; /* Used once IS_INITIALISED */
                BSDP_Edge  *used;   /* Used once IS_USED        */
          } edge;
          C4_Score  start_score;
          C4_Score  end_score;
              gint  start_mailbox;
              gint  end_mailbox;
} BSDP_Node;

typedef C4_Score (*BSDP_ConfirmEdgeFunc)(gpointer src_data,
                                         gpointer edge_data,
                                         gpointer dst_data,
                                         gpointer user_data);

typedef C4_Score (*BSDP_UpdateEdgeFunc)(gpointer src_data,
                                        gpointer edge_data,
                                        gpointer dst_data,
                                        gpointer user_data,
                                        C4_Score prev_score,
                                        gint last_updated);

typedef C4_Score (*BSDP_ConfirmTerminalFunc)(gpointer node_data,
                                             gpointer user_data);

typedef C4_Score (*BSDP_UpdateTerminalFunc)(gpointer node_data,
                                            gpointer user_data,
                                            C4_Score prev_score,
                                            gint last_updated);

typedef gboolean (*BSDP_NodeIsUsedFunc)(gpointer node_data,
                                        gpointer user_data);
typedef void (*BSDP_UseNodeFunc)(gpointer node_data,
                                 gpointer user_data);

typedef void (*BSDP_DestroyFunc)(gpointer data);

typedef struct {
         gint  ref_count;
     C4_Score  score;
    BSDP_Edge *bsdp_edge;
    BSDP_Node *src;
} BSDP_Potential;

typedef struct {
    PQueue *src_edge_pqueue;
    PQueue *dst_edge_pqueue;
} BSDP_Edge_Filter;

typedef struct {
            BSDP_ArgumentSet  *bas;
        BSDP_ConfirmEdgeFunc   confirm_edge_func;
    BSDP_ConfirmTerminalFunc   confirm_start_func;
    BSDP_ConfirmTerminalFunc   confirm_end_func;
         BSDP_UpdateEdgeFunc   update_edge_func;
     BSDP_UpdateTerminalFunc   update_start_func;
     BSDP_UpdateTerminalFunc   update_end_func;
    /**/
            BSDP_DestroyFunc   destroy_node_data_func;
            BSDP_DestroyFunc   destroy_edge_data_func;
                    gpointer   user_data;
                   GPtrArray  *node_list;
                      PQueue  *node_pqueue;
                  RecycleBin  *edge_recycle;
                       guint   path_count;
                   PQueueSet  *pqueue_set;
            BSDP_Edge_Filter **edge_filter;
                  RecycleBin  *potential_recycle;
} BSDP;

/**/

BSDP *BSDP_create(BSDP_ConfirmEdgeFunc confirm_edge_func,
                  BSDP_ConfirmTerminalFunc confirm_start_func,
                  BSDP_ConfirmTerminalFunc confirm_end_func,
                  BSDP_UpdateEdgeFunc update_edge_func,
                  BSDP_UpdateTerminalFunc update_start_func,
                  BSDP_UpdateTerminalFunc update_end_func,
                  BSDP_DestroyFunc destroy_node_data_func,
                  BSDP_DestroyFunc destroy_edge_data_func,
                  gpointer user_data);
void BSDP_destroy(BSDP *bsdp);

gint BSDP_add_node(BSDP *bsdp,
                   gpointer node_data,
                   C4_Score node_score,
                   gboolean is_valid_start,
                   gboolean is_valid_end,
                   C4_Score start_bound,
                   C4_Score end_bound);
void BSDP_add_edge(BSDP *bsdp, gpointer edge_data,
                   gint src_node_id, gint dst_node_id,
                   C4_Score bound_score);

void  BSDP_initialise(BSDP *bsdp, C4_Score threshold);
/* This should be called after all the nodes and edges have been added
 */

typedef struct {
    C4_Score  score;
   GPtrArray *node_list;
} BSDP_Path;

BSDP_Path *BSDP_next_path(BSDP *bsdp, C4_Score threshold);
     void  BSDP_Path_destroy(BSDP_Path *bsdp_path);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* INCLUDED_BSDP_H */

