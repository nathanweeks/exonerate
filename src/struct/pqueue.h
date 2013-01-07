/****************************************************************\
*                                                                *
*  Priority queue library using pairing heaps.                   *
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

#ifndef INCLUDED_PQUEUE_H
#define INCLUDED_PQUEUE_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <glib.h>
#include "recyclebin.h"

typedef struct PQueue_Node {
              gpointer  data;
    struct PQueue_Node *left;
    struct PQueue_Node *next;
    struct PQueue_Node *prev;
} PQueue_Node;

typedef void (*PQueue_Free_Func)(gpointer data, gpointer user_data);
typedef gboolean (*PQueue_Compare_Func)(gpointer low, gpointer high,
                                        gpointer user_data);

typedef struct {
    RecycleBin *pq_recycle;
    RecycleBin *node_recycle;
} PQueueSet;

typedef struct {
            PQueue_Node  *root;      /* The root node       */
                   gint   total;     /* Number of members   */
    PQueue_Compare_Func   comp_func; /* Comparison function */
              PQueueSet  *set;
               gpointer   user_data;
} PQueue;

PQueueSet *PQueueSet_create(void);
     void  PQueueSet_destroy(PQueueSet *pq_set);

     PQueue *PQueue_create(PQueueSet *pq_set,
                           PQueue_Compare_Func comp_func,
                           gpointer user_data);
       void  PQueue_destroy(PQueue *pq, PQueue_Free_Func free_func,
                            gpointer user_data);
PQueue_Node *PQueue_push(PQueue *pq, gpointer data);
       void  PQueue_raise(PQueue *pq, PQueue_Node *pqn);
       void  PQueue_change(PQueue *pq, PQueue_Node *pqn);
   gpointer  PQueue_pop(PQueue *pq);

 PQueue *PQueue_join(PQueue *a, PQueue *b);
gpointer PQueue_remove(PQueue *pq, PQueue_Node *pqn);
typedef gboolean (*PQueue_Traverse_Func)(gpointer data,
                                         gpointer user_data);
/* Returns TRUE to stop the traversal */

void PQueue_traverse(PQueue *pq, PQueue_Traverse_Func tf,
                     gpointer user_data);

#define PQueue_total(pq) ((pq)->total)

#define PQueue_top(pq) \
        ((pq)?(((pq)->root)?((pq)->root->data):(NULL)):(NULL))

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* INCLUDED_PQUEUE_H */

