/****************************************************************\
*                                                                *
*  Library for Queueing Multi-Threaded Jobs                      *
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

#ifndef INCLUDED_JOBQUEUE_H
#define INCLUDED_JOBQUEUE_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <glib.h>

#ifdef USE_PTHREADS
#include <pthread.h>
#endif /* USE_PTHREADS */

#include "pqueue.h"

typedef void (*JobQueue_Func)(gpointer job_data);

typedef struct JobQueue_Task {
    JobQueue_Func  job_func;
         gpointer  job_data;
             gint  priority;
} JobQueue_Task;

typedef struct {
#ifdef USE_PTHREADS
           gboolean  is_complete;
               gint  thread_total;
               gint  running_count;
    pthread_mutex_t  queue_lock;
          pthread_t *thread;
             PQueue *pq;
          PQueueSet *pq_set;
#endif /* USE_PTHREADS */
} JobQueue;

JobQueue *JobQueue_create(gint thread_total);
    void  JobQueue_destroy(JobQueue *jq);
    void  JobQueue_submit(JobQueue *jq, JobQueue_Func job_func,
                          gpointer job_data, gint priority);
    void  JobQueue_complete(JobQueue *jq);

/* Lower priority jobs are run first */

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* INCLUDED_JOBQUEUE_H */


/**/

