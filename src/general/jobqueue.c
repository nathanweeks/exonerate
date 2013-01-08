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

#include "jobqueue.h"
#include <unistd.h> /* For usleep() */

#ifdef USE_PTHREADS
static JobQueue_Task *JobQueue_Task_create(JobQueue_Func job_func,
                                           gpointer job_data, gint priority){
    register JobQueue_Task *task = g_new(JobQueue_Task, 1);
    task->job_func = job_func;
    task->job_data = job_data;
    task->priority = priority;
    return task;
    }

static void JobQueue_Task_destroy(JobQueue_Task *task){
    g_free(task);
    return;
    }

static void *JobQueue_thread_func(void *data){
    register JobQueue *jq = data;
    register JobQueue_Task *task;
    do {
        pthread_mutex_lock(&jq->queue_lock);
        task = PQueue_pop(jq->pq);
        if(task){
            jq->running_count++;
            pthread_mutex_unlock(&jq->queue_lock);
            task->job_func(task->job_data);
            JobQueue_Task_destroy(task);
            pthread_mutex_lock(&jq->queue_lock);
            jq->running_count--;
            pthread_mutex_unlock(&jq->queue_lock);
       } else {
            pthread_mutex_unlock(&jq->queue_lock);
            usleep(1000); /* wait before checking job queue again */
            }
    } while(!jq->is_complete);
    return NULL;
    }

static gboolean JobQueue_Task_compare(gpointer low, gpointer high,
                                      gpointer user_data){
    register JobQueue_Task *low_task = (JobQueue_Task*)low,
                           *high_task = (JobQueue_Task*)high;
    return low_task->priority < high_task->priority;
    }
#endif /* USE_PTHREADS */

JobQueue *JobQueue_create(gint thread_total){
    register JobQueue *jq = g_new(JobQueue, 1);
#ifdef USE_PTHREADS
    register gint i;
    jq->is_complete = FALSE;
    jq->thread_total = thread_total;
    jq->running_count = 0;
    pthread_mutex_init(&jq->queue_lock, NULL);
    jq->thread = g_new0(pthread_t, thread_total);
    jq->pq_set = PQueueSet_create();
    jq->pq = PQueue_create(jq->pq_set, JobQueue_Task_compare, NULL);
    for(i = 0; i < thread_total; i++)
        pthread_create(&jq->thread[i], NULL, JobQueue_thread_func, jq);
#endif /* USE_PTHREADS */
    return jq;
    }

#ifdef USE_PTHREADS
static void JobQueue_Task_free_func(gpointer data, gpointer user_data){
    register JobQueue_Task *task = data;
    JobQueue_Task_destroy(task);
    return;
    }
#endif

void JobQueue_destroy(JobQueue *jq){
#ifdef USE_PTHREADS
    g_assert(!PQueue_total(jq->pq));
    PQueue_destroy(jq->pq, JobQueue_Task_free_func, NULL);
    PQueueSet_destroy(jq->pq_set);
    pthread_mutex_destroy(&jq->queue_lock);
    g_free(jq->thread);
#endif /* USE_PTHREADS */
    g_free(jq);
    return;
    }

void JobQueue_submit(JobQueue *jq, JobQueue_Func job_func,
                     gpointer job_data, gint priority){
#ifdef USE_PTHREADS
    register JobQueue_Task *task;
    g_assert(jq);
    pthread_mutex_lock(&jq->queue_lock);
    task = JobQueue_Task_create(job_func, job_data, priority);
    PQueue_push(jq->pq, task);
    pthread_mutex_unlock(&jq->queue_lock);
#else /* USE_PTHREADS */
    job_func(job_data); /* when no threads available, just run the job */
#endif /* USE_PTHREADS */
    return;
    }
/* waits until less than (thread_total * 2) jobs are queued
 */

void JobQueue_complete(JobQueue *jq){
#ifdef USE_PTHREADS
    register gint i, job_total;
    do { /* wait for job queue to empty */
        pthread_mutex_lock(&jq->queue_lock);
        job_total = PQueue_total(jq->pq) + jq->running_count;
        pthread_mutex_unlock(&jq->queue_lock);
        usleep(1000); /* wait before checking if queue is empty again */
    } while(job_total);
    jq->is_complete = TRUE;
    /* wait for threads to finish */
    for(i = 0; i < jq->thread_total; i++)
        pthread_join(jq->thread[i], NULL);
#endif /* USE_PTHREADS */
    return;
    }
/* waits until the queue is empty and no jobs are running,
 * to catch any jobs which are added to the queue by runnning jobs
 */

/**/

