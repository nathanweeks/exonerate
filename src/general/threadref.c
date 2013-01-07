/****************************************************************\
*                                                                *
*  Basic library for thread-safe reference counting.             *
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

#include "threadref.h"

ThreadRef *ThreadRef_create(void){
    register ThreadRef *threadref = g_new(ThreadRef, 1);
    threadref->ref_count = 1;
#ifdef USE_PTHREADS
    pthread_mutex_init(&threadref->ref_lock, NULL);
#endif /* USE_PTHREADS */
    return threadref;
    }

ThreadRef *ThreadRef_destroy(ThreadRef *threadref){
    ThreadRef_lock(threadref);
    if(--threadref->ref_count){
        ThreadRef_unlock(threadref);
        return threadref;
        }
    ThreadRef_unlock(threadref);
#ifdef USE_PTHREADS
    pthread_mutex_destroy(&threadref->ref_lock);
#endif /* USE_PTHREADS */
    g_free(threadref);
    return NULL;
    }

ThreadRef *ThreadRef_share(ThreadRef *threadref){
    ThreadRef_lock(threadref);
    threadref->ref_count++;
    ThreadRef_unlock(threadref);
    return threadref;
    }

gint ThreadRef_get_count(ThreadRef *threadref){
    register gint ref_count;
    ThreadRef_lock(threadref);
    ref_count = threadref->ref_count;
    ThreadRef_unlock(threadref);
    return ref_count;
    }

void ThreadRef_lock(ThreadRef *threadref){
#ifdef USE_PTHREADS
    pthread_mutex_lock(&threadref->ref_lock);
#endif /* USE_PTHREADS */
    return;
    }

void ThreadRef_unlock(ThreadRef *threadref){
#ifdef USE_PTHREADS
    pthread_mutex_unlock(&threadref->ref_lock);
#endif /* USE_PTHREADS */
    return;
    }

