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

#ifndef INCLUDED_THREADREF_H
#define INCLUDED_THREADREF_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <glib.h>

#ifdef USE_PTHREADS
#include <pthread.h>
#endif /* USE_PTHREADS */


typedef struct {
               gint ref_count;
#ifdef USE_PTHREADS
    pthread_mutex_t ref_lock;
#endif /* USE_PTHREADS */
} ThreadRef;

ThreadRef  *ThreadRef_create(void);
ThreadRef  *ThreadRef_destroy(ThreadRef *threadref);
ThreadRef  *ThreadRef_share(ThreadRef *threadref);
     gint   ThreadRef_get_count(ThreadRef *threadref);
     void   ThreadRef_lock(ThreadRef *threadref);
     void   ThreadRef_unlock(ThreadRef *threadref);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* INCLUDED_THREADREF_H */

