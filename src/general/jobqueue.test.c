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

#include <unistd.h> /* For usleep() */
#include "jobqueue.h"

static void test_run_func(gpointer job_data){
    register gchar *name = job_data;
    g_message("running job [%s]", name);
    usleep(1000000); /* test job length */
    return;
    }

int main(void){
    register gint i;
    register JobQueue *jq = JobQueue_create(2);
    gchar *test_str[10] = {"one", "two", "three", "four", "five",
                           "six", "seven", "eight", "nine", "ten"};

    /* Initialise threads as not using Argument library */
#ifdef USE_PTHREADS
    if(!g_thread_supported())
        g_thread_init(NULL);
#endif /* USE_PTHREADS */

    for(i = 0; i < 10; i++)
        JobQueue_submit(jq, test_run_func, test_str[i], 0);
    JobQueue_complete(jq);
    JobQueue_destroy(jq);
    return 0;
    }

