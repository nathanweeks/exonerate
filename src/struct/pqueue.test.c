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

#include "pqueue.h"
#include <string.h> /* For strcmp() */


typedef struct {
   gchar *word;
    gint  real_priority;
    gint  priority1;
    gint  priority2;
    gint  priority3;
} Data;

#define DATA_TOTAL 12

static gboolean test_comp_low_int(gpointer low, gpointer high,
                                  gpointer user_data){
    register Data *low_d = (Data*)low,
                  *high_d = (Data*)high;
    return low_d->real_priority < high_d->real_priority;
    }

static gboolean test_traverse_func(gpointer data, gpointer user_data){
    register Data *d = (Data*)data;
    g_print("Traverse [%s]\n", d->word);
    return FALSE;
    }

int main(void){
    register guint i;
    register PQueueSet *pq_set = PQueueSet_create();
    register PQueue *pq = PQueue_create(pq_set, test_comp_low_int,
                                        NULL);
    register Data *dn;
    PQueue_Node *acc[DATA_TOTAL];
    Data data[DATA_TOTAL] = {
        { "REMOVE-1", -1, 300, 200,  50},
        { "jumps", -1, 39, 38, 32}, { "dog",   -1, 72, 61, 60},
        { "brown", -1, 99, 40,  8}, { "over",  -1, 68, 35, 33},
        { "REMOVE-2", -1,  25,  20,  15},
        { "fox",   -1, 20, 19, 17}, { "quick", -1, 99, 19,  3},
        { "lazy",  -1, 99, 67, 51}, { "the",   -1, 3,   2,  1},
        { "a",     -1, 82, 41, 34},
        { "REMOVE-3", -1,  40,  30,  20}
    };
    gint r[3] = {0, 5, 11};
    gchar *answer[DATA_TOTAL-3] = {
        "the", "quick", "brown", "fox", "jumps",
        "over", "a", "lazy", "dog"};
    g_print("Adding words with initial priorities\n");
    for(i = 0; i < DATA_TOTAL; i++){
        g_print("Adding [%s]\n", data[i].word);
        data[i].real_priority = data[i].priority1;
        acc[i] = PQueue_push(pq, (gpointer)&data[i]);
        }
    g_print("Set intermediate priorities\n");
    for(i = 0; i < DATA_TOTAL; i++){
        g_print("Raising [%s] from %d to %d\n",
              data[i].word, data[i].priority1, data[i].priority2);
        data[i].real_priority = data[i].priority2;
        PQueue_raise(pq, acc[i]);
        }
    g_print("Set final priorities\n");
    for(i = 0; i < DATA_TOTAL; i++){
        g_print("Raising [%s] from %d to %d\n",
              data[i].word, data[i].priority2, data[i].priority3);
        data[i].real_priority = data[i].priority3;
        PQueue_raise(pq, acc[i]);
        }
    g_print("Testing traverse\n");
    PQueue_traverse(pq, test_traverse_func, NULL);
    for(i = 0; i < 3; i++){
        g_print("Testing remove of \"%s\"\n", data[r[i]].word);
        PQueue_remove(pq, acc[r[i]]);
        g_print("Testing traverse\n");
        PQueue_traverse(pq, test_traverse_func, NULL);
        }
    g_print("Popping data\n");
    for(i = 0; i < (DATA_TOTAL-3); i++){
        dn = PQueue_pop(pq);
        g_print("Word: [%s] Answer: [%s]\n", dn->word, answer[i]);
        g_assert(!strcmp(dn->word, answer[i]));
        }
    PQueue_destroy(pq, NULL, NULL);
    PQueueSet_destroy(pq_set);
    return 0;
    }


