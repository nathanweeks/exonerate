/****************************************************************\
*                                                                *
*  Library for PCR simulation                                    *
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

#include "pcr.h"
#include "submat.h"
#include "sequence.h"

#include "exonerate_util.h"
#include <stdlib.h> /* For free() */
#include <string.h> /* For strlen() */

/**/

static PCR_Match *PCR_Match_allocate(PCR *pcr){
    register PCR_Match *pcr_match;
    if(pcr->match_recycle){
        pcr_match = pcr->match_recycle;
        /* pcr_match->pcr_probe is used as the next pointer */
        pcr->match_recycle = (PCR_Match*)pcr_match->pcr_probe;
    } else {
        pcr_match = (PCR_Match *)malloc(sizeof(PCR_Match));
        }
    return pcr_match;
    }

static PCR_Match *PCR_Match_create(PCR_Probe *pcr_probe, gint position,
                                   gint mismatch){
    register PCR *pcr = pcr_probe->pcr_primer->pcr_experiment->pcr;
    register PCR_Match *pcr_match = PCR_Match_allocate(pcr);
    pcr_match->pcr_probe = pcr_probe;
    pcr_match->position = position;
    pcr_match->mismatch = mismatch;
    return pcr_match;
    }

static void PCR_Match_destroy(PCR_Match *pcr_match){
    register PCR *pcr
        = pcr_match->pcr_probe->pcr_primer->pcr_experiment->pcr;
    pcr_match->pcr_probe = (PCR_Probe*)pcr->match_recycle;
    pcr->match_recycle = pcr_match;
    return;
    }

/**/

static PCR_Probe *PCR_Probe_create(PCR_Primer *pcr_primer,
       Sequence_Strand strand, gchar *probe, gint mismatch){
    register PCR *pcr = pcr_primer->pcr_experiment->pcr;
    PCR_Probe *pcr_probe = (PCR_Probe *)malloc(sizeof(PCR_Probe));
    pcr_probe->pcr_primer = pcr_primer;
    pcr_probe->strand = strand;
    pcr_probe->mismatch = mismatch;
    return pcr_probe;
    }

static void PCR_Probe_register_hit(PCR_Probe *pcr_probe,
                                   Sequence *sequence, guint seq_pos){
    register PCR_Match *pcr_match, *prev_match;
    register PCR_Primer *pcr_primer = pcr_probe->pcr_primer;
    register PCR_Experiment *pcr_experiment
                           = pcr_primer->pcr_experiment;
    register PCR *pcr = pcr_experiment->pcr;
    register gint i, product_length, mismatch = pcr_probe->mismatch;
    register SListNode *node;
    register gint match_start;
    if(pcr_probe->strand == Sequence_Strand_FORWARD){
        match_start = seq_pos - pcr_primer->probe_len + 1;
    } else {
        match_start = seq_pos - pcr_primer->length    + 1;
        }
    if(match_start < 0)
        return; /* Too close to sequence start */
    if((match_start + pcr_primer->length) > sequence->len)
        return; /* Too close to sequence end */
    /* Count number of mismatches in rest of primer */
    g_assert(pcr_primer);
    if(pcr_probe->strand == Sequence_Strand_FORWARD){
        for(i = pcr_primer->probe_len; i < pcr_primer->length; i++){
            if(pcr_primer->forward[i]
            != Sequence_get_symbol(sequence, match_start+i)){
                mismatch++;
                if(mismatch > pcr->mismatch_threshold)
                    return; /* Match failed to extend */
                }
            }
    } else {
        g_assert(pcr_primer);
        for(i = 0; i < (pcr_primer->length-pcr_primer->probe_len); i++){
            if(pcr_primer->revcomp[i]
            != Sequence_get_symbol(sequence, match_start+i)){
                mismatch++;
                if(mismatch > pcr->mismatch_threshold)
                    return; /* Match failed to extend */
                }
            }
        }
    /* Clear any matches longer in range */
    while(!SList_isempty(pcr_experiment->match_list)){
        prev_match = SList_first(pcr_experiment->match_list);
        product_length = match_start
                       - prev_match->position
                       + pcr_probe->pcr_primer->length;
        if(product_length <= pcr_experiment->max_product_len)
            break;
        prev_match = SList_pop(pcr_experiment->match_list);
        PCR_Match_destroy(prev_match);
        }
    /* Create a match for this probe */
    pcr_match = PCR_Match_create(pcr_probe, match_start, mismatch);
    /* Report hit with any previous matches in range */
    SList_for_each(pcr_experiment->match_list, node){
        prev_match = node->data;
        product_length = match_start
                       - prev_match->position
                       + pcr_probe->pcr_primer->length;
        if(product_length < pcr_experiment->min_product_len)
            break;
        /* Primers must be on different strands */
        if(prev_match->pcr_probe->strand
        != pcr_match->pcr_probe->strand)
            /* Primers must be facing each other */
            if((prev_match->pcr_probe->strand
                == Sequence_Strand_FORWARD)
            && (pcr_match->pcr_probe->strand
                == Sequence_Strand_REVCOMP))
                pcr->report_func(sequence,
                                 prev_match, pcr_match, product_length,
                                 pcr->user_data);
        }
    /* Record this match */
    SList_queue(pcr_experiment->match_list, pcr_match);
    return;
    }

/**/

static PCR_Sensor *PCR_Sensor_create(PCR *pcr){
    PCR_Sensor *pcr_sensor = (PCR_Sensor *)malloc(sizeof(PCR_Sensor));
    pcr_sensor->owned_probe_list = SList_create(pcr->slist_set);
    pcr_sensor->borrowed_sensor_list
                                    = SList_create(pcr->slist_set);
    pcr->sensor_count++;
    return pcr_sensor;
    }

static void PCR_Sensor_destroy(PCR_Sensor *pcr_sensor, PCR *pcr){
    register PCR_Probe *pcr_probe;
    if(pcr_sensor->owned_probe_list){ /* Maybe removed in merge */
        while(!SList_isempty(pcr_sensor->owned_probe_list)){
            pcr_probe = SList_pop(pcr_sensor->owned_probe_list);
            free(pcr_probe);
            }
        SList_destroy(pcr_sensor->owned_probe_list);
        }
    SList_destroy(pcr_sensor->borrowed_sensor_list);
    free(pcr_sensor);
    pcr->sensor_count--;
    return;
    }

static gpointer PCR_FSM_merge_PCR_Sensor(gpointer a, gpointer b,
                                           gpointer user_data){
    register PCR *pcr = user_data;
    register PCR_Sensor *ps_a = a, *ps_b = b;
    g_assert(ps_a);
    g_assert(ps_b);
    g_assert(SList_isempty(ps_a->borrowed_sensor_list));
    g_assert(SList_isempty(ps_b->borrowed_sensor_list));
    ps_a->owned_probe_list = SList_join(ps_a->owned_probe_list,
                                        ps_b->owned_probe_list);
    SList_destroy(ps_b->owned_probe_list);
    ps_b->owned_probe_list = NULL;
    PCR_Sensor_destroy(ps_b, pcr);
    return ps_a;
    }
/* Merge two PCR_Sensors which have the same sequence */

static gpointer PCR_FSM_combine_PCR_Sensor(gpointer a, gpointer b,
                                             gpointer user_data){
    register PCR_Sensor *ps_a = a, *ps_b = b;
    g_assert(ps_a);
    g_assert(ps_b);
    /**/
    SList_queue(ps_a->borrowed_sensor_list, ps_b);
    return ps_a;
    }
/* Join two PCR_Sensors when one is a subsequence of the other */

/**/

static void PCR_Primer_add_probe(PCR_Primer *pcr_primer,
        Sequence_Strand strand, gchar *primer, gint mismatch){
    register PCR_Probe *pcr_probe = PCR_Probe_create(pcr_primer,
              strand, primer, mismatch);
    register PCR_Sensor *pcr_sensor = PCR_Sensor_create(
                           pcr_primer->pcr_experiment->pcr);
    SList_queue(pcr_sensor->owned_probe_list, pcr_probe);
    g_ptr_array_add(pcr_primer->probe_list, pcr_probe);
    FSM_add(pcr_primer->pcr_experiment->pcr->fsm, primer,
            pcr_primer->probe_len, pcr_sensor);
    return;
    }

static gboolean PCR_Primer_Wordhood_traverse_func(gchar *word,
                    gint score, gpointer user_data){
    register PCR_Primer *pcr_primer = user_data;
    register gint mismatch = pcr_primer->probe_len-score;
    PCR_Primer_add_probe(pcr_primer, Sequence_Strand_FORWARD,
        word, mismatch);
    Sequence_revcomp_in_place(word, pcr_primer->probe_len);
    PCR_Primer_add_probe(pcr_primer, Sequence_Strand_REVCOMP,
        word, mismatch);
    Sequence_revcomp_in_place(word, pcr_primer->probe_len);
    return FALSE;
    }

static PCR_Primer *PCR_Primer_create(PCR_Experiment *pcr_experiment,
                                     gchar *primer){
    PCR_Primer *pcr_primer = (PCR_Primer*)malloc(sizeof(PCR_Primer)* 1);
    g_assert(pcr_experiment);
    g_assert(primer);
    pcr_primer->pcr_experiment = pcr_experiment;
    pcr_primer->length = strlen(primer);
    if((pcr_experiment->pcr->seed_length)
    && (pcr_primer->length > pcr_experiment->pcr->seed_length))
        pcr_primer->probe_len = pcr_experiment->pcr->seed_length;
    else
        pcr_primer->probe_len = pcr_primer->length;
    pcr_primer->forward = strdup(primer);
    pcr_primer->revcomp = strdup(primer);
    Sequence_revcomp_in_place(pcr_primer->revcomp, pcr_primer->length);
    strup(pcr_primer->forward);
    strup(pcr_primer->revcomp);
    pcr_primer->probe_list = g_ptr_array_new();
    WordHood_traverse(pcr_experiment->pcr->wordhood,
                      PCR_Primer_Wordhood_traverse_func,
                      pcr_primer->forward,
                      pcr_primer->probe_len,
                      pcr_primer);
    return pcr_primer;
    }

static void PCR_Primer_destroy(PCR_Primer *pcr_primer){
    register gint i;
    for(i = 0; i < pcr_primer->probe_list->len; i++)
        free(pcr_primer->probe_list->pdata[i]);
    free(pcr_primer->forward);
    free(pcr_primer->revcomp);
    g_ptr_array_free(pcr_primer->probe_list, TRUE);
    g_free(pcr_primer);
    }

static gsize PCR_Primer_memory_usage(PCR_Primer *pcr_primer){
    register gsize primer_size = sizeof(gchar) * pcr_primer->length;
    return sizeof(PCR_Primer)
         + (primer_size * 2)
         + (pcr_primer->probe_list->len
            * (sizeof(PCR_Probe) + sizeof(gpointer)));
    }

/**/

static PCR_Experiment *PCR_Experiment_create(PCR *pcr, gchar *id,
                 gchar *primer_a, gchar *primer_b,
                 gint min_product_len, gint max_product_len){
    PCR_Experiment *pcr_experiment = (PCR_Experiment*)malloc(sizeof(PCR_Experiment)* 1);
    /*
    g_message("Adding [%s][%s][%s][%d][%d]",
            id, primer_a, primer_b, min_product_len, max_product_len);
    */
    g_assert(id);
    g_assert(primer_a);
    g_assert(primer_b);
    g_assert(min_product_len > 0);
    g_assert(max_product_len >= min_product_len);
    pcr_experiment->pcr = pcr;
    pcr_experiment->id = strdup(id);
    pcr_experiment->primer_a = PCR_Primer_create(pcr_experiment,
                                                 primer_a);
    pcr_experiment->primer_b = PCR_Primer_create(pcr_experiment,
                                                 primer_b);
    pcr_experiment->min_product_len = min_product_len;
    pcr_experiment->max_product_len = max_product_len;
    pcr_experiment->match_list = SList_create(pcr->slist_set);
    pcr_experiment->product_count = 0;
    return pcr_experiment;
    }

static void PCR_Experiment_clear(PCR_Experiment *pcr_experiment){
    register PCR_Match *pcr_match;
    while(!SList_isempty(pcr_experiment->match_list)){
        pcr_match = SList_pop(pcr_experiment->match_list);
        PCR_Match_destroy(pcr_match);
        }
    }

static void PCR_Experiment_destroy(PCR_Experiment *pcr_experiment){
    g_assert(pcr_experiment);
    PCR_Experiment_clear(pcr_experiment);
    PCR_Primer_destroy(pcr_experiment->primer_a);
    PCR_Primer_destroy(pcr_experiment->primer_b);
    SList_destroy(pcr_experiment->match_list);
    free(pcr_experiment->id);
    g_free(pcr_experiment);
    return;
    }

static gsize PCR_Experiment_memory_usage(
             PCR_Experiment *pcr_experiment){
    return sizeof(PCR_Experiment)
         + sizeof(gpointer) /* for pcr->experment_list */
         + (sizeof(gchar) * strlen(pcr_experiment->id))
         + PCR_Primer_memory_usage(pcr_experiment->primer_a)
         + PCR_Primer_memory_usage(pcr_experiment->primer_b);
    }

/**/

PCR *PCR_create(PCR_ReportFunc report_func, gpointer user_data,
                gint mismatch_threshold, gint seed_length){
    PCR *pcr = (PCR*)malloc(sizeof(PCR)* 1);
    register Submat *submat = Submat_create("iupac-identity");
    register WordHood_Alphabet *wha;
    pcr->slist_set = SListSet_create();
    pcr->fsm = FSM_create("ACGT", PCR_FSM_merge_PCR_Sensor,
                                  PCR_FSM_combine_PCR_Sensor, pcr);
    pcr->experiment_list = g_ptr_array_new();
    pcr->mismatch_threshold = mismatch_threshold;
    pcr->seed_length = seed_length;
    pcr->report_func = report_func;
    pcr->is_prepared = FALSE;
    wha = WordHood_Alphabet_create_from_Submat(
                            "GARTKWDCSMVYBHN", "ACGT", submat, FALSE);
    pcr->wordhood = WordHood_create(wha, pcr->mismatch_threshold, TRUE);
    WordHood_Alphabet_destroy(wha);
    pcr->experiment_memory_usage = 0;
    pcr->sensor_count = 0;
    pcr->user_data = user_data;
    Submat_destroy(submat);
    pcr->match_recycle = NULL;
    return pcr;
    }

void PCR_destroy(PCR *pcr){
    register gint i;
    for(i = 0; i < pcr->experiment_list->len; i++)
        PCR_Experiment_destroy(pcr->experiment_list->pdata[i]);
    g_ptr_array_free(pcr->experiment_list, TRUE);
    FSM_destroy(pcr->fsm);
    WordHood_destroy(pcr->wordhood);
    SListSet_destroy(pcr->slist_set);
    free(pcr);
    return;
    }

/**/

static gsize PCR_memory_usage(PCR *pcr){
    return sizeof(PCR)
         + pcr->experiment_memory_usage
         + SListSet_memory_usage(pcr->slist_set)
         + FSM_memory_usage(pcr->fsm)
         + (sizeof(PCR_Sensor) * pcr->sensor_count);
    }

gsize PCR_add_experiment(PCR *pcr, gchar *id,
                        gchar *primer_a, gchar *primer_b,
                        gint min_product_len, gint max_product_len){
    register PCR_Experiment *pcr_experiment;
    g_assert(pcr);
    g_assert(!pcr->is_prepared);
    pcr_experiment = PCR_Experiment_create(pcr, id, primer_a, primer_b,
                                   min_product_len, max_product_len);
    g_ptr_array_add(pcr->experiment_list, pcr_experiment);
    pcr->experiment_memory_usage
         += PCR_Experiment_memory_usage(pcr_experiment);
    return PCR_memory_usage(pcr);
    }

void PCR_prepare(PCR *pcr){
    g_assert(pcr);
    g_assert(!pcr->is_prepared);
    FSM_compile(pcr->fsm);
    pcr->is_prepared = TRUE;
    return;
    }

static void PCR_Sensor_register_hit_recur(PCR_Sensor *pcr_sensor,
                                 Sequence *sequence, guint seq_pos){
    register SListNode *node;
    register PCR_Probe *pcr_probe;
    register PCR_Sensor *borrowed_sensor;
    SList_for_each(pcr_sensor->owned_probe_list, node){
        pcr_probe = node->data;
        PCR_Probe_register_hit(pcr_probe, sequence, seq_pos);
        }
    SList_for_each(pcr_sensor->borrowed_sensor_list, node){
        borrowed_sensor = node->data;
        PCR_Sensor_register_hit_recur(borrowed_sensor,
                                        sequence, seq_pos);
        }
    return;
    }

static void PCR_FSM_traverse_func(guint seq_pos,
                                  gpointer node_data,
                                  gpointer user_data){
    register PCR_Sensor *pcr_sensor = node_data;
    register Sequence *sequence = user_data;
    PCR_Sensor_register_hit_recur(pcr_sensor, sequence, seq_pos);
    return;
    }

void PCR_simulate(PCR *pcr, Sequence *sequence){
    register gint i;
    register gchar *str = Sequence_get_str(sequence);
    g_assert(pcr->is_prepared);
    FSM_traverse(pcr->fsm, str,
                 PCR_FSM_traverse_func, sequence);
    g_free(str);
    for(i = 0; i < pcr->experiment_list->len; i++)
        PCR_Experiment_clear(pcr->experiment_list->pdata[i]);
    return;
    }

/**/

/* FIXME: Add wordhood_abandon and match_abandon thresholds.
 */

