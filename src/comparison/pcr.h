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

#ifndef INCLUDED_PCR_H
#define INCLUDED_PCR_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <glib.h>
#include "fsm.h"
#include "sequence.h"
#include "wordhood.h"
#include "slist.h"

typedef struct {
    SList *owned_probe_list;
    SList *borrowed_sensor_list;
} PCR_Sensor;
/* pcr->fsm nodes contain PCR_Sensor objects.
 *
 * This is to allow both cases where:
 *  (a) Multiple probes are present with the same sequence
 *  (b) One probe is a subsequence of another
 */

typedef struct {
    struct PCR_Primer *pcr_primer;
      Sequence_Strand  strand;
                 gint  mismatch;
} PCR_Probe;
/* A derived primer to facilitate approximate primer matching */

typedef struct {
    PCR_Probe *pcr_probe;
         gint  position; /* Position at probe end */
         gint  mismatch;
} PCR_Match;
/* A primer annealed to a position on a sequence */

typedef struct PCR_Primer {
    struct PCR_Experiment *pcr_experiment;
                     gint  length;
                     gint  probe_len;
                    gchar *forward;
                    gchar *revcomp;
                GPtrArray *probe_list;
} PCR_Primer;

typedef struct PCR_Experiment {
    struct PCR *pcr;
         gchar *id;
    PCR_Primer *primer_a;
    PCR_Primer *primer_b;
          gint  min_product_len;
          gint  max_product_len;
         SList *match_list;
          gint  product_count;
} PCR_Experiment;
/* A pair of primers to be simulated in a PCR reaction */

typedef gboolean (*PCR_ReportFunc)(Sequence *sequence,
                PCR_Match *match_a, PCR_Match *match_b,
                gint product_length, gpointer user_data);
/* Return TRUE to stop the simulation */

typedef struct PCR {
          SListSet *slist_set; /* Lists of PCR_Probes */
               FSM *fsm;
         GPtrArray *experiment_list;
              gint  mismatch_threshold;
              gint  seed_length;
          WordHood *wordhood;
    PCR_ReportFunc  report_func;
          gboolean  is_prepared;
             gsize  experiment_memory_usage;
              gint  sensor_count;
          gpointer  user_data;
         PCR_Match *match_recycle;
} PCR;

 PCR *PCR_create(PCR_ReportFunc report_func,
                 gpointer user_data, gint mismatch_threshold,
                 gint seed_length);
void  PCR_destroy(PCR *pcr);

gsize PCR_add_experiment(PCR *pcr, gchar *id,
                         gchar *primer_a, gchar *primer_b,
                         gint min_product_len, gint max_product_len);
/* Returns the memory currently being used by PCR */

void PCR_prepare(PCR *pcr);
/* Called after last call to PCR_add_primer_pair() */

void PCR_simulate(PCR *pcr, Sequence *sequence);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* INCLUDED_IPCRESS_H */

