/****************************************************************\
*                                                                *
*  Module of splice site and intron modelling                    *
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

#ifndef INCLUDED_INTRON_H
#define INCLUDED_INTRON_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "c4.h"
#include "splice.h"
#include "argument.h"
#include "phase.h"

typedef struct {
                  gint  min_intron;
                  gint  max_intron;
    SplicePredictorSet *sps;
              C4_Score  intron_open_penalty;
} Intron_ArgumentSet;

Intron_ArgumentSet *Intron_ArgumentSet_create(Argument *arg);

typedef struct {
                    gint  curr_intron_start;
    SplicePrediction_Set *sps;
} Intron_ChainData;

/* FIXME: can remove chain data and just use curr_{qy,tg}_intron_start */

typedef struct {
    Intron_ArgumentSet *ias;
      Intron_ChainData *query_data;
      Intron_ChainData *target_data;
                Region *curr_region;
} Intron_Data;

Intron_Data *Intron_Data_create(void);
       void  Intron_Data_destroy(Intron_Data *id);

C4_Model *Intron_create(gchar *suffix, gboolean on_query,
                        gboolean on_target, gboolean is_forward);

void Intron_ChainData_init_splice_prediction(Intron_ChainData *icd,
        Sequence *s, SplicePredictorSet *sps, SpliceType type,
        gboolean use_single);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* INCLUDED_INTRON_H */

