/****************************************************************\
*                                                                *
*  Module for various ungapped alignment models                  *
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

#ifndef INCLUDED_UNGAPPED_H
#define INCLUDED_UNGAPPED_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "c4.h"
#include "heuristic.h"
#include "submat.h"
#include "translate.h"
#include "sequence.h"
#include "alignment.h"
#include "argument.h"
#include "intron.h"
#include "frameshift.h"
#include "match.h"
#include "hspset.h"

typedef struct {
          Heuristic_Data  heuristic_data; /* Inherit */
                Sequence *query;
                Sequence *target;
       Match_ArgumentSet *mas;
                   Match *match_list[Match_Type_TOTAL];
             Intron_Data *intron_data;
         Frameshift_Data *frameshift_data;
} Ungapped_Data;

#define Ungapped_Data_get_query(ud) (ud)->query
#define Ungapped_Data_get_target(ud) (ud)->target
#define Ungapped_Data_get_dna_submat(ud) ((ud)->mas->dna_submat)
#define Ungapped_Data_get_protein_submat(ud) ((ud)->mas->protein_submat)
#define Ungapped_Data_get_translate(ud) (ud)->mas->translate
#define Ungapped_Data_get_Intron_Data(ud) (ud)->intron_data
#define Ungapped_Data_get_Frameshift_Data(ud) (ud)->frameshift_data

void Ungapped_Data_init(Ungapped_Data *ud,
                        Sequence *query, Sequence *target,
                        Match_Type match_type);
Ungapped_Data *Ungapped_Data_create(Sequence *query, Sequence *target,
                        Match_Type match_type);

void Ungapped_Data_clear(Ungapped_Data *ud);
void Ungapped_Data_destroy(Ungapped_Data *ud);

C4_Model *Ungapped_create(Match_Type match_type);

Alignment *Ungapped_Alignment_create(C4_Model *model,
                                     Ungapped_Data *ud, HSP *hsp);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* INCLUDED_UNGAPPED_H */

