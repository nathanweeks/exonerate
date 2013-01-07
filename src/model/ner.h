/****************************************************************\
*                                                                *
*  Model for alignments with non-equivalenced regions            *
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

#ifndef INCLUDED_NER_H

#ifdef __cplusplus
extern "C"  {
#endif /* __cplusplus */

#include "c4.h"
#include "affine.h"
#include "argument.h"

typedef struct {
        gint min_ner;
        gint max_ner;
    C4_Score ner_open_penalty;
} NER_ArgumentSet;

NER_ArgumentSet *NER_ArgumentSet_create(Argument *arg);

typedef struct {
    Affine_Data  ad; /* inherit */
NER_ArgumentSet *nas;
           gint  curr_ner_id;
           gint  curr_ner_query_start;
           gint  curr_ner_target_start;
} NER_Data;

#define NER_Data_get_query(nd) Affine_Data_get_query(&((nd)->ad))
#define NER_Data_get_target(nd) Affine_Data_get_target(&((nd)->ad))
#define NER_Data_get_dna_submat(nd) Affine_Data_get_dna_submat(&((nd)->ad))
#define NER_Data_get_protein_submat(nd) \
        Affine_Data_get_protein_submat(&((nd)->ad))

NER_Data *NER_Data_create(Sequence *query, Sequence *target);
void NER_Data_destroy(NER_Data *gnd);

C4_Model *NER_create(Alphabet_Type query_type,
                     Alphabet_Type target_type);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* INCLUDED_NER_H */

