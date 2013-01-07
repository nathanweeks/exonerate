/****************************************************************\
*                                                                *
*  Module for various affine gapped models                       *
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

#ifndef INCLUDED_AFFINE_H
#define INCLUDED_AFFINE_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "c4.h"
#include "ungapped.h"
#include "argument.h"

typedef struct {
    C4_Score gap_open;
    C4_Score gap_extend;
    C4_Score codon_gap_open;
    C4_Score codon_gap_extend;
} Affine_ArgumentSet;

Affine_ArgumentSet *Affine_ArgumentSet_create(Argument *arg);

typedef enum {
      Affine_Model_Type_GLOBAL,
      Affine_Model_Type_BESTFIT,
      Affine_Model_Type_LOCAL,
      Affine_Model_Type_OVERLAP,
      Affine_Model_Type_UNKNOWN
} Affine_Model_Type;

gchar *Affine_Model_Type_get_name(Affine_Model_Type type);

typedef struct {
         Ungapped_Data  ud; /* inherit */
    Affine_ArgumentSet *aas;
} Affine_Data;

#define Affine_Data_get_query(ad) Ungapped_Data_get_query(&((ad)->ud))
#define Affine_Data_get_target(ad) Ungapped_Data_get_target(&((ad)->ud))
#define Affine_Data_get_dna_submat(ad) \
        Ungapped_Data_get_dna_submat(&((ad)->ud))
#define Affine_Data_get_protein_submat(ad) \
        Ungapped_Data_get_protein_submat(&((ad)->ud))
#define Affine_Data_get_translate(ad) \
        Ungapped_Data_get_translate(&((ad)->ud))
#define Affine_Data_get_Intron_Data(ad) \
        Ungapped_Data_get_Intron_Data(&((ad)->ud))
#define Affine_Data_get_Frameshift_Data(ad) \
        Ungapped_Data_get_Frameshift_Data(&((ad)->ud))

void Affine_Data_init(Affine_Data *ad,
                      Sequence *query, Sequence *target,
                      gboolean translate_both);
Affine_Data *Affine_Data_create(Sequence *query, Sequence *target,
                                gboolean translate_both);
void Affine_Data_clear(Affine_Data *ad);
void Affine_Data_destroy(Affine_Data *ad);

C4_Model *Affine_create(Affine_Model_Type type,
                        Alphabet_Type query_type,
                        Alphabet_Type target_type,
                        gboolean translate_both);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* INCLUDED_AFFINE_H */

