/****************************************************************\
*                                                                *
*  Protein <-> DNA comparison model                              *
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

#ifndef INCLUDED_PROTEIN2DNA_H
#define INCLUDED_PROTEIN2DNA_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "c4.h"
#include "sequence.h"
#include "affine.h"
#include "frameshift.h"

typedef struct {
    Affine_Data  ad; /* inherit */
} Protein2DNA_Data;

/**/

#define Protein2DNA_Data_get_submat(p2dd) \
        Affine_Data_get_protein_submat(&((p2dd)->ad))
#define Protein2DNA_Data_get_translate(p2dd) \
        Affine_Data_get_translate(&((p2dd)->ad))
#define Protein2DNA_Data_get_Intron_Data(p2dd) \
        Affine_Data_get_Intron_Data(&((p2dd)->ad))
#define Protein2DNA_Data_get_Frameshift_Data(p2dd) \
        Affine_Data_get_Frameshift_Data(&((p2dd)->ad))
#define Protein2DNA_Data_get_query(p2dd) \
        Affine_Data_get_query(&((p2dd)->ad))
#define Protein2DNA_Data_get_target(p2dd) \
        Affine_Data_get_target(&((p2dd)->ad))

void Protein2DNA_Data_init(Protein2DNA_Data *p2dd,
                           Sequence *query, Sequence *target);
Protein2DNA_Data *Protein2DNA_Data_create(Sequence *query,
                                          Sequence *target);
void Protein2DNA_Data_clear(Protein2DNA_Data *p2dd);
void Protein2DNA_Data_destroy(Protein2DNA_Data *p2dd);

/**/

C4_Model *Protein2DNA_create(Affine_Model_Type type);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* INCLUDED_PROTEIN2DNA_H */

