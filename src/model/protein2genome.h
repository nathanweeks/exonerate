/****************************************************************\
*                                                                *
*  Protein <-> Genome comparison model                           *
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

#ifndef INCLUDED_PROTEIN2GENOME_H
#define INCLUDED_PROTEIN2GENOME_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "c4.h"
#include "sequence.h"
#include "affine.h"
#include "protein2dna.h"

typedef struct {
    Protein2DNA_Data p2dd;
} Protein2Genome_Data;

/**/

#define Protein2Genome_Data_get_submat(p2gd) \
        Protein2DNA_Data_get_submat(&((p2gd)->p2dd))
#define Protein2Genome_Data_get_translate(p2gd) \
        Protein2DNA_Data_get_translate(&((p2gd)->p2dd))
#define Protein2Genome_Data_get_Intron_Data(p2gd) \
        Protein2DNA_Data_get_Intron_Data(&((p2gd)->p2dd))
#define Protein2Genome_Data_get_query(p2gd) \
        Protein2DNA_Data_get_query(&((p2gd)->p2dd))
#define Protein2Genome_Data_get_target(p2gd) \
        Protein2DNA_Data_get_target(&((p2gd)->p2dd))

Protein2Genome_Data *Protein2Genome_Data_create(Sequence *query,
                                                Sequence *target);
void Protein2Genome_Data_destroy(Protein2Genome_Data *p2gd);

/**/

C4_Model *Protein2Genome_create(Affine_Model_Type type);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* INCLUDED_PROTEIN2GENOME_H */

