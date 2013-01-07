/****************************************************************\
*                                                                *
*  Genome <-> Genome comparison model                            *
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

#ifndef INCLUDED_GENOME2GENOME_H
#define INCLUDED_GENOME2GENOME_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "c4.h"
#include "sequence.h"
#include "affine.h"
#include "cdna2genome.h"

typedef struct {
    CDNA2Genome_Data cd2gd;
} Genome2Genome_Data;

#define Genome2Genome_Data_get_dna_submat(g2gd) \
     CDNA2Genome_Data_get_dna_submat(&((g2gd)->cd2gd))
#define Genome2Genome_Data_get_protein_submat(g2gd) \
     CDNA2Genome_Data_get_protein_submat(&((g2gd)->cd2gd))
#define Genome2Genome_Data_get_translate(g2gd) \
     CDNA2Genome_Data_get_translate(&((g2gd)->cd2gd))
#define Genome2Genome_Data_get_Intron_Data(g2gd) \
     CDNA2Coding_Data_get_Intron_Data(&((g2gd)->cd2gd))

/**/

Genome2Genome_Data *Genome2Genome_Data_create(Sequence *query,
                                              Sequence *target);
void Genome2Genome_Data_destroy(Genome2Genome_Data *p2gd);

/**/

C4_Model *Genome2Genome_create(void);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* INCLUDED_GENOME2GENOME_H */

