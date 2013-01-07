/****************************************************************\
*                                                                *
*  cDNA <-> genomic alignment model                              *
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

#ifndef INCLUDED_CDNA2GENOME_H
#define INCLUDED_CDNA2GENOME_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "c4.h"
#include "sequence.h"
#include "affine.h"
#include "coding2genome.h"
#include "est2genome.h"
#include "argument.h"

typedef union { /* NB. Must be a union for multiple-inheritance hack */
    Coding2Genome_Data  c2gd;
       EST2Genome_Data  e2gd;
} CDNA2Genome_Data;

#define CDNA2Genome_Data_get_dna_submat(cd2gd) \
         EST2Genome_Data_get_submat(&((cd2gd)->e2gd))
#define CDNA2Genome_Data_get_protein_submat(cd2gd) \
         Coding2Genome_Data_get_submat(&((cd2gd)->c2gd))
#define CDNA2Genome_Data_get_translate(cd2gd) \
         Coding2Genome_Data_get_translate(&((cd2gd)->c2gd))
#define CDNA2Genome_Data_get_Intron_Data(cd2gd) \
         Coding2Genome_Data_get_Intron_Data(&((cd2gd)->c2gd))

/**/

void CDNA2Genome_Data_init(CDNA2Genome_Data *cd2gd,
                           Sequence *query, Sequence *target);
CDNA2Genome_Data *CDNA2Genome_Data_create(Sequence *query,
                                          Sequence *target);
void CDNA2Genome_Data_clear(CDNA2Genome_Data *cd2gd);
void CDNA2Genome_Data_destroy(CDNA2Genome_Data *cd2gd);

/**/

C4_Model *CDNA2Genome_create(void);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* INCLUDED_CDNA2GENOME_H */

