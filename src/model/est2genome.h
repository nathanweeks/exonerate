/****************************************************************\
*                                                                *
*  Module for EST <-> genome alignments                          *
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

#ifndef INCLUDED_EST2GENOME_H
#define INCLUDED_EST2GENOME_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <ctype.h> /* For tolower() */

#include "c4.h"
#include "sequence.h"
#include "heuristic.h"
#include "splice.h"
#include "argument.h"
#include "affine.h"

typedef struct {
    Affine_Data ad;
} EST2Genome_Data;

#define EST2Genome_Data_get_submat(e2gd) \
        Affine_Data_get_dna_submat(&((e2gd)->ad))
#define EST2Genome_Data_get_Intron_Data(e2gd) \
        Affine_Data_get_Intron_Data(&((e2gd)->ad))

void EST2Genome_Data_init(EST2Genome_Data *e2gd,
                          Sequence *query, Sequence *target);
EST2Genome_Data *EST2Genome_Data_create(Sequence *query,
                                        Sequence *target);


void EST2Genome_Data_clear(EST2Genome_Data *e2gd);
void EST2Genome_Data_destroy(EST2Genome_Data *e2gd);

C4_Model *EST2Genome_create(void);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* INCLUDED_EST2GENOME_H */

