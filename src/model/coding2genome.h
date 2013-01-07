/****************************************************************\
*                                                                *
*  Coding <-> Genome comparison model                            *
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

#ifndef INCLUDED_CODING2GENOME_H
#define INCLUDED_CODING2GENOME_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "c4.h"
#include "sequence.h"
#include "affine.h"
#include "coding2coding.h"

typedef struct {
    Coding2Coding_Data c2cd;
} Coding2Genome_Data;

#define Coding2Genome_Data_get_submat(g2gd) \
     Coding2Coding_Data_get_submat(&((g2gd)->c2cd))
#define Coding2Genome_Data_get_translate(g2gd) \
     Coding2Coding_Data_get_translate(&((g2gd)->c2cd))
#define Coding2Genome_Data_get_Intron_Data(g2gd) \
     Coding2Coding_Data_get_Intron_Data(&((g2gd)->c2cd))

/**/

void Coding2Genome_Data_init(Coding2Genome_Data *c2gd,
                             Sequence *query, Sequence *target);
Coding2Genome_Data *Coding2Genome_Data_create(Sequence *query,
                                              Sequence *target);
void Coding2Genome_Data_clear(Coding2Genome_Data *c2gd);
void Coding2Genome_Data_destroy(Coding2Genome_Data *c2gd);


/**/

C4_Model *Coding2Genome_create(void);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* INCLUDED_CODING2GENOME_H */

