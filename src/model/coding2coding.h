/****************************************************************\
*                                                                *
*  Coding DNA comparison model                                   *
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

#ifndef INCLUDED_CODING2CODING_H
#define INCLUDED_CODING2CODING_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "c4.h"
#include "sequence.h"
#include "affine.h"

typedef struct {
    Affine_Data  ad; /* inherit */
} Coding2Coding_Data;

/**/

#define Coding2Coding_Data_get_submat(c2dd) \
        Affine_Data_get_protein_submat(&((c2dd)->ad))
#define Coding2Coding_Data_get_translate(c2dd) \
        Affine_Data_get_translate(&((c2dd)->ad))
#define Coding2Coding_Data_get_Frameshift_Data(c2dd) \
        Affine_Data_get_Frameshift_Data(&((c2dd)->ad))
#define Coding2Coding_Data_get_Intron_Data(c2dd) \
        Affine_Data_get_Intron_Data(&((c2dd)->ad))

void Coding2Coding_Data_init(Coding2Coding_Data *c2cd,
                             Sequence *query, Sequence *target);
Coding2Coding_Data *Coding2Coding_Data_create(Sequence *query,
                                              Sequence *target);
void Coding2Coding_Data_clear(Coding2Coding_Data *c2cd);
void Coding2Coding_Data_destroy(Coding2Coding_Data *c2cd);

/**/

C4_Model *Coding2Coding_create(void);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* INCLUDED_CODING2CODING_H */

