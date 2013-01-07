/****************************************************************\
*                                                                *
*  Comparison : A module for pairwise sequence comparisons       *
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

#ifndef INCLUDED_COMPARISON_H
#define INCLUDED_COMPARISON_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <glib.h>
#include "argument.h"
#include "alphabet.h"
#include "sequence.h"
#include "translate.h"
#include "submat.h"
#include "hspset.h"
#include "threadref.h"

typedef struct Comparison_Param {
                 ThreadRef *thread_ref;
             Alphabet_Type  query_type;
             Alphabet_Type  target_type;
                 HSP_Param *dna_hsp_param;
                 HSP_Param *protein_hsp_param;
                 HSP_Param *codon_hsp_param;
   struct Comparison_Param *mirror;
} Comparison_Param;

Comparison_Param *Comparison_Param_create(Alphabet_Type query_type,
                                          Alphabet_Type target_type,
                                          HSP_Param *dna_hsp_param,
                                          HSP_Param *protein_hsp_param,
                                          HSP_Param *codon_hsp_param);
void Comparison_Param_destroy(Comparison_Param *param);
Comparison_Param *Comparison_Param_share(Comparison_Param *param);
Comparison_Param *Comparison_Param_swap(Comparison_Param *param);

/**/

typedef struct {
           ThreadRef *thread_ref;
    Comparison_Param *param;
            Sequence *query;
            Sequence *target;
              HSPset *dna_hspset;
              HSPset *protein_hspset;
              HSPset *codon_hspset;
} Comparison;

Comparison *Comparison_create(Comparison_Param *param,
                              Sequence *query, Sequence *target);
Comparison *Comparison_share(Comparison *comparison);
void Comparison_destroy(Comparison *comparison);
void Comparison_print(Comparison *comparison);
HSPset_ArgumentSet *Comparison_Param_get_HSPSet_Argument_Set(
                    Comparison_Param *param);

gboolean  Comparison_has_hsps(Comparison *comparison);
    void  Comparison_finalise(Comparison *comparison);
    void  Comparison_swap(Comparison *comparison);
    void  Comparison_revcomp(Comparison *comparison);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* INCLUDED_COMPARISON_H */

