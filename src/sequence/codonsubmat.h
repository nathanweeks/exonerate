/****************************************************************\
*                                                                *
*  Codon Substitution Matrix Object                              *
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

#ifndef INCLUDED_CODONSUBMAT_H
#define INCLUDED_CODONSUBMAT_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <limits.h>

#include <glib.h>
#include <limits.h>
#include <submat.h>

#ifndef ALPHABETSIZE
#define ALPHABETSIZE (1 << CHAR_BIT)
#endif /* ALPHABETSIZE */

#define CODON_ALPHABETSIZE_ACGT  (4*4*4)
#define CODON_ALPHABETSIZE_ACGTN (5*5*5)

typedef struct {
         gint  ref_count;
       guchar  base_index[ALPHABETSIZE];
         gint  codon_index[5][5][5];
         gint  codon_submat[CODON_ALPHABETSIZE_ACGTN]
                           [CODON_ALPHABETSIZE_ACGTN];
} CodonSubmat;

CodonSubmat *CodonSubmat_create(void);
       void  CodonSubmat_destroy(CodonSubmat *cs);
CodonSubmat *CodonSubmat_share(CodonSubmat *cs);
       gint  CodonSubmat_max_score(CodonSubmat *cs);
       void  CodonSubmat_add_nucleic(CodonSubmat *cs, Submat *nucleic);

#define CodonSubmat_get_base(codonsubmat, base) \
       ((codonsubmat)->base_index[(base)])

#define CodonSubmat_get_codon(codonsubmat, codon)           \
       ((codonsubmat)->codon_index                          \
            [(CodonSubmat_get_base(codonsubmat, codon[0]))] \
            [(CodonSubmat_get_base(codonsubmat, codon[1]))] \
            [(CodonSubmat_get_base(codonsubmat, codon[2]))] )

#define CodonSubmat_get_codon_base(codonsubmat, base_a, base_b, base_c) \
       ((codonsubmat)->codon_index                                      \
            [(CodonSubmat_get_base(codonsubmat, base_a))]               \
            [(CodonSubmat_get_base(codonsubmat, base_b))]               \
            [(CodonSubmat_get_base(codonsubmat, base_c))]               )

#define CodonSubmat_lookup_codon(codonsubmat, codon_a, codon_b)       \
       ((codonsubmat)->codon_submat[(codon_a)][(codon_b)])

#define CodonSubmat_lookup(codonsubmat, codon_a, codon_b)       \
       ((codonsubmat)->codon_submat                             \
            [(CodonSubmat_get_codon((codonsubmat), (codon_a)))] \
            [(CodonSubmat_get_codon((codonsubmat), (codon_b)))] )

#define CodonSubmat_lookup_base(codonsubmat, base_a1, base_a2, base_a3, \
                                             base_b1, base_b2, base_b3) \
       ((codonsubmat)->codon_submat                                     \
            [(CodonSubmat_get_codon_base((codonsubmat),                 \
                                    (base_a1), (base_a2), (base_a3)))]  \
            [(CodonSubmat_get_codon_base((codonsubmat),                 \
                                    (base_b1), (base_b2), (base_b3)))])

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* INCLUDED_SUBMAT_H */

