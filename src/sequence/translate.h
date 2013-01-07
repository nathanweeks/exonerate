/****************************************************************\
*                                                                *
*  Nucleotide Translation Code                                   *
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

#ifndef INCLUDED_TRANSLATE_H
#define INCLUDED_TRANSLATE_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/**/

#include <glib.h>
#include "argument.h"

typedef struct {
    gchar *genetic_code;
} Translate_ArgumentSet;

Translate_ArgumentSet *Translate_ArgumentSet_create(Argument *arg);

#define Translate_AA_SET_SIZE          40
#define Translate_NT_SET_SIZE       (1<<4)
#define Translate_PIMA_SET_SIZE        18
#define Translate_ALPHABET_SIZE     (1<<8)
#define Translate_TRANSLATION_SIZE   4096

typedef struct {
        guint  ref_count;
       guchar *nt;
       guchar *aa;
       guchar *code;
       guchar  nt2d[Translate_ALPHABET_SIZE];
       guchar  aa2d[Translate_ALPHABET_SIZE];
         gint  aamask[Translate_AA_SET_SIZE];
       guchar  trans[Translate_TRANSLATION_SIZE];
    GPtrArray *revtrans[Translate_ALPHABET_SIZE];
} Translate;

Translate *Translate_create(gboolean use_pima);
Translate *Translate_share(Translate *t);
     void  Translate_destroy(Translate *t);
     gint  Translate_sequence(Translate *t, gchar *dna,
                  gint dna_length, gint frame, gchar *aaseq,
                  guchar *filter);
/* Sequence_translate frame is  1, 2, 3 for forward,
 *                          or -1,-2,-3 for reverse.
 * Sufficient space for translated sequence must provided in aaseq
 * Translated seq is NULL terminated.
 * Length of translated sequence is returned
 */

typedef void (*Translate_reverse_func)(gchar *dna, gint length,
                                       gpointer user_data);
void Translate_reverse(Translate *t, gchar *aaseq, gint length,
                       Translate_reverse_func trf, gpointer user_data);

#define Translate_codon(t, codon)                                     \
    ((gchar)((t)->aa[(t)->trans[((t)->nt2d[(guchar)(codon)[0]])       \
                               |((t)->nt2d[(guchar)(codon)[1]]<<4)    \
                               |((t)->nt2d[(guchar)(codon)[2]]<<8)]]))

#define Translate_base(t, a, b, c)                          \
         ((t)->aa[(t)->trans[((t)->nt2d[(guchar)(a)])       \
                            |((t)->nt2d[(guchar)(b)]<<4)    \
                            |((t)->nt2d[(guchar)(c)]<<8)]])

/**/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* INCLUDED_TRANSLATE_H */

