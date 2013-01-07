/****************************************************************\
*                                                                *
*  Simple Alphabet Object                                        *
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

#ifndef INCLUDED_ALPHABET_H
#define INCLUDED_ALPHABET_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <stdio.h>
#include <glib.h>

#include "argument.h"

#ifndef ALPHABET_SIZE
#define ALPHABET_SIZE (CHAR_BIT << 8)
#endif /* ALPHABET_SIZE */

typedef struct {
    gboolean use_aa_tla;
} Alphabet_ArgumentSet;

Alphabet_ArgumentSet *Alphabet_ArgumentSet_create(Argument *arg);

/**/

typedef enum {
    Alphabet_Type_DNA,
    Alphabet_Type_PROTEIN,
    Alphabet_Type_UNKNOWN
} Alphabet_Type;

Alphabet_Type Alphabet_Type_guess(gchar *str);

gchar *Alphabet_Argument_parse_alphabet_type(gchar *arg_string, gpointer data);

typedef struct {
             gint  ref_count;
         gboolean  is_soft_masked;
           guchar *member;
    Alphabet_Type  type;
           guchar *is_member;
           guchar *unmasked;
           guchar *masked;
           guchar *is_valid;
           guchar *complement;
           guchar *clean;
           guchar *non_ambig;
} Alphabet;

Alphabet *Alphabet_create(Alphabet_Type type,
                          gboolean is_soft_masked);
    void  Alphabet_destroy(Alphabet *alphabet);
Alphabet *Alphabet_share(Alphabet *alphabet);

/**/

typedef enum {
    Alphabet_Filter_Type_UNMASKED,
    Alphabet_Filter_Type_MASKED,
    Alphabet_Filter_Type_IS_VALID,
    Alphabet_Filter_Type_COMPLEMENT,
    Alphabet_Filter_Type_CLEAN,
    Alphabet_Filter_Type_CLEAN_ACGTN,
    Alphabet_Filter_Type_NON_AMBIG,
    Alphabet_Filter_Type_NUMBER_OF_TYPES
} Alphabet_Filter_Type;

guchar *Alphabet_get_filter_by_type(Alphabet *alphabet,
                                    Alphabet_Filter_Type filter_type);

gchar *Alphabet_Filter_Type_get_name(Alphabet_Filter_Type filter_type);

#define Alphabet_is_masked(alphabet, symbol) \
         ((alphabet)->unmasked[(symbol)] != alphabet->masked[(symbol)])

#define Alphabet_get_unmasked(alphabet, symbol) \
            ((alphabet)->unmasked[(symbol)_])
#define Alphabet_get_masked(alphabet, symbol) \
            ((alphabet)->masked[(symbol)])

#define Alphabet_symbol_is_valid(alphabet, symbol) \
            ((alphabet)->is_valid[(symbol)] != '-')
#define Alphabet_symbol_is_member(alphabet, symbol) \
            ((alphabet)->non_ambig[(symbol)] != '-')

#define Alphabet_symbol_complement(alphabet, symbol) \
            ((alphabet)->complement[(symbol)])
#define Alphabet_symbol_clean(alphabet, symbol) \
            ((alphabet)->clean[(symbol)])

#define Alphabet_is_DNA(alphabet) \
    ((alphabet)->type == Alphabet_Type_DNA)

/**/

gchar *Alphabet_Type_get_name(Alphabet_Type type);
Alphabet_Type Alphabet_name_get_type(gchar *name);

gchar *Alphabet_aa2tla(gchar aa);
/* Give three-letter-abbreviation for amino acid.
 * eg 'M' -> "Met" etc */

gchar *Alphabet_nt2ambig(gchar nt);
/* Returns the ambiguity bases for the given nucleotide
 * The returned string does not need to be freed.
 */


/**/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* INCLUDED_ALPHABET_H */


