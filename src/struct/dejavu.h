/****************************************************************\
*                                                                *
*  Deja-vu library for fast linear space repeat finding          *
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

#ifndef INCLUDED_DEJAVU_H
#define INCLUDED_DEJAVU_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <glib.h>
#include <limits.h>

#ifndef ALPHABET_SIZE
#define ALPHABET_SIZE (1 << CHAR_BIT)
#endif /* ALPHABET_SIZE */

typedef struct {
       gchar *seq;
        gint  len;
        gint *next;
        gint  symbol_list[ALPHABET_SIZE];
        gint  symbol_list_len;
} DejaVu;

typedef void (*DejaVu_TraverseFunc)(gint first_pos, gint curr_pos,
                                    gint length, gchar *seq, gint len,
                                    gpointer user_data);

DejaVu *DejaVu_create(gchar *seq, gint len);
/* seq is used, without a copy being made */

void DejaVu_info(DejaVu *dv);
void DejaVu_destroy(DejaVu *dv);
void DejaVu_traverse(DejaVu *dv,
                     gint min_wordlen, gint max_wordlen,
                     DejaVu_TraverseFunc dvtf, gpointer user_data,
                     guchar *filter, gint verbosity);
/* No filter applied when filter is NULL.
 * Filter is Alphabet size.
 * symbols are excluded when (filter[seq[i]] == '-')
 */

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* INCLUDED_DEJAVU_H */

