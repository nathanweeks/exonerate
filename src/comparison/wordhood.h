/****************************************************************\
*                                                                *
*  Library for word-neighbourhood generation                     *
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

#ifndef INCLUDED_WORDHOOD_H
#define INCLUDED_WORDHOOD_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <glib.h>
#include <stdio.h>
#include <limits.h>

#include "submat.h"
#include "codonsubmat.h"

#ifndef ALPHABETSIZE
#define ALPHABETSIZE (1<<(CHAR_BIT))
#endif /* ALPHABETSIZE */

typedef struct WordHood_Alphabet {
         gint   ref_count;
         gint   advance;
     gpointer   user_data;
         gint   input_index[ALPHABETSIZE];
         gint   output_index[ALPHABETSIZE];
    GPtrArray  *member_list;
       Submat  *submat;
  CodonSubmat  *codon_submat;
         gint (*score_func)(struct WordHood_Alphabet *wha,
                            gchar *seq_a, gchar *seq_b);
     gboolean (*is_valid_func)(struct WordHood_Alphabet *wha,
                               gchar *seq);
         gint (*index_func)(struct WordHood_Alphabet *wha, gchar *seq);
} WordHood_Alphabet;
/* Invalid input gets -1 from the index_func
 * Output is possible for all of member_list
 */

WordHood_Alphabet *WordHood_Alphabet_create_from_Submat(
                        gchar *input_alphabet,
                        gchar *output_alphabet, Submat *submat,
                        gboolean case_sensitive_input);
WordHood_Alphabet *WordHood_Alphabet_create_from_CodonSubmat(
                    CodonSubmat *cs, gboolean case_sensitive_input);
WordHood_Alphabet *WordHood_Alphabet_share(WordHood_Alphabet *wha);
void WordHood_Alphabet_destroy(WordHood_Alphabet *wha);

/* Using different input and output alphabets allows
 * generation of a neighbourhood covering redundant symbols
 *
 */

typedef gboolean (*WordHood_Traverse_Func)(gchar *word,
                  gint score, gpointer user_data);
/* Return TRUE to stop the traversal */

typedef struct {
         WordHood_Alphabet *wha;
                      gint  threshold;
                  gboolean  use_dropoff;
                     gchar *orig_word;
                     gchar *curr_word;
                      gint *depth_threshold;
                      gint  word_pos;
                      gint  curr_score;
                      gint  curr_len;
                      gint  alloc_len;
} WordHood;

WordHood *WordHood_create(WordHood_Alphabet *wha,
                          gint threshold, gboolean use_dropoff);

/* When use_dropoff is FALSE:
 *   - each word will have a score of at least the threshold
 * When use_dropoff is TRUE:
 *   - each word will have a score within the threshold
 *     of the score given by comparison to itself
 *
 * Thus threshold=0, use_dropoff=TRUE will yield the minimum wordhood
 * (ie. just valid words from the query).
 */

void WordHood_destroy(WordHood *wh);
void WordHood_info(WordHood *wh);

void WordHood_traverse(WordHood *wh, WordHood_Traverse_Func whtf,
                       gchar *word, gint len, gpointer user_data);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* INCLUDED_WORDHOOD_H */

