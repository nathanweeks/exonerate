/****************************************************************\
*                                                                *
*  Match : A module for pairwise symbol comparison               *
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

#ifndef INCLUDED_MATCH_H
#define INCLUDED_MATCH_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <glib.h>
#include "argument.h"

#include "alphabet.h"
#include "sequence.h"
#include "translate.h"
#include "submat.h"
#include "codonsubmat.h"

typedef struct {
     gboolean  softmasked_query;
     gboolean  softmasked_target;
       Submat *dna_submat;
       Submat *protein_submat;
  /* CodonSubmat *codon_submat; */
    Translate *translate;
} Match_ArgumentSet;

Match_ArgumentSet *Match_ArgumentSet_create(Argument *arg);

/**/

typedef enum {
    Match_Type_DNA2DNA,
    Match_Type_PROTEIN2PROTEIN,
    Match_Type_DNA2PROTEIN,
    Match_Type_PROTEIN2DNA,
    Match_Type_CODON2CODON,
    Match_Type_TOTAL
} Match_Type;

Match_Type Match_Type_find(Alphabet_Type query_type,
                           Alphabet_Type target_type,
                           gboolean translate_both);
gchar *Match_Type_get_name(Match_Type type);

gchar *Match_Type_get_score_macro(Match_Type type);

/**/

typedef gint Match_Score;
#define MATCH_IMPOSSIBLY_LOW_SCORE -987654321


typedef struct Match_Strand {
             Alphabet  *alphabet;
                guint   advance;
             gboolean   is_translated;
          Match_Score (*self_func)(struct Match_Strand *strand,
                                   Sequence *seq, guint pos);
             gboolean (*mask_func)(struct Match_Strand *strand,
                                   Sequence *seq, guint pos);
         struct Match  *match;
} Match_Strand;
/* get_raw returns the raw sequence from that position
 *
 * self_func returns score of comparison to self at given position.
 *
 * mask_func returns TRUE if position contains any (soft)masked symbols.
 *
 */

void Match_Strand_get_raw(Match_Strand *strand, Sequence *sequence,
                          guint pos, gchar *result);

typedef struct Match {
    Match_ArgumentSet  *mas;
           Match_Type   type;
         Match_Strand  *query;
         Match_Strand  *target;
             Alphabet  *comparison_alphabet;
          Match_Score (*score_func)(struct Match *match,
                                    Sequence *query, Sequence *target,
                                    guint query_pos, guint target_pos);
                void  (*display_func)(struct Match *match,
                                    Sequence *query, Sequence *target,
                                    guint query_pos, guint target_pos,
                                    gchar *display_str);
          Match_Score (*split_score_func)(struct Match *match,
                                         Sequence *query, Sequence *target,
                                         guint qp1, guint qp2, guint qp3,
                                         guint tp1, guint tp2, guint tp3);
                void  (*split_display_func)(struct Match *match,
                                            Sequence *query, Sequence *target,
                                            guint qp1, guint qp2, guint qp3,
                                            guint tp1, guint tp2, guint tp3,
                                            gchar *display_str);
         struct Match  *mirror;
} Match;
/* score_func returns the score for comparison at a given position
 * display_func returns the equivalence symbol for a given position
 *
 * Match is a singleton object : created by first call to Match_find()
 * and destroyed by call to Match_destroy_all()
 */

Match *Match_find(Match_Type type);
 void  Match_destroy_all(void);

Match *Match_swap(Match *match);

Match_Score Match_max_score(Match *match);


/**/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* INCLUDED_MATCH_H */

