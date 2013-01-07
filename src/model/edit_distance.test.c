/****************************************************************\
*                                                                *
*  edit distance model                                           *
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

#include "edit_distance.h"
#include "alignment.h"
#include "optimal.h"
#include "sequence.h"

int Argument_main(Argument *arg){
    register C4_Model *edit_distance = EditDistance_create();
    register Alphabet *alphabet = Alphabet_create(Alphabet_Type_DNA,
                                                  FALSE);
    register Sequence
        *query = Sequence_create("qy", NULL,
                     "gtgcactacgtacgtnatcgtgcttnaacgcg"
                     "tacgtgatngtgcttgaacgtacgtacgtgatcg"
                     "tgcttga", 0,
                      Sequence_Strand_UNKNOWN, alphabet),
        *target = Sequence_create("tg", NULL,
                     "actacgtacgtgatcgtgcaacgcactacg"
                     "tacgtgancttgaacgcactacgtacgtgatcg"
                     "tgcntgaacgn", 0,
                      Sequence_Strand_UNKNOWN, alphabet);
    register gchar *qy = Sequence_get_str(query),
                   *tg = Sequence_get_str(target);
    register EditDistance_Data *edd = EditDistance_Data_create(qy, tg);
/**/
    register C4_Score score;
    register Alignment *alignment;
    register Optimal *optimal = Optimal_create(edit_distance,
                                               NULL,
                                               Optimal_Type_SCORE
                                              |Optimal_Type_PATH,
                                               FALSE);
/**/
    register Region *region = Region_create(0, 0,
                                     edd->query_len, edd->target_len);
    score = Optimal_find_score(optimal, region, edd, NULL);
    g_message("Score is [%d] (expect -23)", score);
    g_assert(score == -23);
/**/
    alignment = Optimal_find_path(optimal, region, edd,
                                  C4_IMPOSSIBLY_LOW_SCORE, NULL);
    g_message("Alignment score is [%d] (expect -23)",
              alignment->score);
    Alignment_display(alignment, query, target, NULL, NULL, NULL, stdout);
    g_assert(score == alignment->score);
/**/
    C4_Model_destroy(edit_distance);
    Alignment_destroy(alignment);
    Optimal_destroy(optimal);
    EditDistance_Data_destroy(edd);
    Sequence_destroy(query);
    Sequence_destroy(target);
    Alphabet_destroy(alphabet);
    Region_destroy(region);
    g_free(qy);
    g_free(tg);
    return 0;
    }

