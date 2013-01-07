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

#include "coding2coding.h"
#include "alignment.h"
#include "optimal.h"
#include "submat.h"

static void test_coding2coding(Sequence *query, Sequence *target){
    register C4_Score score;
    register C4_Model *model = Coding2Coding_create();
    register Coding2Coding_Data *c2cd = Coding2Coding_Data_create(
                                               query, target);
    register Alignment *alignment;
    register Optimal *optimal = Optimal_create(model, NULL,
                                               Optimal_Type_SCORE
                                              |Optimal_Type_PATH,
                                              FALSE);
    Region region;
    Region_init_static(&region, 0, 0, query->len, target->len);
    score = Optimal_find_score(optimal, &region, c2cd, NULL);
    g_message("Score is [%d]", score);
    g_assert(score == 169);
    alignment = Optimal_find_path(optimal, &region, c2cd,
                                  C4_IMPOSSIBLY_LOW_SCORE, NULL);
    g_message("Alignment score is [%d]", alignment->score);
    Alignment_display(alignment, query, target,
                     Coding2Coding_Data_get_submat(c2cd),
                     Coding2Coding_Data_get_submat(c2cd),
                     Coding2Coding_Data_get_translate(c2cd), stdout);
    g_assert(score == alignment->score);
    Coding2Coding_Data_destroy(c2cd);
    C4_Model_destroy(model);
    return;
    }

int Argument_main(Argument *arg){
    register Alphabet *alphabet = Alphabet_create(Alphabet_Type_DNA,
                                                  FALSE);
    register Sequence
        *qy = Sequence_create("qy", NULL,
"AGCCCAGCCAAGCACTGTCAGGAATCCTGTGAAGCAGCTCCAGCTATGTGTGAAGAAG"
"AGGACAGCACTGCCTTGGTGTGTGACAATGGCTCTGGGCTCTGTAAGGCCGGCTTTGCT", 0,
                     Sequence_Strand_UNKNOWN, alphabet),
        *tg = Sequence_create("tg", NULL,
   "AGCCCAGCCAAACACTGTCAGGAATCCTGT"
"NNN"
"GAAGCAGCTCCAGCTATGTGTGAAGAAG"
"AGGACAGCACTGCCTTGGTGTGTGACAATGGC"
"NN"
"TCTGGGCTCTGTAAGGCCGGCTTTGCT", 0,
                     Sequence_Strand_UNKNOWN, alphabet);
    Match_ArgumentSet_create(arg);
    Affine_ArgumentSet_create(arg);
    Frameshift_ArgumentSet_create(arg);
    Viterbi_ArgumentSet_create(arg);
    Argument_process(arg, "coding2coding.test", NULL, NULL);
    g_message("Testing coding2coding");
    test_coding2coding(qy, tg);
    Sequence_destroy(qy);
    Sequence_destroy(tg);
    Alphabet_destroy(alphabet);
    return 0;
    }

