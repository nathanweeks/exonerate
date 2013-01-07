/****************************************************************\
*                                                                *
*  Coding <-> Genome comparison model                            *
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

#include "coding2genome.h"
#include "alignment.h"
#include "optimal.h"
#include "submat.h"

static void test_coding2genome(Sequence *query, Sequence *target){
    register C4_Score score;
    register C4_Model *model = Coding2Genome_create();
    register Coding2Genome_Data *c2gd = Coding2Genome_Data_create(
                                               query, target);
    register Alignment *alignment;
    register Optimal *optimal = Optimal_create(model, NULL,
                                               Optimal_Type_SCORE
                                              |Optimal_Type_PATH,
                                              FALSE);
    register Region *region = Region_create(0, 0,
                                            query->len, target->len);
    score = Optimal_find_score(optimal, region, c2gd, NULL);
    g_message("Score is [%d]", score);
    /* g_assert(score == 151); */
    alignment = Optimal_find_path(optimal, region, c2gd,
                                  C4_IMPOSSIBLY_LOW_SCORE, NULL);
    g_message("Alignment score is [%d]", alignment->score);
    Alignment_display(alignment, query, target,
                     Coding2Genome_Data_get_submat(c2gd),
                     Coding2Genome_Data_get_submat(c2gd),
                     Coding2Genome_Data_get_translate(c2gd), stdout);
    g_assert(score == alignment->score);
    Coding2Genome_Data_destroy(c2gd);
    C4_Model_destroy(model);
    Region_destroy(region);
    Optimal_destroy(optimal);
    return;
    }

int Argument_main(Argument *arg){
    register Alphabet *alphabet = Alphabet_create(Alphabet_Type_DNA,
                                                  FALSE);
    register Sequence
        *qy = Sequence_create("qyr", NULL,
"AGCCCAGCCAAGCACTGTCAGGAATCCTGT"
"GAAGCAGCTCCAGCTATGTGTGAAGAA"
"GAGGACAGCACTGCCTTGGTGTGTGACAATG"
/* "GTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNAG" */
"GCTCTGGGCTCTGTAAGGCCGGCTTTGCT", 0,
                     Sequence_Strand_UNKNOWN, alphabet),
        *tg = Sequence_create("tgr", NULL,
"AGCCCAGCCAAGCACTGTCAGGAATCCTG"
"GTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNAG"
"T"
"GAAGCAGCTCCAGCTATGTGTGAAGAA"
"GAGGACAGCACTGCCTTGGTGTGTGACAATGGC"
"TCTGGGCTCTGTAAGGCCGGCTTTGCT", 0,
                     Sequence_Strand_UNKNOWN, alphabet);
    Match_ArgumentSet_create(arg);
    Affine_ArgumentSet_create(arg);
    Viterbi_ArgumentSet_create(arg);
    Frameshift_ArgumentSet_create(arg);
    Intron_ArgumentSet_create(arg);
    Argument_process(arg, "coding2genome.test", NULL, NULL);
    g_message("Testing coding2genome");
    g_message("seq lengths [%d,%d]", qy->len, tg->len);
    test_coding2genome(qy, tg);
    Sequence_destroy(qy);
    Sequence_destroy(tg);
    Alphabet_destroy(alphabet);
    return 0;
    }


