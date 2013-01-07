/****************************************************************\
*                                                                *
*  Model for alignments with non-equivalenced regions            *
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

#include "ner.h"
#include "alignment.h"
#include "optimal.h"

int Argument_main(Argument *arg){
    register C4_Score score;
    register C4_Model *ner;
    register Alignment *alignment;
    register Optimal *optimal;
    register Alphabet *dna_alphabet
           = Alphabet_create(Alphabet_Type_DNA, FALSE);
    register Sequence
        *query = Sequence_create("qy", NULL,
       "TTTTATCTTCCCAAGAGNCCCCATNNNGCGA"
       "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
       "AAAAAAAAAAAAAA"
       "GTGATTGAAATGTGGATGAAACATTTC", 0,
       Sequence_Strand_UNKNOWN, dna_alphabet),
        *target = Sequence_create("tg", NULL,
       "TTTTATCTTCCCAAGAGCCCCATGAGGCGA"
       "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
       "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
       "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
       "GTGANTGAAATGTGGATGAACATTTC", 0,
       Sequence_Strand_UNKNOWN, dna_alphabet);
    register NER_Data *nd;
/**/
    Region region;
    Match_ArgumentSet_create(arg);
    Affine_ArgumentSet_create(arg);
    NER_ArgumentSet_create(arg);
    Argument_process(arg, "ner.test", NULL, NULL);
    ner = NER_create(Alphabet_Type_DNA, Alphabet_Type_DNA);
    nd = NER_Data_create(query, target);
    optimal = Optimal_create(ner, NULL,
                             Optimal_Type_SCORE|Optimal_Type_PATH,
                             FALSE);
/**/
    Region_init_static(&region, 0, 0, query->len, target->len);
    score = Optimal_find_score(optimal, &region, nd, NULL);
    g_message("Score is [%d]", score);
    g_assert(score == 208);
/**/
    alignment = Optimal_find_path(optimal, &region, nd,
                                  C4_IMPOSSIBLY_LOW_SCORE, NULL);
    g_message("Alignment score is [%d]", alignment->score);
    Alignment_display(alignment, query, target,
            NER_Data_get_dna_submat(nd),
            NER_Data_get_protein_submat(nd), NULL, stdout);
    g_assert(score == alignment->score);
/**/
    Alignment_destroy(alignment);
    Optimal_destroy(optimal);
    NER_Data_destroy(nd);
/**/
    C4_Model_destroy(ner);
    Alphabet_destroy(dna_alphabet);
    Sequence_destroy(query);
    Sequence_destroy(target);
    return 0;
    }

