/****************************************************************\
*                                                                *
*  Module for various ungapped alignment models                  *
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

#include "ungapped.h"
#include "alignment.h"
#include "argument.h"
#include "submat.h"
#include "optimal.h"
#include "alignment.h"

static void run_ungapped_test(Sequence *query, Sequence *target,
        gboolean translate_both, gint crib){
    register Match_Type match_type
           = Match_Type_find(query->alphabet->type,
                             target->alphabet->type,
                             translate_both);
    register Ungapped_Data *ud = Ungapped_Data_create(query, target,
                                                      match_type);
    register C4_Model *ungapped = Ungapped_create(match_type);
    register C4_Score score;
    register Alignment *alignment;
    register Optimal *optimal = Optimal_create(ungapped, NULL,
                                               Optimal_Type_SCORE
                                              |Optimal_Type_PATH,
                                              FALSE);
    g_message("Using Match_Type [%s]",
            Match_Type_get_name(match_type));
    Region region;
    Region_init_static(&region, 0,0, query->len, target->len);
    score = Optimal_find_score(optimal, &region, ud, NULL);
    g_message("Score is [%d] crib [%d]", score, crib);
    g_assert(score == crib);
    alignment = Optimal_find_path(optimal, &region, ud,
                                  C4_IMPOSSIBLY_LOW_SCORE, NULL);
    g_message("Alignment score is [%d] expect [%d]",
               alignment->score, score);
    g_assert(score == alignment->score);
    g_message("Using model [%s]", ungapped->name);
    Alignment_display(alignment, query, target,
              Ungapped_Data_get_dna_submat(ud),
              Ungapped_Data_get_protein_submat(ud),
              Ungapped_Data_get_translate(ud), stdout);
    Alignment_destroy(alignment);
    Optimal_destroy(optimal);
    C4_Model_destroy(ungapped);
    Ungapped_Data_destroy(ud);
    return;
    }

int Argument_main(Argument *arg){
    register Alphabet
     *dna_alphabet = Alphabet_create(Alphabet_Type_DNA, FALSE),
     *protein_alphabet = Alphabet_create(Alphabet_Type_PROTEIN, FALSE);
    register Sequence *dna1 = Sequence_create("dna1", NULL,
             "CGATCAGCTAGCTAGCTACGATCGATCGAT", 0,
             Sequence_Strand_UNKNOWN, dna_alphabet);
    register Sequence *dna2 = Sequence_create("dna2", NULL,
             "CGATACGATCGCTCTGAGATCTCGACTCAG", 0,
             Sequence_Strand_UNKNOWN, dna_alphabet);
    register Sequence *protein1 = Sequence_create("protein1", NULL,
             "DFCPIECFLNHILCIPEF", 0,
             Sequence_Strand_UNKNOWN, protein_alphabet);
    register Sequence *protein2 = Sequence_create("protein2", NULL,
             "DFHLISCPIRYLICIEFP", 0,
             Sequence_Strand_UNKNOWN, protein_alphabet);
    Match_ArgumentSet_create(arg);
    Argument_process(arg, "ungapped.test", NULL, NULL);
    run_ungapped_test(dna1, dna2, FALSE, 47);
    run_ungapped_test(protein1, protein2, FALSE, 36);
    run_ungapped_test(dna1, protein1, FALSE, 7);
    run_ungapped_test(protein2, dna2, FALSE, 12);
    run_ungapped_test(dna1, dna2, TRUE, 22);
    Sequence_destroy(dna1);
    Sequence_destroy(dna2);
    Sequence_destroy(protein1);
    Sequence_destroy(protein2);
    Alphabet_destroy(dna_alphabet);
    Alphabet_destroy(protein_alphabet);
    return 0;
    }

