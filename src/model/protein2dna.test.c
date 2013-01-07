/****************************************************************\
*                                                                *
*  Protein <-> DNA comparison model                              *
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

#include "protein2dna.h"
#include "alignment.h"
#include "optimal.h"
#include "frameshift.h"

static void test_alignment(C4_Model *model,
                           Sequence *query, Sequence *target,
                           Protein2DNA_Data *p2dd){
    register C4_Score score;
    register Alignment *alignment;
    register Optimal *optimal = Optimal_create(
                                model, NULL,
                                Optimal_Type_SCORE|Optimal_Type_PATH,
                                FALSE);
    Region region;
    Region_init_static(&region, 0, 0, query->len, target->len);
    score = Optimal_find_score(optimal, &region, p2dd, NULL);
    g_message("Score is [%d]", score);
    g_assert(score == 134);
    alignment = Optimal_find_path(optimal, &region, p2dd,
                                  C4_IMPOSSIBLY_LOW_SCORE, NULL);
    g_message("Alignment score is [%d]", alignment->score);
    Alignment_display(alignment, query, target,
                      NULL,
                      Protein2DNA_Data_get_submat(p2dd),
                      Protein2DNA_Data_get_translate(p2dd), stdout);
    g_assert(score == alignment->score);
    Alignment_destroy(alignment);
    Optimal_destroy(optimal);
    return;
    }

static void test_protein2dna(Sequence *query, Sequence *target){
    register C4_Model *protein2dna = Protein2DNA_create(Affine_Model_Type_LOCAL);
    register Protein2DNA_Data *p2dd
           = Protein2DNA_Data_create(query, target);
    test_alignment(protein2dna, query, target, p2dd);
    Protein2DNA_Data_destroy(p2dd);
    C4_Model_destroy(protein2dna);
    return;
    }

int Argument_main(Argument *arg){
    register Alphabet *dna_alphabet
           = Alphabet_create(Alphabet_Type_DNA, FALSE),
             *protein_alphabet
           = Alphabet_create(Alphabet_Type_PROTEIN, FALSE);
    register Sequence
        *dna = Sequence_create("dna", NULL,
                     "ATGGCTGACCAGCTGACTGAGGAGCAGATT"
                     "GCAGAGTTCNAAGGAGGCCTTCTCCCTCTTT"
                     "GACAAGGATGGA"
                     "NNACTGTCCATAATTGC" "TGGTACTTCAGCGGTCGATGG"
                     "GATGGCACTCTGACCACC", 0,
                     Sequence_Strand_UNKNOWN, dna_alphabet),
        *protein = Sequence_create("protein", NULL,
                     "NNNNNNMADQLTEQIAEFKEAFSLFDKDG"
                     "TVHNC" "X" "WYFSGRW"
                     "DGTITT", 0,
                     Sequence_Strand_UNKNOWN, protein_alphabet);
    Match_ArgumentSet_create(arg);
    Affine_ArgumentSet_create(arg);
    Frameshift_ArgumentSet_create(arg);
    Argument_process(arg, "dna2protein.test", NULL, NULL);
    g_message("Testing protein2dna");
    test_protein2dna(protein, dna);
    Sequence_destroy(dna);
    Sequence_destroy(protein);
    Alphabet_destroy(dna_alphabet);
    Alphabet_destroy(protein_alphabet);
    return 0;
    }

