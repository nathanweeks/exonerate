/****************************************************************\
*                                                                *
*  Protein <-> Genome comparison model                           *
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

#include "protein2genome.h"
#include "alignment.h"
#include "optimal.h"

static void test_alignment(C4_Model *model,
                           Sequence *query, Sequence *target,
                           Protein2Genome_Data *p2gd){
    register C4_Score score;
    register Alignment *alignment;
    register Optimal *optimal = Optimal_create(model,
                                               NULL,
                                               Optimal_Type_SCORE
                                              |Optimal_Type_PATH,
                                              FALSE);
    register Region *region = Region_create(0, 0,
                                           query->len, target->len);
    score = Optimal_find_score(optimal, region, p2gd, NULL);
    g_message("Score is [%d] (expect 125)", score);
    g_assert(score == 125);
    alignment = Optimal_find_path(optimal, region, p2gd,
                                  C4_IMPOSSIBLY_LOW_SCORE, NULL);
    g_message("Alignment score is [%d]", alignment->score);
    Alignment_display(alignment, query, target,
                      NULL,
                      Protein2Genome_Data_get_submat(p2gd),
                      Protein2Genome_Data_get_translate(p2gd), stdout);
    g_assert(score == alignment->score);
    Alignment_destroy(alignment);
    Optimal_destroy(optimal);
    Region_destroy(region);
    return;
    }

static void test_protein2genome(Sequence *query, Sequence *target){
    register C4_Model *protein2genome =
        Protein2Genome_create(Affine_Model_Type_LOCAL);
    register Protein2Genome_Data *p2gd
           = Protein2Genome_Data_create(query, target);
    test_alignment(protein2genome, query, target,
                   p2gd);
    Protein2Genome_Data_destroy(p2gd);
    C4_Model_destroy(protein2genome);
    return;
    }

int Argument_main(Argument *arg){
    register Alphabet *dna_alphabet
           = Alphabet_create(Alphabet_Type_DNA, FALSE),
             *protein_alphabet
           = Alphabet_create(Alphabet_Type_PROTEIN, FALSE);
    register Sequence
        *protein = Sequence_create("protein", NULL,
                     "MADQLTEQIAEFKEAFSLFDKDGDGTITT", 0,
                     Sequence_Strand_UNKNOWN, protein_alphabet),
        *genome = Sequence_create("genome", NULL,
                     "ATGGCTGACCAGCTGACTGAGCAGATT"
                     /* "GCAGAGTTCAAG" */
                     "GCAGAGTTCAA"
                     "GTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNAG"
                     /* "GAGGCCTTCTCCCTCTTT" */
                     "GGAGGCCTTCTCCCTCTTT"
                     "GACAAGGATGGAGATGGCACTATTACCACC", 0,
                     Sequence_Strand_UNKNOWN, dna_alphabet);
    Match_ArgumentSet_create(arg);
    Affine_ArgumentSet_create(arg);
    Frameshift_ArgumentSet_create(arg);
    Intron_ArgumentSet_create(arg);
    Viterbi_ArgumentSet_create(arg);
    Argument_process(arg, "protein2genome.test", NULL, NULL);
    g_message("Testing protein2genome");
    test_protein2genome(protein, genome);
    Sequence_destroy(protein);
    Sequence_destroy(genome);
    Alphabet_destroy(dna_alphabet);
    Alphabet_destroy(protein_alphabet);
    return 0;
    }

