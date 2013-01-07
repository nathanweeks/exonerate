/****************************************************************\
*                                                                *
*  Module for EST <-> genome alignments                          *
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

#include "est2genome.h"

#include "alignment.h"
#include "optimal.h"

int Argument_main(Argument *arg){
    register C4_Model *est2genome;
    register C4_Score score;
    register Alignment *alignment;
    register Optimal *optimal;
    register gchar *query_seq =
                  "CGATCGATCGNATCGATCGATC"
                  "CATCTATCTAGCGAGCGATCTA",
                   *target_seq =
                  "CGATCGATCGATCGATCGATC"
              "GTNNNNNNNNNNNNNNNNNNNN"
              "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
              "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
              "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
              "NNNNNNNNNNNNNNNNNNNNNNNNNNNAG"
           /* "CTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNAC" */
                   "CATCTATCTANNNGCGAGCGATCTA";
    register Alphabet *dna_alphabet = Alphabet_create(
                                       Alphabet_Type_DNA, FALSE);
    register Sequence *query = Sequence_create("query", NULL,
                               query_seq, 0, Sequence_Strand_UNKNOWN,
                               dna_alphabet),
                      *target = Sequence_create("target", NULL,
                              target_seq, 0, Sequence_Strand_UNKNOWN,
                              dna_alphabet);
/**/
    register Region *region = Region_create(0, 0,
                               query->len, target->len);
    register EST2Genome_Data *e2gd;
    Match_ArgumentSet_create(arg);
    Affine_ArgumentSet_create(arg);
    Intron_ArgumentSet_create(arg);
    Viterbi_ArgumentSet_create(arg);
    Argument_process(arg, "est2genome.test", NULL, NULL);
    est2genome = EST2Genome_create();
    g_message("Has [%d] shadows", est2genome->shadow_list->len);
    optimal = Optimal_create(est2genome, NULL,
                             Optimal_Type_SCORE|Optimal_Type_PATH,
                             FALSE);
    e2gd = EST2Genome_Data_create(query, target);
    score = Optimal_find_score(optimal, region, e2gd, NULL);
    g_message("Score is [%d]", score);
    g_assert(score == 157);
/**/
    alignment = Optimal_find_path(optimal, region, e2gd,
                                  C4_IMPOSSIBLY_LOW_SCORE, NULL);
    g_message("Alignment score is [%d]", alignment->score);
    Alignment_display(alignment, query, target,
                      EST2Genome_Data_get_submat(e2gd),
                      EST2Genome_Data_get_submat(e2gd),
                      NULL, stdout);
    g_assert(score == alignment->score);
/**/
    Alphabet_destroy(dna_alphabet);
    C4_Model_destroy(est2genome);
    Alignment_destroy(alignment);
    Optimal_destroy(optimal);
    EST2Genome_Data_destroy(e2gd);
    Sequence_destroy(query);
    Sequence_destroy(target);
    Region_destroy(region);
    return 0;
    }

