/****************************************************************\
*                                                                *
*  Module for various affine gapped models                       *
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

#include "affine.h"
#include "alignment.h"
#include "optimal.h"

static void test_model(Affine_Model_Type type, C4_Score expected_score){
    register C4_Model *affine = Affine_create(type,
                                              Alphabet_Type_PROTEIN,
                                              Alphabet_Type_PROTEIN,
                                              FALSE);
/**/
    register C4_Score score;
    register Alignment *alignment;
    register gchar *name = NULL;
    register Optimal *optimal;
    register Alphabet *alphabet = Alphabet_create(Alphabet_Type_PROTEIN,
                                                  FALSE);
    register Sequence
        *query = Sequence_create("qy", NULL,
                     "MEEPQSDPSVEPPLSQETFSDLWKLL", 0,
                     Sequence_Strand_UNKNOWN, alphabet),
        *target = Sequence_create("tg", NULL,
                     "PENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP"
                     "EHSCETFDIWKWCPIECDFLNVISEPNEPIPSQ", 0,
                     Sequence_Strand_UNKNOWN, alphabet);
    register Affine_Data *ad = Affine_Data_create(query, target, FALSE);
      /*
      "TVPEPIVTEPTTITEPEVPEKEEPKAEVEKTKKAKGSKPKKASKPRNPA"
      "SHPTYEEMIKDAIVSLKEKNGSSQYAIA",
      "VVEEQAAPETVKDEANPPAKSGKAKKETKAKKPAAPRKRSATPTHPPYF"
      "EMIKDAIVTLKERTGSSQHAIT",
      */
      /* Paracelsus challenge sequences: */
      /*
      "MTYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE",
      "MTKKAILALNTAKFLRTQAAVLAAKLEKLGAQEANDNAVDLEDTADDLYKTLLVLA",
      */
    Region region;
/**/
    switch(type){
        case Affine_Model_Type_GLOBAL:
            name = "global affine";
            break;
        case Affine_Model_Type_BESTFIT:
            name = "bestfit affine";
            break;
        case Affine_Model_Type_LOCAL:
            name = "local affine";
            break;
        case Affine_Model_Type_OVERLAP:
            name = "overlap affine";
            break;
        default:
            g_error("Model Type not supported [%d]", type);
            break;
        }
    optimal = Optimal_create(affine, NULL,
                             Optimal_Type_SCORE|Optimal_Type_PATH,
                             FALSE);
/**/
    Region_init_static(&region, 0, 0, ad->ud.query->len,
                                      ad->ud.target->len);
    score = Optimal_find_score(optimal, &region, ad, NULL);
    g_message("Score for [%s] [%d]", name, score);
    g_assert(score == expected_score);
/**/
    alignment = Optimal_find_path(optimal, &region, ad,
                                  C4_IMPOSSIBLY_LOW_SCORE, NULL);
    g_message("Alignment score is [%d]", alignment->score);
    Alignment_display(alignment, query, target,
            Affine_Data_get_dna_submat(ad),
            Affine_Data_get_protein_submat(ad), NULL, stdout);
    Alignment_display_vulgar(alignment, query, target, stdout);
    g_assert(score == alignment->score);
/**/
    C4_Model_destroy(affine);
    Alignment_destroy(alignment);
    Optimal_destroy(optimal);
    Affine_Data_destroy(ad);
    Sequence_destroy(query);
    Sequence_destroy(target);
    Alphabet_destroy(alphabet);
    return;
    }

int Argument_main(Argument *arg){
    Match_ArgumentSet_create(arg);
    Affine_ArgumentSet_create(arg);
    Argument_process(arg, "affine.test", NULL, NULL);
    test_model(Affine_Model_Type_GLOBAL, -151);
    test_model(Affine_Model_Type_BESTFIT, 18);
    test_model(Affine_Model_Type_LOCAL, 32);
    test_model(Affine_Model_Type_OVERLAP, 18);
    return 0;
    }



