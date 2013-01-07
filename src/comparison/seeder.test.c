/****************************************************************\
*                                                                *
*  Seeder : A module for seeding pairwise alignments             *
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

#include "seeder.h"
#include "argument.h"

static void Seeder_Test_report_func(Comparison *comparison,
                                    gpointer user_data){
    register gint *count = user_data;
    register gint i;
    register HSP *hsp;
    g_message("Have comparison with [%d] DNA hsps",
              comparison->dna_hspset->hsp_list->len);
    for(i = 0; i < comparison->dna_hspset->hsp_list->len; i++){
        hsp = comparison->dna_hspset->hsp_list->pdata[i];
        HSP_print(hsp, "hsp");
        }
    (*count) += comparison->dna_hspset->hsp_list->len;
    return;
    }

gint Argument_main(Argument *arg){
    register Match *match;
    register HSP_Param *dna_hsp_param;
    register Comparison_Param *comparison_param;
    register Seeder *seeder;
    register Alphabet *dna_alphabet = Alphabet_create(Alphabet_Type_DNA,
                                                  FALSE);
    register Sequence *query = Sequence_create("qy", NULL,
        "ACGTAACCGGTTAGCT", 0,
        Sequence_Strand_FORWARD, dna_alphabet),
                  *target = Sequence_create("tg", NULL,
        "NNNNNNNNNNNNNNNNNACGTAACCGGTTAGCTNNNNNNNNNNNNNNNNNNN"
        "NNACGTAACCGGTTAGCTNNNNNNNNACGTAACCGGTTAGCTNNNNNNNNNN", 0,
        Sequence_Strand_FORWARD, dna_alphabet);
    gint count = 0;
    Match_ArgumentSet_create(arg);
    HSPset_ArgumentSet_create(arg);
    Seeder_ArgumentSet_create(arg);
    Argument_process(arg, "seeder.test", NULL, NULL);
    match = Match_find(Match_Type_DNA2DNA);
    dna_hsp_param = HSP_Param_create(match, TRUE);
    comparison_param = Comparison_Param_create(Alphabet_Type_DNA,
                                               Alphabet_Type_DNA,
                                               dna_hsp_param,
                                               NULL, NULL);
    seeder = Seeder_create(1, comparison_param, 0,
                           Seeder_Test_report_func, &count);
    g_message("Seeder test");
    Seeder_add_query(seeder, query);
    Seeder_add_target(seeder, target);
    /**/
    g_message("final count [%d]", count);
    g_assert(count == 3);
    /**/
    Alphabet_destroy(dna_alphabet);
    Sequence_destroy(query);
    Sequence_destroy(target);
    HSP_Param_destroy(dna_hsp_param);
    Comparison_Param_destroy(comparison_param);
    Seeder_destroy(seeder);
    return 0;
    }

