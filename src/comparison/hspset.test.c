/****************************************************************\
*                                                                *
*  Library for HSP sets (high-scoring segment pairs)             *
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

#include "hspset.h"
#include "submat.h"
#include "translate.h"

typedef struct {
    guint query_start;
    guint target_start;
} TestHSPseed;

static void test_hsp_set(Match_Type match_type,
                         gchar *qy_seq, gchar *tg_seq,
                         TestHSPseed *seed, gint seed_total){
    register Sequence *query = Sequence_create("qy", NULL, qy_seq, 0,
                                Sequence_Strand_UNKNOWN, NULL),
                      *target = Sequence_create("tg", NULL, tg_seq, 0,
                                Sequence_Strand_UNKNOWN, NULL);
    register Match *match = Match_find(match_type);
    register HSP_Param *hsp_param = HSP_Param_create(match, TRUE);
    register HSPset *hsp_set = HSPset_create(query, target, hsp_param);
    register gint i;
    for(i = 0; i < seed_total; i++)
        HSPset_seed_hsp(hsp_set, seed[i].query_start,
                                 seed[i].target_start);
    HSPset_finalise(hsp_set);
    HSPset_print(hsp_set);
    HSPset_destroy(hsp_set);
    HSP_Param_destroy(hsp_param);
    Sequence_destroy(query);
    Sequence_destroy(target);
    return;
    }

gint Argument_main(Argument *arg){
    register gchar
    *ntnt_qy = "AAAAGTGAGAGAGAGAGAGAGGCGAAAAAAAAAACCCCCCCCCCACCCCGCGA",
/*                  ||||||| | ||||||||||          |||||||||| ||||||| */
    *ntnt_tg = "TTTTGTGAGAGTGTGAGAGAGGCGTTTTTTTTTTCCCCCCCCCCTCCCCGCCT";
/*              0         1         2         3         4         5  */
/*              01234567890123456789012345678901234567890123456789012*/
    register gchar *aant_qy = "PNKDEGSCPIECDFLCRHQYISDP",
     *aant_tg =
     "ACGTACGTACGTACGAGTGCGTGCCCCCTTNNNTGTGACTACATCTGCAAAACGTACGTACGT";
    HSPset_ArgumentSet_create(arg);
    Match_ArgumentSet_create(arg);
    Argument_process(arg, "hspset.test", NULL, NULL);
    TestHSPseed ntnt_seed[2] = { {  8,  8},
                                 { 36, 36} },
                ntaa_seed[1] = { { 8, 24} };
    g_message("d2d:");
    test_hsp_set(Match_Type_DNA2DNA,
                 ntnt_qy, ntnt_tg, (TestHSPseed*)ntnt_seed, 2);
    g_message("p2d:");
    test_hsp_set(Match_Type_PROTEIN2DNA,
                 aant_qy, aant_tg, (TestHSPseed*)ntaa_seed, 1);
    return 0;
    }

