/****************************************************************\
*                                                                *
*  Library for PCR simulation                                    *
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

#include "pcr.h"

static gboolean pcr_test_report_func(Sequence *sequence,
              PCR_Match *match_a, PCR_Match *match_b,
              gint product_length, gpointer user_data){
    g_message("Have match");
    return FALSE;
    }

gint Argument_main(Argument *arg){
    register PCR *pcr = PCR_create(pcr_test_report_func, NULL, 1, 0);
    register Alphabet *alphabet = Alphabet_create(Alphabet_Type_DNA,
                                                  FALSE);
    register Sequence *sequence = Sequence_create("testseq",
            NULL,
            "NNNNACGTNNNAACCNNNNNNNNNNNNNNNNNNNNNNNNNCCGGNNGGTTNN", 0,
            Sequence_Strand_FORWARD, alphabet);
    /* PCR_add_experiment(pcr, "test1", "ACGT", "AAKC", 10, 15); */
    PCR_add_experiment(pcr, "test1", "ACGT", "AACC", 10, 15);
    PCR_add_experiment(pcr, "test2", "CCGG", "GNGT", 20, 25);
    PCR_add_experiment(pcr, "test3", "CCG", "NGT", 20, 25);
    PCR_add_experiment(pcr, "test3", "CC", "GT", 20, 25);
    PCR_prepare(pcr);
    PCR_simulate(pcr, sequence);
    PCR_destroy(pcr);
    Sequence_destroy(sequence);
    Alphabet_destroy(alphabet);
    return 0;
    }

