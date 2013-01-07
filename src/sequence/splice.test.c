/****************************************************************\
*                                                                *
*  Library for Splice site prediction.                           *
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

#include "splice.h"
#include <string.h> /* For strlen() */

static void check_splice_prediction(SplicePredictor *sp, gchar *seq,
                                    gint crib){
    gfloat score;
    register gint len = strlen(seq);
    register gint pos = SplicePredictor_predict(sp, seq, len, 0, len, &score);
    g_message("pos [%d] crib [%d] score [%f]", pos, crib, score);
    g_assert(pos == crib);
    g_assert(score > 0);
    return;
    }

gint Argument_main(Argument *arg){
    register SplicePredictorSet *sps = SplicePredictorSet_create();
    g_message("Checking splice site predition [%f,%f,%f,%f]",
            SplicePredictor_get_max_score(sps->ss5_forward),
            SplicePredictor_get_max_score(sps->ss5_reverse),
            SplicePredictor_get_max_score(sps->ss3_forward),
            SplicePredictor_get_max_score(sps->ss3_reverse));
    check_splice_prediction(sps->ss5_forward,
     /*                     v                       */
      "TGGCGCTCCTGGCCTCTGCCCGTAAGCACTTGGTGGGACTGGG", 21);
     /*               (21) |GT... 5ss fwd           */
    check_splice_prediction(sps->ss3_forward,
     /*                  v                          */
      "CAGCCCTGCTCTTTCCTCAGGAGCTTCAGAGGCCGAGGATGCC", 18);
     /*       3ss fwd ...AG| (18)                   */
    check_splice_prediction(sps->ss5_reverse,
     /*                    v                        */
      "CCCAGTCCCACCAAGTGCTTACGGGCAGAGGCCAGGAGCGCCA", 20);
     /*        3ss rev ... AC| (20)                 */
    check_splice_prediction(sps->ss3_reverse,
     /*                       v                     */
      "GGCATCCTCGGCCTCTGAAGCTCCTGAGGAAAGAGCAGGGCTG", 23);
     /*                      |CT... 5ss rev (23)    */
    SplicePredictorSet_destroy(sps);
    return 0;
    }

