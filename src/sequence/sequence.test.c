/****************************************************************\
*                                                                *
*  Simple Sequence Object                                        *
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

#include "sequence.h"

gint Argument_main(Argument *arg){
    register Alphabet *alphabet = Alphabet_create(Alphabet_Type_DNA,
                                                  TRUE);
    register gchar *seq = "CGTATGACGtagctagctGCAGATAGC";
    register Sequence *s = Sequence_create("testseq", NULL, seq, 0,
                                           Sequence_Strand_FORWARD,
                                           alphabet);
    register Sequence *s2, *s3, *s4;
    register gchar *result;
    register Translate *translate = Translate_create(FALSE);
    s2 = Sequence_revcomp(s);
    s3 = Sequence_translate(s2, translate, 1);
    s4 = Sequence_mask(s3);
    /**/
    result = Sequence_get_str(s4);
    g_message("result [%s]", result);
    g_free(result);
    /**/
    Sequence_destroy(s);
    Sequence_destroy(s2);
    Sequence_destroy(s3);
    Sequence_destroy(s4);
    Alphabet_destroy(alphabet);
    Translate_destroy(translate);
    return 0;
    }

