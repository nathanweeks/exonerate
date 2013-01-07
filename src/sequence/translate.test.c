/****************************************************************\
*                                                                *
*  Nucleotide Translation Code                                   *
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

#include "translate.h"

static void reverse_func(gchar *dna, gint length, gpointer user_data){
    register gint *rt_count = user_data;
    g_message("Reverse translated [%s]", dna);
    (*rt_count)++;
    return;
    }

static void check_translate(Translate *t, gchar *codon, gchar crib){
    register gchar aa = Translate_codon(t, codon);
    g_message("Translate [%s] -> [%c] (crib:[%c]", codon, aa, crib);
    g_assert(aa == crib);
    return;
    }

static void check_reverse_translate(Translate *t, gchar *word, gint crib){
    int rt_count = 0;
    g_message("reverse translate [%s]", word);
    Translate_reverse(t, word, 5, reverse_func, &rt_count);
    g_message("produced [%d] reverse translations\n", rt_count);
    g_assert(rt_count == crib);
    return;
    }

gint Argument_main(Argument *arg){
    register Translate *t;
    Translate_ArgumentSet_create(arg);
    Argument_process(arg, "translate.test", NULL, NULL);
    t = Translate_create(FALSE);
    check_translate(t, "ATG", 'M');
    check_translate(t, "GGX", 'G');
    check_translate(t, "NGT", 'X');
    check_translate(t, "TGA", '*');
    /**/
    check_reverse_translate(t, "S*MXG",   72);
    check_reverse_translate(t, "SSSSS", 7776);
    check_reverse_translate(t, "MMMMM",    1);
    Translate_destroy(t);
    return 0;
    }

