/****************************************************************\
*                                                                *
*  Codon Substitution Matrix Object                              *
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

#include "argument.h"
#include "codonsubmat.h"

gint Argument_main(Argument *arg){
    register CodonSubmat *cs = CodonSubmat_create();
    register guchar *alphabet = (guchar*)"ACGTN";
    register gint a, b, c, d, e, f;
    register gint max;
    guchar codon_a[4], codon_b[4];
    codon_a[3] = codon_b[3] = '\0';
    for(a = 0; a < 5; a++){
        codon_a[0] = alphabet[a];
        for(b = 0; b < 5; b++){
            codon_a[1] = alphabet[b];
            for(c = 0; c < 5; c++){
                codon_a[2] = alphabet[c];
                for(d = 0; d < 5; d++){
                    codon_b[0] = alphabet[d];
                    for(e = 0; e < 5; e++){
                        codon_b[1] = alphabet[e];
                        for(f = 0; f < 5; f++){
                            codon_b[2] = alphabet[f];
                            g_message("Score [%s]:[%s] [%d]",
                                    codon_a, codon_b,
                                    CodonSubmat_lookup(cs, codon_a,
                                                           codon_b));
                            }
                        }
                    }
                }
            }
        }
    max = CodonSubmat_max_score(cs);
    CodonSubmat_destroy(cs);
    g_message("Max value = [%d]", max);
    g_assert(max == 19);
    return 0;
    }

