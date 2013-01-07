/****************************************************************\
*                                                                *
*  Simple Alphabet Object                                        *
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

#include "alphabet.h"

gint Argument_main(Argument *arg){
    register Alphabet *dna_alphabet
           = Alphabet_create(Alphabet_Type_DNA, FALSE);
    Alphabet_destroy(dna_alphabet);
    return 0;
    }

