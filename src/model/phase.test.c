/****************************************************************\
*                                                                *
*  Model for phased introns.                                     *
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
#include "phase.h"

gint Argument_main(Argument *arg){
    register Match *match;
    register C4_Model *phase;
    Match_ArgumentSet_create(arg);
    Argument_process(arg, "phase.test", NULL, NULL);
    match = Match_find(Match_Type_DNA2PROTEIN);
    phase = Phase_create("test", match, TRUE, FALSE);
    C4_Model_destroy(phase);
    Match_destroy_all();
    return 0;
    }

