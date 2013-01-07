/****************************************************************\
*                                                                *
*  C4 dynamic programming library - DP layout code               *
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

#include <glib.h>

#include "optimal.h"
#include "match.h"

int Argument_main(Argument *arg){
    register C4_Model *model = C4_Model_create("test");
    register Optimal *optimal;
    register Match *match;
    Match_ArgumentSet_create(arg);
    Argument_process(arg, "optimal.test", NULL, NULL);
    match = Match_find(Match_Type_DNA2DNA);
    C4_Model_add_transition(model, "start to end", NULL, NULL, 1, 1,
                            NULL, C4_Label_MATCH, match);
    C4_Model_close(model);
    optimal = Optimal_create(model, "optimal test",
                             Optimal_Type_SCORE, FALSE);
    /**/
    Optimal_destroy(optimal);
    C4_Model_destroy(model);
    return 0;
    }

