/****************************************************************\
*                                                                *
*  C4 dynamic programming library - alignment code               *
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

#include "alignment.h"
#include "match.h"

gint Argument_main(Argument *arg){
    register Region *region = Region_create_blank();
    register C4_Model *model = C4_Model_create("test");
    register Alignment *alignment;
    register Match *match;
    Match_ArgumentSet_create(arg);
    Argument_process(arg, "alignment.test", NULL, NULL);
    match = Match_find(Match_Type_DNA2DNA);
    C4_Model_add_transition(model, "start to end",
                            NULL, NULL, 1, 1, NULL, C4_Label_MATCH,
                            match);
    C4_Model_close(model);
    alignment = Alignment_create(model, region, 0);
    g_message("alignment test");
    Alignment_destroy(alignment);
    C4_Model_destroy(model);
    Region_destroy(region);
    return 0;
    }
/* Proper testing is done with the model tests
 */

