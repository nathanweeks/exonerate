/****************************************************************\
*                                                                *
*  Module for splice site and intron modelling                   *
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

#include "intron.h"

int Argument_main(Argument *arg){
    register C4_Model *intron = Intron_create("test", TRUE, FALSE,
                                              TRUE);
    g_warning("test does nothing");
    C4_Model_destroy(intron);
    return 0;
    }

