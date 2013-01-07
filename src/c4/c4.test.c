/****************************************************************\
*                                                                *
*  C4 dynamic programming library - code for models              *
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

#include "c4.h"

int main(void){
    register C4_Model *model = C4_Model_create("test model");
    C4_Model_destroy(model);
    return 0;
    }

/* Testing for the main c4 module is done mostly by the model tests
 */

