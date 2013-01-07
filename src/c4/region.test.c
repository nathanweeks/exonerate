/****************************************************************\
*                                                                *
*  C4 dynamic programming library - code for regions             *
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
    register Region *region = Region_create(0, 0, 10, 10);
    Region_destroy(region);
    return 0;
    }

