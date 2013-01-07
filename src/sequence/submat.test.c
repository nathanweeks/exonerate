/****************************************************************\
*                                                                *
*  Substitution Matrix Object                                    *
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

#include "submat.h"

int main(int argc, char **argv){
    register Submat *s = Submat_create(NULL);
    g_message("test AA [%d], AT [%d], AN [%d] NN [%d]",
            Submat_lookup(s, 'A', 'A'),
            Submat_lookup(s, 'A', 'T'),
            Submat_lookup(s, 'A', 'N'),
            Submat_lookup(s, 'N', 'N'));
    Submat_destroy(s);
    return 0;
    }

