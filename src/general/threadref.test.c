/****************************************************************\
*                                                                *
*  Basic library for thread-safe reference counting.             *
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

#include <stdio.h>
#include "threadref.h"

int main(void){
    register ThreadRef *tr = ThreadRef_create(), *tr2;
    tr2 = ThreadRef_share(tr);
    ThreadRef_destroy(tr2);
    ThreadRef_lock(tr);
    ThreadRef_unlock(tr);
    ThreadRef_destroy(tr);
    return 0;
    }

