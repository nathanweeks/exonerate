/****************************************************************\
*                                                                *
*  A simple bitarray data structure                              *
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

#include "bitarray.h"

int main(void){
    register BitArray *ba = BitArray_create();
    register gint i;
    for(i = 0; i < 16; i++)
        BitArray_append(ba, i, 4);
    BitArray_info(ba);
    BitArray_destroy(ba);
    return 0;
    }

