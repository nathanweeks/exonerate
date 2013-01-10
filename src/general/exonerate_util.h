/****************************************************************\
*                                                                *
*  General routines for exonerate.                               *
*                                                                *
*  Nathan Weeks         mailto:weeks@iastate.edu                 *
*  Copyright (C) 2013. All Rights Reserved.                      *
*                                                                *
*  This source code is distributed under the terms of the        *
*  GNU General Public License, version 3. See the file COPYING   *
*  or http://www.gnu.org/licenses/gpl.txt for details            *
*                                                                *
*  If you use this code, please keep this notice intact.         *
*                                                                *
\****************************************************************/

#include <ctype.h>

#ifndef EXONERATE_UTIL_H
#define EXONERATE_UTIL_H

static inline void strup (char *string) {
   char *c = string;
   while (*c)
       *(c++) = toupper((int)*c);
}

static inline void strdown (char *string) {
   char *c = string;
   while (*c)
       *(c++) = tolower((int)*c);
}

#endif
