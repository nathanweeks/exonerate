/****************************************************************\
*                                                                *
*  Library for basic line-by-line parsing of text files.         *
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
#include "lineparse.h"

int main(void){
    register gchar *path =
        g_strconcat(SOURCE_ROOT_DIR, G_DIR_SEPARATOR_S,
                    "src", G_DIR_SEPARATOR_S,
                    "general", G_DIR_SEPARATOR_S,
                    __FILE__, NULL);
    register FILE *fp = fopen(path, "r");
    register LineParse *lp;
    register int i;
    if(!fp)
        g_error("Could not open file [%s]", path);
    lp = LineParse_create(fp);
    while(LineParse_word(lp) != EOF){
        g_print("Line");
        for(i = 0; i < lp->word->len; i++)
            g_print("-[%s]", (gchar*)lp->word->pdata[i]);
        g_print("\n");
        }
    LineParse_destroy(lp);
    fclose(fp);
    g_free(path);
    return 0;
    }

