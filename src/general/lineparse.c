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

#include "lineparse.h"
#include <ctype.h> /* For isspace() */

LineParse  *LineParse_create(FILE *fp){
    register LineParse *lp = g_new(LineParse, 1);
    lp->fp = fp;
    lp->line = g_string_sized_new(128);
    lp->word = g_ptr_array_new();
    return lp;
    }

void LineParse_destroy(LineParse *lp){
    g_ptr_array_free(lp->word, TRUE);
    g_string_free(lp->line, TRUE);
    g_free(lp);
    return;
    }

gint LineParse_line(LineParse *lp){
    register int ch;
    g_string_truncate(lp->line, 0);
    while((ch = getc(lp->fp)) != EOF){ /* read a line */
        if(ch == '\n')
            return lp->line->len;
        g_string_append_c(lp->line, ch);
        }
    return lp->line->len?lp->line->len:EOF;
    }

gint LineParse_word(LineParse *lp){
    register guchar *prev, *ptr;
    g_ptr_array_set_size(lp->word, 0);
    switch(LineParse_line(lp)){
        case   0: return 0;
        case EOF: return EOF;
         default: break;
        }
    /* skip start */
    for(ptr = (guchar*)lp->line->str; isspace(*ptr); ptr++);
    prev = ptr;
    while(*ptr){
        if(isspace(*ptr)){
            *ptr = '\0';
            do ptr++; while(isspace(*ptr));
            if(!*ptr)
                break;
            g_ptr_array_add(lp->word, prev); /* add a word */
            prev = ptr;
            }
        ptr++;
        }
    if(prev != ptr)
        g_ptr_array_add(lp->word, prev); /* add final word */
    return lp->word->len;
    }

gint LineParse_seek(LineParse *lp, gchar *pattern){
    register gchar *ptr = pattern;
    register gint ch, total = 0;
    while((ch = getc(lp->fp)) != EOF){
        if(*ptr == ch)
            ptr++;
        else /* Allow mismatch on pattern[0] */
            ptr = pattern + ((*pattern == ch)?1:0);
        total++;
        if(!*ptr)
            return total;
        }
    return EOF;
    }
/* Moves lp->fp to end of next occurrence of pattern.
   Returns number of characters skipped
           or EOF if pattern is not found.
*/

