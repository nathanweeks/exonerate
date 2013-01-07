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

#ifndef INCLUDED_LINEPARSE_H
#define INCLUDED_LINEPARSE_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <glib.h>
#include <stdio.h>

typedef struct {
     FILE  *fp;
  GString  *line;
GPtrArray  *word;
} LineParse;

LineParse  *LineParse_create(FILE *fp);
     void   LineParse_destroy(LineParse *lp);
     gint   LineParse_line(LineParse *lp);
     gint   LineParse_word(LineParse *lp);
     gint   LineParse_seek(LineParse *lp, gchar *pattern);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* INCLUDED_LINEPARSE_H */

