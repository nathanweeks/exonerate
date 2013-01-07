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

#ifndef INCLUDED_SUBMAT_H
#define INCLUDED_SUBMAT_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <limits.h>

#include <glib.h>

#ifndef ALPHABETSIZE
#define ALPHABETSIZE (1<<CHAR_BIT)
#endif /* ALPHABETSIZE */

#define SUBMAT_ALPHABETSIZE 24

typedef gint SubmatMatrix[SUBMAT_ALPHABETSIZE][SUBMAT_ALPHABETSIZE];

typedef struct {
            gint  ref_count;
    SubmatMatrix  matrix;
          guchar *index;
} Submat;

Submat *Submat_create(gchar *path);

/* FIXME: should separate Submat_create(type) and Submat_load(path)
 */

    void  Submat_destroy(Submat *s);
  Submat *Submat_share(Submat *s);
    gint  Submat_max_score(Submat *s);

/* If the path is {blosum62,pam250,nucleic,edit,identity,iupac-identity}
 * then a built-in matrix is used.
 */

#define Submat_lookup(submat, symbol_a, symbol_b)              \
       ((submat)->matrix[(submat)->index[(guchar)(symbol_a)]]  \
                        [(submat)->index[(guchar)(symbol_b)]])

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* INCLUDED_SUBMAT_H */

