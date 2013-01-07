/****************************************************************\
*                                                                *
*  C4 dynamic programming library - SDP Lookahead Object         *
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

#ifndef INCLUDED_LOOKAHEAD_H
#define INCLUDED_LOOKAHEAD_H 
#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <limits.h>  /* For CHAR_BIT */
#include <strings.h> /* For ffs() */
#include <glib.h>

#ifdef ffsll
    typedef long long int Lookahead_Mask;
#define Lookahead_Mask_ffs ffsll
#else /* ffsll */
    typedef int Lookahead_Mask;
#define Lookahead_Mask_ffs ffs
#endif /* ffsll */

#define Lookahead_Mask_WIDTH (sizeof(Lookahead_Mask)*(CHAR_BIT))
    /* It is assumed that Lookahead_Mask_WIDTH is always greather
     * than the model max advance.
     */

/**/

typedef void (*Lookahead_FreeFunc)(gpointer data, gpointer user_data);

typedef struct {
                  gint  reset_pos;
                  gint  pos;
                  gint  max_advance;
              gpointer *index[Lookahead_Mask_WIDTH];
        Lookahead_Mask  mask;
    Lookahead_FreeFunc  free_func;
              gpointer  user_data;
} Lookahead;

Lookahead *Lookahead_create(gint reset_pos, gint max_advance,
                            Lookahead_FreeFunc free_func,
                            gpointer user_data);
     void  Lookahead_destroy(Lookahead *lookahead);
     void  Lookahead_reset(Lookahead *lookahead);

 gpointer  Lookahead_get(Lookahead *lookahead, gint advance);
     void  Lookahead_set(Lookahead *lookahead, gint advance,
                         gpointer data);

     gint  Lookahead_next(Lookahead *lookahead);
     void  Lookahead_move(Lookahead *lookahead, gint pos);
/* Lookahead_next() will assume that the curr is set and clear it
 *                  and then return new position
 * Lookahead_move() will set pos as curr clearing any data inbetween
 */

/* FIXME: optimisation : replace get/set with macros ? */

/**/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* INCLUDED_LOOKAHEAD_H */

