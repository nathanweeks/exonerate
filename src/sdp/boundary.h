/****************************************************************\
*                                                                *
*  C4 dynamic programming library - Boundary Object              *
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

#ifndef INCLUDED_BOUNDARY_H
#define INCLUDED_BOUNDARY_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <glib.h>
#include "region.h"

/**/

typedef struct {
    gint  query_pos;
    gint  length;
    gint  seed_id;
} Boundary_Interval;

typedef struct {
         gint  target_pos;
    GPtrArray *interval_list;
} Boundary_Row;

typedef struct {
         gint  ref_count;
    GPtrArray *row_list;
} Boundary;

/**/

Boundary *Boundary_create(void);
    void  Boundary_destroy(Boundary *boundary);
Boundary *Boundary_share(Boundary *boundary);

/**/

Boundary_Row *Boundary_add_row(Boundary *boundary, gint target_pos);
Boundary_Row *Boundary_get_last_row(Boundary *boundary);
        void  Boundary_remove_empty_last_row(Boundary *boundary);

/**/

    void  Boundary_reverse(Boundary *boundary);
Boundary *Boundary_select(Boundary *boundary, Region *region);

/**/

void Boundary_Row_prepend(Boundary_Row *boundary_row,
                          gint query_pos, gint seed_id);

void Boundary_insert(Boundary *boundary, Boundary *insert);

/**/

void Boundary_print_gnuplot(Boundary *boundary, gint colour);

/**/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* INCLUDED_BOUNDARY_H */

