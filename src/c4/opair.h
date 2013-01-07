/****************************************************************\
*                                                                *
*  C4 dynamic programming library - code for optimal pairs       *
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

#ifndef INCLUDED_OPAIR_H
#define INCLUDED_OPAIR_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <glib.h>
#include "c4.h"
#include "alignment.h"
#include "optimal.h"
#include "subopt.h"

typedef struct {
    gpointer   user_data;
     Optimal  *optimal;
      Region  *region;
      SubOpt  *subopt;
    C4_Score   prev_score;
} OPair;

OPair *OPair_create(Optimal *optimal, SubOpt *subopt,
                    gint query_length, gint target_length,
                    gpointer user_data);
void OPair_destroy(OPair *opair);
/* user_data is provided for the Optimal functions */

Alignment *OPair_next_path(OPair *hpair, C4_Score threshold);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* INCLUDED_OPAIR_H */

