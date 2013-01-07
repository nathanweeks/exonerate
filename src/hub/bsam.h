/****************************************************************\
*                                                                *
*  BSAM: Big Sequence Alignment Manager                          *
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

#ifndef INCLUDED_BSAM_H
#define INCLUDED_BSAM_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <glib.h>
#include "comparison.h"

/**/

typedef struct {
    Comparison_Param *comparison_param;
                gint  saturate_threshold;
                gint  verbosity;
} BSAM;

BSAM *BSAM_create(Comparison_Param *comparison_param,
                  gint saturate_threshold, gint verbosity);
void  BSAM_destroy(BSAM *bsam);

Comparison *BSAM_compare(BSAM *bsam, Sequence *query, Sequence *target);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* INCLUDED_BSAM_H */

