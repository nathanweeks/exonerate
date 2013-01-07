/****************************************************************\
*                                                                *
*  Model for phased introns.                                     *
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

#ifndef INCLUDED_PHASE_H
#define INCLUDED_PHASE_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "c4.h"
#include "sequence.h"
#include "match.h"

C4_Model *Phase_create(gchar *suffix, Match *match,
                       gboolean on_query, gboolean on_target);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* INCLUDED_PHASE_H */

