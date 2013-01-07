/****************************************************************\
*                                                                *
*  Module for frameshift modelling                               *
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

#ifndef INCLUDED_FRAMESHIFT_H
#define INCLUDED_FRAMESHIFT_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "c4.h"
#include "argument.h"

typedef struct {
    C4_Score frameshift_penalty;
} Frameshift_ArgumentSet;

typedef struct {
    Frameshift_ArgumentSet *fas;
} Frameshift_Data;

Frameshift_Data *Frameshift_Data_create(void);
           void  Frameshift_Data_destroy(Frameshift_Data *fd);

Frameshift_ArgumentSet *Frameshift_ArgumentSet_create(Argument *arg);

void Frameshift_add(C4_Model *model, C4_State *match_state,
                    gchar *suffix, gboolean apply_to_query);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* INCLUDED_FRAMESHIFT_H */

