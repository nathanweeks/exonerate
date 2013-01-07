/****************************************************************\
*                                                                *
*  edit distance model                                           *
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

#ifndef INCLUDED_EDIT_DISTANCE_H
#define INCLUDED_EDIT_DISTANCE_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "c4.h"

typedef struct {
    gchar *query;
    gchar *target;
     gint  query_len;
     gint  target_len;
} EditDistance_Data;

EditDistance_Data *EditDistance_Data_create(gchar *query,
                                            gchar *target);
void EditDistance_Data_destroy(EditDistance_Data *edd);

C4_Model *EditDistance_create(void);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* INCLUDED_EDIT_DISTANCE_H */

