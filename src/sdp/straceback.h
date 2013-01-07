/****************************************************************\
*                                                                *
*  C4 dynamic programming library - Scheduler Traceback          *
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

#ifndef INCLUDED_STRACEBACK_H
#define INCLUDED_STRACEBACK_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <glib.h>

#include "recyclebin.h"
#include "c4.h"
#include "region.h"

/**/

typedef struct STraceback_Cell {
                      gint  ref_count;
             C4_Transition *transition;
                      gint  length;
    struct STraceback_Cell *prev;
} STraceback_Cell;

typedef struct {
          gint  ref_count;
      gboolean  is_forward;
    RecycleBin *cell_recycle;
      C4_Model *model;
} STraceback;

     STraceback *STraceback_create(C4_Model *model, gboolean is_forward);
     STraceback *STraceback_share(STraceback *straceback);
           void  STraceback_destroy(STraceback *straceback);
STraceback_Cell *STraceback_add(STraceback *straceback,
                                C4_Transition *transition, gint length,
                                STraceback_Cell *prev);

STraceback_Cell *STraceback_Cell_share(STraceback_Cell *cell);
           void  STraceback_Cell_destroy(STraceback_Cell *cell,
                                         STraceback *straceback);

typedef struct {
    C4_Transition *transition;
             gint  length;
} STraceback_Operation;

typedef struct {
    GPtrArray *operation_list;
} STraceback_List;

STraceback_List *STraceback_List_create(STraceback *straceback,
                                        STraceback_Cell *source);
           void  STraceback_List_destroy(STraceback_List *stlist);

/**/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* INCLUDED_STRACEBACK_H */

