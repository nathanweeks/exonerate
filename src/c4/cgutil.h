/****************************************************************\
*                                                                *
*  Utilities to facilitate code generation DP implementations    *
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

#ifndef INCLUDED_CGUTIL_H
#define INCLUDED_CGUTIL_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <glib.h>
#include <stdio.h>

#include "c4.h"
#include "codegen.h"

void CGUtil_print_header(Codegen *codegen, C4_Model *model);
void CGUtil_print_footer(Codegen *codegen);
void CGUtil_compile(Codegen *codegen, C4_Model *model);

void CGUtil_prep(Codegen *codegen, C4_Model *model, gboolean is_init);

/*
void CGUtil_model_init(Codegen *codegen);
void CGUtil_model_end(Codegen *codegen);
*/


/**/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* INCLUDED_CODEGEN_H */

