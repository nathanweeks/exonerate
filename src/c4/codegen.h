/****************************************************************\
*                                                                *
*  Code generation module for C4                                 *
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

#ifndef INCLUDED_CODEGEN_H
#define INCLUDED_CODEGEN_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <glib.h>
#include <stdio.h>

#include "argument.h"

typedef struct {
    gboolean use_compiled;
} Codegen_ArgumentSet;

Codegen_ArgumentSet *Codegen_ArgumentSet_create(Argument *arg);

gchar *Codegen_clean_path_component(gchar *name);

gboolean Codegen_file_exists(gchar *path);
gboolean Codegen_directory_exists(gchar *path);
/* FIXME: Use glib-2 functions instead when changing over */

typedef struct {
     FILE *fp;             /* FILE descriptor for the code */
    gchar *name;           /* Name of the file/function    */
    gchar *code_path;      /* Path to the .c file          */
    gchar *object_path;    /* Path to the .o file          */
     gint  indent;         /* Current indentation level    */
    gchar *func_prototype;
    gchar *func_return;
} Codegen;

Codegen *Codegen_create(gchar *directory, gchar *name);
   void  Codegen_destroy(Codegen *codegen);
/* If directory is NULL, the default is used ( ~/.c4_plugins )
 */

void  Codegen_indent(Codegen *codegen, gint indent_change);
void  Codegen_printf(Codegen *codegen, gchar *format, ...);
void  Codegen_compile(Codegen *c,
                      gchar *add_ccflags, gchar *add_ldflags);

/**/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* INCLUDED_CODEGEN_H */

