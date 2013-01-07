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

#include <glib.h>
#include "codegen.h"

gint Argument_main(Argument *arg){
    register gchar *module = "testmodule";
    register Codegen *c = Codegen_create(NULL, module);
    Codegen_printf(c, "%s\n\n%s\n",
          "#include <stdio.h>",
          "int test_func(){");
    Codegen_indent(c, 1);
    Codegen_printf(c, "%s\n%s\n%s\n%s",
          "register int total = 1+2+3+4;",
          "printf(\"testing C4 plugin [%d]\\n\", total);",
          "return total;",
          "}");
    Codegen_indent(c, -1);
    Codegen_compile(c, NULL, NULL);
    Codegen_destroy(c);
    return 0;
    }

