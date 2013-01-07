/****************************************************************\
*                                                                *
*  Library for command line argument processing                  *
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

#include "argument.h"

static gchar *file_argument_handler(gchar *arg_string, gpointer data){
    g_message("file handler received [%s]", arg_string);
    return NULL;
    }

int Argument_main(Argument *arg){
    register ArgumentSet *as_a = ArgumentSet_create("Arg set one");
    register ArgumentSet *as_b = ArgumentSet_create("Arg set two");
    register gint i;
    GPtrArray *argument_list;
    gboolean bool_A, bool_B;
    ArgumentSet_add_option(as_a, 'f', "file", "path name",
                           "Input file for blah or blah and blah",
                           "/dev/zero",
                           file_argument_handler,
                           NULL);
    ArgumentSet_add_option(as_a, 'o', "output", "output path",
                           "output or whatever",
                           "/dev/null",
                           file_argument_handler,
                           NULL);
    ArgumentSet_add_option(as_a, 'A', "bool_A", NULL,
                           "whatever",
                           "FALSE",
                           Argument_parse_boolean, &bool_A);
    ArgumentSet_add_option(as_a, 'B', "bool_B", NULL,
                           "whatever",
                           "FALSE",
                           Argument_parse_boolean, &bool_B);
    ArgumentSet_add_option(as_a, 'l', "list", "path list",
                           "list of paths or something",
                           "one",
                           NULL, &argument_list);
    Argument_absorb_ArgumentSet(arg, as_a);
    Argument_absorb_ArgumentSet(arg, as_b);
    Argument_process(arg, "blah",
        "A program for stuff\n"
        "Guy St.C. Slater.  guy@ebi.ac.uk  December 2000.\n", NULL);
    for(i = 0; i < argument_list->len; i++)
        g_message("List [%d] = [%s]",
                i, (gchar*)argument_list->pdata[i]);
    return 0;
    }

