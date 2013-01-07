/****************************************************************\
*                                                                *
*  Library for reading large and/or split files                  *
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

#include <stdio.h>
#include "compoundfile.h"

int main(void){
    register gchar *path_a =
        g_strconcat(SOURCE_ROOT_DIR, G_DIR_SEPARATOR_S,
                    "src", G_DIR_SEPARATOR_S,
                    "general", G_DIR_SEPARATOR_S,
                    "compoundfile.h", NULL);
    register gchar *path_b =
        g_strconcat(SOURCE_ROOT_DIR, G_DIR_SEPARATOR_S,
                    "src", G_DIR_SEPARATOR_S,
                    "general", G_DIR_SEPARATOR_S,
                    "compoundfile.c", NULL);
    register gchar *path_c =
        g_strconcat(SOURCE_ROOT_DIR, G_DIR_SEPARATOR_S,
                    "src", G_DIR_SEPARATOR_S,
                    "general", G_DIR_SEPARATOR_S,
                    "compoundfile.test.c", NULL);
    register CompoundFile *cf;
    register gint ch, total = 0;
    register GPtrArray *path_list = g_ptr_array_new();
    register CompoundFile_Pos total_length;
    register CompoundFile_Location *cfl_start, *cfl_stop;
    g_ptr_array_add(path_list, path_a);
    g_ptr_array_add(path_list, path_b);
    g_ptr_array_add(path_list, path_c);
    cf = CompoundFile_create(path_list, TRUE);
    g_message("Compound file:");
    while((ch = CompoundFile_getc(cf)) != EOF){
        g_print("%c", ch);
        total++;
        }
    g_message("--");
    g_message("Compound file again:");
    CompoundFile_rewind(cf);
    while((ch = CompoundFile_getc(cf)) != EOF)
        g_print("%c", ch);
    g_message("--");
    /* Read central half of compound file */
    total_length = CompoundFile_get_length(cf);
    cfl_start = CompoundFile_Location_from_pos(cf, total_length/4);
    cfl_stop = CompoundFile_Location_from_pos(cf, (total_length/4)*3);
    CompoundFile_set_limits(cf, cfl_start, cfl_stop);
    g_print("File truncated from [");
    CompoundFile_Pos_print(stdout, total_length);
    g_print("] to [");
    CompoundFile_Pos_print(stdout, CompoundFile_get_length(cf));
    g_print("]\n");
    CompoundFile_Location_destroy(cfl_start);
    CompoundFile_Location_destroy(cfl_stop);
    g_message("Middle half:");
    while((ch = CompoundFile_getc(cf)) != EOF)
        g_print("%c", ch);
    g_message("\n--");
    /**/
    CompoundFile_destroy(cf);
    g_ptr_array_free(path_list, TRUE);
    g_free(path_a);
    g_free(path_b);
    g_free(path_c);
    return 0;
    }

