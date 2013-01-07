/****************************************************************\
*                                                                *
*  Analysis module for exonerate                                 *
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

#include "analysis.h"

int Argument_main(Argument *arg){
    register Analysis *analysis;
    register gchar *query_path, *target_path;
    register Alphabet_Type query_type, target_type;
    register GPtrArray *query_path_list = g_ptr_array_new(),
                       *target_path_list = g_ptr_array_new();
    if(arg->argc == 5){
        query_path = arg->argv[1];
        target_path = arg->argv[2];
        query_type = (arg->argv[3][0] == 'd')?Alphabet_Type_DNA
                                             :Alphabet_Type_PROTEIN;
        target_type = (arg->argv[4][0] == 'd')?Alphabet_Type_DNA
                                              :Alphabet_Type_PROTEIN;
        g_ptr_array_add(query_path_list, query_path);
        g_ptr_array_add(target_path_list, target_path);
        g_message("Creating Analysis with [%s][%s]",
                   query_path, target_path);
        analysis = Analysis_create(query_path_list, query_type, 0, 0,
                                   target_path_list, target_type, 0, 0,
                                   3);
        Analysis_process(analysis);
        Analysis_destroy(analysis);
        g_ptr_array_free(query_path_list, TRUE);
        g_ptr_array_free(target_path_list, TRUE);
    } else {
        g_warning("Test [%s] does nothing without arguments", __FILE__);
        g_message("usage: <query_path> <target_path>\n"
         "                <query_type> <target_type> (type==[d|p])\n");
        }
    return 0;
    }

