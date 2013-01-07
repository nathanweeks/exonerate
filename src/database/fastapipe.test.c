/****************************************************************\
*                                                                *
*  Library for All-vs-All Fasta Database Comparisons             *
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

#include "fastapipe.h"

static void test_init_func(gpointer user_data){
    register gint *count = user_data;
    g_message("Initiate pipeline");
    (*count) = 0;
    return;
    }

static void test_prep_func(gpointer user_data){
    g_message("Prepare pipeline");
    return;
    }

static void test_term_func(gpointer user_data){
    g_message("Terminate pipeline");
    return;
    }

static gboolean test_query_func(FastaDB_Seq *fdbs,
                                gpointer user_data){
    register gint *count = user_data;
    g_message("Have query [%s]", fdbs->seq->id);
    if((++(*count)) < 2) /* Try to use 2 queries per pipeline */
        return FALSE;
    return TRUE;
    }

static gboolean test_target_func(FastaDB_Seq *fdbs,
                                 gpointer user_data){
    g_message("Have target [%s]", fdbs->seq->id);
    return TRUE;
    }

gint Argument_main(Argument *arg){
    register FastaPipe *fasta_pipe;
    register gchar *query_path, *target_path;
    register GPtrArray *query_path_list = g_ptr_array_new(),
                       *target_path_list = g_ptr_array_new();
    register Alphabet *alphabet;
    register FastaDB *query_fdb, *target_fdb;
    gint count;
    if(arg->argc == 3){
        query_path = arg->argv[1];
        target_path = arg->argv[2];
        g_ptr_array_add(query_path_list, query_path);
        g_ptr_array_add(target_path_list, target_path);
        g_message("Processing test FastaPipe with:\nQ:[%s]\nT:[%s]",
                  query_path, target_path);
        alphabet = Alphabet_create(Alphabet_Type_UNKNOWN, FALSE);
        query_fdb = FastaDB_open_list(query_path_list, alphabet);
        target_fdb = FastaDB_open_list(target_path_list, alphabet);
        fasta_pipe = FastaPipe_create(query_fdb, target_fdb,
                test_init_func, test_prep_func, test_term_func,
                test_query_func, test_target_func, FastaDB_Mask_ALL,
                FALSE, TRUE);
        Alphabet_destroy(alphabet);
        while(FastaPipe_process(fasta_pipe, &count)){
            g_message("Processing pipeline");
            }
        FastaPipe_destroy(fasta_pipe);
        FastaDB_close(query_fdb);
        FastaDB_close(target_fdb);
    } else {
        g_warning("Test [%s] does nothing without arguments", __FILE__);
        }
    g_ptr_array_free(query_path_list, TRUE);
    g_ptr_array_free(target_path_list, TRUE);
    return 0;
    }

