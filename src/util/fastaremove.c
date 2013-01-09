/****************************************************************\
*                                                                *
*  fastaremove : remove sequences from a fasta format file       *
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

#include <stdlib.h> /* For exit() */
#include <string.h> /* For strcmp () */
#include <strings.h> /* For strcasecmp() */
#include <ctype.h>  /* For isspace() */

#include "argument.h"
#include "fastadb.h"

static gboolean fasta_remove_traverse_func(FastaDB_Seq *fdbs,
                                           gpointer user_data){
    register GTree *id_tree = user_data;
    if(!g_tree_lookup(id_tree, fdbs->seq->id))
        FastaDB_Seq_print(fdbs, stdout, FastaDB_Mask_ID
                                       |FastaDB_Mask_DEF
                                       |FastaDB_Mask_SEQ);
    return FALSE;
    }

static void id_tree_prepare(FILE *fp, GTree *id_tree,
                            GStringChunk *id_alloc){
    register gint ch;
    register GString *id = g_string_sized_new(64);
    register gchar *temp;
    while((ch = getc(fp)) != EOF){
        if(ch == '\n'){
            temp = g_string_chunk_insert(id_alloc, id->str);
            g_tree_insert(id_tree, temp, temp);
            g_string_truncate(id, 0);
        } else {
            if(!isspace(ch))
                g_string_append_c(id, ch);
            }
        }
    g_string_free(id, TRUE);
    return;
    }

static gint id_tree_comp_func(gconstpointer a, gconstpointer b){
    return strcmp((gchar*)a, (gchar*)b);
    }

int Argument_main(Argument *arg){
    register FastaDB *fdb;
    register GTree *id_tree = g_tree_new(id_tree_comp_func);
    register GStringChunk *id_alloc = g_string_chunk_new(4096);
    register FILE *remove_fp;
    register ArgumentSet *as
           = ArgumentSet_create("Sequence Input Options");
    gchar *query_path, *remove_path;
    ArgumentSet_add_option(as, 'f', "fasta", "path",
        "Fasta input file", NULL,
        Argument_parse_string, &query_path);
    ArgumentSet_add_option(as, 'r', "remove", "path | stdin",
        "List of sequences to remove", NULL,
        Argument_parse_string, &remove_path);
    Argument_absorb_ArgumentSet(arg, as);
    Argument_process(arg, "fastaremove",
        "Remove sequences from a fasta file\n"
        "Guy St.C. Slater. guy@ebi.ac.uk. 2000-2003.\n", NULL);
    if(!strcasecmp(remove_path, "stdin")){
        remove_fp = stdin;
    } else {
        remove_fp = fopen(remove_path, "r");
        if(!remove_fp) {
            fprintf(stderr, "Could not open list for removal [%s]",
                     remove_path);
            exit(EXIT_FAILURE);
        }
        }
    id_tree_prepare(remove_fp, id_tree, id_alloc);
    fclose(remove_fp);
    fdb = FastaDB_open(query_path, NULL);
    FastaDB_traverse(fdb, FastaDB_Mask_ALL,
                     fasta_remove_traverse_func, id_tree);
    FastaDB_close(fdb);
    g_tree_destroy(id_tree);
    g_string_chunk_free(id_alloc);
    return 0;
    }

