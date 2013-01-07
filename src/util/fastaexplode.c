/****************************************************************\
*                                                                *
*  fastaexplode : break a fasta file into individual sequences   *
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

#include "argument.h"
#include "fastadb.h"

static gboolean fasta_explode_traverse_func(FastaDB_Seq *fdbs,
                                            gpointer user_data){
    register gchar *dir_path = user_data;
    register gchar *output_path = g_strconcat(dir_path,
        G_DIR_SEPARATOR_S, fdbs->seq->id, ".fa", NULL);
    register FILE *fp = fopen(output_path, "r");
    if(fp){
        fclose(fp);
        g_error("File [%s] already exists", output_path);
        }
    fp = fopen(output_path, "w");
    if(!fp)
        g_error("Could not open [%s] to write output", output_path);
    FastaDB_Seq_print(fdbs, fp, FastaDB_Mask_ID
                               |FastaDB_Mask_DEF
                               |FastaDB_Mask_SEQ);
    g_free(output_path);
    fclose(fp);
    return FALSE;
    }

int Argument_main(Argument *arg){
    register FastaDB *fdb;
    register ArgumentSet *as
           = ArgumentSet_create("Sequence Input Options");
    gchar *query_path, *dir_path;
    ArgumentSet_add_option(as, 'f', "fasta", "path",
        "Fasta input file", NULL,
        Argument_parse_string, &query_path);
    ArgumentSet_add_option(as, 'd', "directory", "path",
        "Output file directory", ".",
        Argument_parse_string, &dir_path);
    Argument_absorb_ArgumentSet(arg, as);
    Argument_process(arg, "fastaexplode",
        "Split a fasta file up into individual sequences\n"
        "Guy St.C. Slater. guy@ebi.ac.uk. 2000-2003.\n", NULL);
    fdb = FastaDB_open(query_path, NULL);
    FastaDB_traverse(fdb, FastaDB_Mask_ALL,
                     fasta_explode_traverse_func, dir_path);
    FastaDB_close(fdb);
    return 0;
    }

