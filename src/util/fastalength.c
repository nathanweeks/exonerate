/****************************************************************\
*                                                                *
*  fastalength : a utility to print fasta sequence length        *
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
#include "fastadb.h"


static gboolean fasta_length_traverse_func(FastaDB_Seq *fdbs,
                                           gpointer user_data){
    g_print("%d %s\n", fdbs->seq->len, fdbs->seq->id);
    return FALSE;
    }

int Argument_main(Argument *arg){
    register FastaDB *fdb;
    register ArgumentSet *as
           = ArgumentSet_create("Sequence Input Options");
    gchar *query_path;
    ArgumentSet_add_option(as, 'f', "fasta", "path",
        "Fasta input file", NULL,
        Argument_parse_string, &query_path);
    Argument_absorb_ArgumentSet(arg, as);
    Argument_process(arg, "fastalength",
        "A utility to report fasta sequence lengths\n"
        "Guy St.C. Slater. guy@ebi.ac.uk. 2000-2003.\n", NULL);
    fdb = FastaDB_open(query_path, NULL);
    FastaDB_traverse(fdb, FastaDB_Mask_ID|FastaDB_Mask_LEN,
                     fasta_length_traverse_func, NULL);
    FastaDB_close(fdb);
    return 0;
    }

