/****************************************************************\
*                                                                *
*  fastahardmask : convert lower case chars to 'N's              *
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

#include <ctype.h>

#include "argument.h"
#include "fastadb.h"

static gboolean fasta_hard_mask_traverse_func(FastaDB_Seq *fdbs,
                                              gpointer user_data){
    register Sequence *s = Sequence_filter(fdbs->seq,
                                           Alphabet_Filter_Type_MASKED);
    Sequence_print_fasta(s, stdout, FALSE);
    Sequence_destroy(s);
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
    Argument_process(arg, "fastahardmask",
        "A utility to convert lower fasta sequence symbols to Ns\n"
        "Guy St.C. Slater. guy@ebi.ac.uk. 2000-2003.\n", NULL);
    fdb = FastaDB_open(query_path, NULL);
    FastaDB_traverse(fdb, FastaDB_Mask_ALL,
                     fasta_hard_mask_traverse_func, NULL);
    FastaDB_close(fdb);
    return 0;
    }

