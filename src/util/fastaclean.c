/****************************************************************\
*                                                                *
*  fastaclean : a utility to clean fasta format files            *
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

typedef struct {
    gboolean clean_protein;
    gboolean clean_acgtn;
} FastaClean_Info;

static gboolean fasta_clean_traverse_func(FastaDB_Seq *fdbs,
                                          gpointer user_data){
    register FastaClean_Info *fci = user_data;
    register Alphabet_Filter_Type filter_type;
    register Sequence *seq;
    if(fci->clean_acgtn)
        filter_type = Alphabet_Filter_Type_CLEAN_ACGTN;
    else
        filter_type = Alphabet_Filter_Type_CLEAN;
    seq = Sequence_filter(fdbs->seq, filter_type);
    Sequence_print_fasta(seq, stdout, FALSE);
    Sequence_destroy(seq);
    return FALSE;
    }

int Argument_main(Argument *arg){
    register FastaDB *fdb;
    register ArgumentSet *as
           = ArgumentSet_create("Sequence Input Options");
    gchar *query_path;
    register Alphabet *alphabet;
    FastaClean_Info fci;
    ArgumentSet_add_option(as, 'f', "fasta", "path",
        "Fasta input file", NULL,
        Argument_parse_string, &query_path);
    ArgumentSet_add_option(as, 'p', "protein", NULL,
        "Clean protein database", "FALSE",
        Argument_parse_boolean, &fci.clean_protein);
    ArgumentSet_add_option(as, 'a', "acgtn", NULL,
        "Only allow [ACGTN] nucleotide symbols", "FALSE",
        Argument_parse_boolean, &fci.clean_acgtn);
    Argument_absorb_ArgumentSet(arg, as);
    Argument_process(arg, "fastaclean",
        "A utility to clean fasta format file symbols\n"
        "Guy St.C. Slater. guy@ebi.ac.uk. 2000-2003.\n", NULL);
    if(fci.clean_protein)
        alphabet = Alphabet_create(Alphabet_Type_PROTEIN, FALSE);
    else
        alphabet = Alphabet_create(Alphabet_Type_DNA, FALSE);
    fdb = FastaDB_open(query_path, alphabet);
    FastaDB_traverse(fdb, FastaDB_Mask_ALL, fasta_clean_traverse_func,
                     &fci);
    FastaDB_close(fdb);
    Alphabet_destroy(alphabet);
    return 0;
    }

