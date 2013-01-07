/****************************************************************\
*                                                                *
*  fastasubseq: extract subsequences from fasta format sequences *
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

int Argument_main(Argument *arg){
    register FastaDB_Seq *fdbs;
    register Sequence *subseq;
    register ArgumentSet *as
           = ArgumentSet_create("Sequence Input Options");
    gchar *query_path;
    gint subseq_start, subseq_length;
    ArgumentSet_add_option(as, 'f', "fasta", "path",
        "Fasta input file", NULL,
        Argument_parse_string, &query_path);
    ArgumentSet_add_option(as, 's', "start", "start",
        "Sequence start position", NULL,
        Argument_parse_int, &subseq_start);
    ArgumentSet_add_option(as, 'l', "length", "length",
        "Subsequence length", NULL,
        Argument_parse_int, &subseq_length);
    Argument_absorb_ArgumentSet(arg, as);
    Argument_process(arg, "fastasubseq",
        "A utility to extract fasta format subsequences\n"
        "Guy St.C. Slater. guy@ebi.ac.uk. 2000-2003.\n", NULL);
    fdbs = FastaDB_get_single(query_path, NULL);
    if(subseq_start < 0)
        g_error("Subsequence must start after sequence start");
    if(subseq_length < 0)
        g_error("Subsequence length must be greater than zero");
    if((subseq_start+subseq_length) > fdbs->seq->len)
        g_error("Subsequence must end before end of [%s](%d)",
                fdbs->seq->id, fdbs->seq->len);
    subseq = Sequence_subseq(fdbs->seq, subseq_start, subseq_length);
    Sequence_print_fasta(subseq, stdout, FALSE);
    Sequence_destroy(subseq);
    FastaDB_Seq_destroy(fdbs);
    return 0;
    }

