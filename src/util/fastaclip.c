/****************************************************************\
*                                                                *
*  fastaclip : clip terminal Ns from fasta format sequences      *
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

static gboolean fasta_clip_traverse_func(FastaDB_Seq *fdbs,
                                              gpointer user_data){
    register gint clip_char = GPOINTER_TO_INT(user_data);
    register gint i, start_clip, end_clip, clip_len;
    register Sequence *clip_seq;
    for(i = 0; i < fdbs->seq->len; i++)
        if(Sequence_get_symbol(fdbs->seq, i) != clip_char)
            break;
    start_clip = i;
    for(i = fdbs->seq->len-1; i >= 0; i--)
        if(Sequence_get_symbol(fdbs->seq, i) != clip_char)
            break;
    end_clip = fdbs->seq->len-i-1;
    clip_len = fdbs->seq->len-start_clip-end_clip;
    clip_seq = Sequence_subseq(fdbs->seq, start_clip, clip_len);
    Sequence_print_fasta(clip_seq, stdout, FALSE);
    Sequence_destroy(clip_seq);
    return FALSE;
    }

int Argument_main(Argument *arg){
    register FastaDB *fdb;
    register ArgumentSet *as
           = ArgumentSet_create("Sequence Input Options");
    gchar *query_path;
    gchar clip_char;
    ArgumentSet_add_option(as, 'f', "fasta", "path",
        "Fasta input file", NULL,
        Argument_parse_string, &query_path);
    ArgumentSet_add_option(as, 'c', "clip", "path",
        "Clip character", "N",
        Argument_parse_char, &clip_char);
    Argument_absorb_ArgumentSet(arg, as);
    Argument_process(arg, "fastaclip",
        "A utility clip terminal Ns from fasta format sequences\n"
        "Guy St.C. Slater. guy@ebi.ac.uk. 2000-2003.\n", NULL);
    fdb = FastaDB_open(query_path, NULL);
    FastaDB_traverse(fdb, FastaDB_Mask_ALL,
                     fasta_clip_traverse_func,
                     GINT_TO_POINTER((gint)clip_char));
    FastaDB_close(fdb);
    return 0;
    }

