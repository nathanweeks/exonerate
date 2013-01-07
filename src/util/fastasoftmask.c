/****************************************************************\
*                                                                *
*  fastasoftmask: merge masked and unmasked files to softmask    *
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

#include <string.h> /* For strcmp() */
#include <ctype.h>  /* For tolower() */

#include "argument.h"
#include "fastadb.h"

static void fasta_softmask_merge(Sequence *unmasked, Sequence *masked){
    register Sequence *softmask_seq;
    register gchar *sm = g_new(gchar, unmasked->len+1);
    register gint i;
    register gchar *ms = Sequence_get_str(unmasked),
                   *us = Sequence_get_str(masked);
    sm[unmasked->len] = '\0';
    for(i = 0; i < unmasked->len; i++){
        if((ms[i] == 'N') || (ms[i] == 'n')
        || (ms[i] == 'X') || (ms[i] == 'x')){
            sm[i] = tolower(us[i]);
        } else {
            sm[i] = us[i];
            }
        }
    softmask_seq = Sequence_create(unmasked->id, unmasked->def,
        sm, unmasked->len, Sequence_Strand_UNKNOWN, masked->alphabet);
    g_free(ms);
    g_free(us);
    g_free(sm);
    Sequence_print_fasta(softmask_seq, stdout, FALSE);
    Sequence_destroy(softmask_seq);
    return;
    }

static void fasta_softmask(gchar *unmasked_path, gchar *masked_path){
    register Alphabet *alphabet_unmasked, *alphabet_masked;
    register FastaDB *fdb_unmasked, *fdb_masked;
    register FastaDB_Seq *fdbs_unmasked = NULL, *fdbs_masked = NULL;
    register FastaDB_Mask mask = FastaDB_Mask_ID
                                |FastaDB_Mask_DEF
                                |FastaDB_Mask_SEQ;
    alphabet_unmasked = Alphabet_create(Alphabet_Type_UNKNOWN, FALSE);
    alphabet_masked = Alphabet_create(Alphabet_Type_UNKNOWN, TRUE);
    fdb_unmasked = FastaDB_open(unmasked_path, alphabet_unmasked);
    fdb_masked = FastaDB_open(masked_path, alphabet_masked);
    while((fdbs_unmasked = FastaDB_next(fdb_unmasked, mask))){
        fdbs_masked = FastaDB_next(fdb_masked, mask);
        if(!fdbs_masked)
            g_error("Could not find masked to match [%s]",
                    fdbs_unmasked->seq->id);
        if(strcmp(fdbs_unmasked->seq->id, fdbs_masked->seq->id))
            g_error("ID mismatch unmasked:[%s] masked:[%s]",
                    fdbs_unmasked->seq->id,
                    fdbs_masked->seq->id);
        if(fdbs_unmasked->seq->len != fdbs_unmasked->seq->len)
            g_error("Length mismatch: unmasked:%s(%d) masked:%s(%d)\n",
                    fdbs_unmasked->seq->id,
                    fdbs_unmasked->seq->len,
                    fdbs_masked->seq->id,
                    fdbs_masked->seq->len);
        fasta_softmask_merge(fdbs_unmasked->seq, fdbs_masked->seq);
        FastaDB_Seq_destroy(fdbs_unmasked);
        FastaDB_Seq_destroy(fdbs_masked);
        }
    fdbs_masked = FastaDB_next(fdb_masked, mask);
    if(fdbs_masked)
        g_error("Could not find unmasked to match [%s]",
                fdbs_masked->seq->id);
    FastaDB_close(fdb_unmasked);
    FastaDB_close(fdb_masked);
    Alphabet_destroy(alphabet_unmasked);
    Alphabet_destroy(alphabet_masked);
    return;
    }

int Argument_main(Argument *arg){
    register ArgumentSet *as
           = ArgumentSet_create("Sequence Input Options");
    gchar *unmasked_path, *masked_path;
    ArgumentSet_add_option(as, 'u', "unmasked", "path",
        "First fasta input file", NULL,
        Argument_parse_string, &unmasked_path);
    ArgumentSet_add_option(as, 'm', "masked", "path",
        "Second fasta input file", NULL,
        Argument_parse_string, &masked_path);
    Argument_absorb_ArgumentSet(arg, as);
    Argument_process(arg, "fastasoftmask",
        "Merge masked and unmasked files as softmasked sequences\n"
        "Guy St.C. Slater. guy@ebi.ac.uk. 2000-2003.\n", NULL);
    fasta_softmask(unmasked_path, masked_path);
    return 0;
    }

