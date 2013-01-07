/****************************************************************\
*                                                                *
*  fastavalidcds: a utility to check for valid CDS sequences     *
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

#include <ctype.h> /* For toupper() */

#include "argument.h"
#include "fastadb.h"

static gboolean is_start_codon(gchar *seq){
    if((toupper(seq[0]) == 'A')
    && (toupper(seq[1]) == 'T')
    && (toupper(seq[2]) == 'G'))
        return TRUE;
    return FALSE;
    }

static gboolean is_stop_codon(gchar *seq){
    if(toupper(seq[0]) != 'T')
        return FALSE;
    if(toupper(seq[1]) == 'A'){
        if((toupper(seq[2]) != 'G') && (toupper(seq[2]) != 'A'))
            return FALSE;
    } else if(toupper(seq[1]) == 'G'){
        if(toupper(seq[2]) != 'A')
            return FALSE;
    } else {
        return FALSE;
        }
    return TRUE;
    }

static gboolean fasta_valid_cds_traverse_func(FastaDB_Seq *fdbs,
                                              gpointer user_data){
    register gint i;
    register gboolean explain = *((gboolean*)user_data);
    register gchar *seq;
    if(fdbs->seq->len % 3){
        if(explain)
            g_warning("%s length (%d) not divisible by 3",
                    fdbs->seq->id, fdbs->seq->len);
        return FALSE;
        }
    if(fdbs->seq->len < 6){
        if(explain)
            g_warning("%s length (%d) too short ( <2 codons )",
                    fdbs->seq->id, fdbs->seq->len);
        return FALSE;
        }
    seq = Sequence_get_str(fdbs->seq);
    if(!is_start_codon(seq)){
        if(explain)
            g_warning("%s missing start codon (has:%c%c%c)",
                    fdbs->seq->id,
                    seq[0], seq[1], seq[2]);
        g_free(seq);
        return FALSE;
        }
    if(!is_stop_codon(seq+fdbs->seq->len-3)){
        if(explain)
            g_warning("%s missing stop codon (has:%c%c%c)",
                    fdbs->seq->id,
                    seq[fdbs->seq->len-3],
                    seq[fdbs->seq->len-2],
                    seq[fdbs->seq->len-1]);
        g_free(seq);
        return FALSE;
        }
    for(i = 3; i < (fdbs->seq->len-3); i+=3)
        if(is_stop_codon(seq+i)){
            if(explain)
                g_warning("%s contains in-frame stop codon"
                          " (pos:%d, codon:%c%c%c)",
                    fdbs->seq->id, i, seq[i], seq[i+1], seq[i+2]);
            g_free(seq);
            return FALSE;
            }
    for(i = 0; i < fdbs->seq->len; i++){
        if((toupper(seq[i]) != 'A')
        && (toupper(seq[i]) != 'C')
        && (toupper(seq[i]) != 'G')
        && (toupper(seq[i]) != 'T')){
            if(explain){
                g_warning("%s contains non-ACGT base:"
                          " [pos: %d base: %c]",
                          fdbs->seq->id, i, seq[i]);
                }
            g_free(seq);
            return FALSE;
            }
        }
    FastaDB_Seq_print(fdbs, stdout, FastaDB_Mask_ID
                                   |FastaDB_Mask_DEF
                                   |FastaDB_Mask_SEQ);
    g_free(seq);
    return FALSE;
    }

int Argument_main(Argument *arg){
    register FastaDB *fdb;
    register ArgumentSet *as
           = ArgumentSet_create("Sequence Input Options");
    gchar *query_path;
    gboolean explain;
    ArgumentSet_add_option(as, 'f', "fasta", "path",
        "Fasta input file", NULL,
        Argument_parse_string, &query_path);
    ArgumentSet_add_option(as, 'e', "explain", NULL,
        "Explain reasons for rejection", "FALSE",
        Argument_parse_boolean, &explain);
    Argument_absorb_ArgumentSet(arg, as);
    Argument_process(arg, "fastavalidcds",
        "A utility to check for valid CDS sequences\n"
        "Guy St.C. Slater. guy@ebi.ac.uk. 2004.\n", NULL);
    fdb = FastaDB_open(query_path, NULL);
    FastaDB_traverse(fdb, FastaDB_Mask_ALL,
                         fasta_valid_cds_traverse_func, &explain);
    FastaDB_close(fdb);
    return 0;
    }

/**/

