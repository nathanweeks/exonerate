/****************************************************************\
*                                                                *
*  fastaannotatecdna : a utility to annotate cDNA with a protein *
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

#include <string.h> /* For strstr() */

#include "argument.h"
#include "fastadb.h"
#include "translate.h"

static gint find_translation(gchar *protein_str, gint protein_len,
                             Sequence *cdna, Translate *translate){
    register Sequence *trans_seq;
    register gchar *trans_str, *result;
    register gint i, pos, count = 0;
    for(i = 1; i <= 3; i++){
        trans_seq = Sequence_translate(cdna, translate, i);
        trans_str = Sequence_get_str(trans_seq);
        result = trans_str;
        while((result = strstr(result, protein_str))){
            pos = result-trans_str;
            g_print("annotation: %s %c %d %d\n",
                    cdna->id,
                    Sequence_get_strand_as_char(cdna),
                    (pos*3)+i,
                    protein_len*3);
            count++;
            result++;
            }
        g_free(trans_str);
        Sequence_destroy(trans_seq);
        }
    return count;
    }

static gint fastaannotatecdna(gchar *cdna_path, gchar *protein_path){
    register FastaDB *cdna_fdb, *protein_fdb;
    register FastaDB_Seq *cdna = NULL, *protein = NULL;
    register gint ret_val = 0, total;
    register FastaDB_Mask mask = FastaDB_Mask_ID|FastaDB_Mask_SEQ;
    register Sequence *rc_seq;
    register gchar *protein_str;
    register Translate *translate = Translate_create(FALSE);
    register Alphabet
        *cdna_alphabet = Alphabet_create(Alphabet_Type_DNA, FALSE),
        *protein_alphabet = Alphabet_create(Alphabet_Type_PROTEIN, FALSE);
    cdna_fdb = FastaDB_open(cdna_path, cdna_alphabet);
    protein_fdb = FastaDB_open(protein_path, protein_alphabet);
    while((cdna = FastaDB_next(cdna_fdb, mask))){
        protein = FastaDB_next(protein_fdb, mask);
        if(!protein){
            g_print("ERROR: fastaannotatecdna: %s: protein: %s is absent\n",
                    protein_path, cdna->seq->id);
            ret_val = 1;
            break;
            }
        if((protein->seq->len*3) > cdna->seq->len)
            g_print("ERROR: fastaannoatecdna: protein [%s](%d) too long for cdna [%s](%d)\n",
                     protein->seq->id, protein->seq->len,
                     cdna->seq->id, cdna->seq->len);
        protein_str = Sequence_get_str(protein->seq);
        /* forward */
        total = find_translation(protein_str, protein->seq->len,
                                 cdna->seq, translate);
        /* revcomp */
        rc_seq = Sequence_revcomp(cdna->seq);
        total += find_translation(protein_str, protein->seq->len,
                                  rc_seq, translate);
        if(total != 1){
            g_print("ERROR: fastaannoatecdna: Found %d locations for protein [%s] in [%s]\n",
                    total, protein->seq->id, cdna->seq->id);
            ret_val = 1;
            break;
            }
        Sequence_destroy(rc_seq);
        /**/
        g_free(protein_str);
        FastaDB_Seq_destroy(cdna);
        FastaDB_Seq_destroy(protein);
        cdna = protein = NULL;
        }
    if(cdna)
        FastaDB_Seq_destroy(cdna);
    if(protein)
        FastaDB_Seq_destroy(protein);
    if(ret_val == 0){
        protein = FastaDB_next(protein_fdb, mask);
        if(protein){
            g_print("ERROR: fastaannoatecdna: %s: cdna: %s absent\n",
                    cdna_path, protein->seq->id);
            ret_val = 1;
            FastaDB_Seq_destroy(protein);
            }
        }
    FastaDB_close(cdna_fdb);
    FastaDB_close(protein_fdb);
    Translate_destroy(translate);
    Alphabet_destroy(cdna_alphabet);
    Alphabet_destroy(protein_alphabet);
    return ret_val;
    }

int Argument_main(Argument *arg){
    register ArgumentSet *as
           = ArgumentSet_create("Sequence Input Options");
    gchar *cdna_path, *protein_path;
    ArgumentSet_add_option(as, 'c', "cdna", "path",
        "cDNA input file", NULL,
        Argument_parse_string, &cdna_path);
    ArgumentSet_add_option(as, 'p', "protein", "path",
        "Protein input file", NULL,
        Argument_parse_string, &protein_path);
    Argument_absorb_ArgumentSet(arg, as);
    Translate_ArgumentSet_create(arg);
    Argument_process(arg, "fastaannotatecdna",
        "Annotate a cDNA using the corresponding protein\n"
        "Guy St.C. Slater. guy@ebi.ac.uk. 2000-2006.\n", NULL);
    return fastaannotatecdna(cdna_path, protein_path);
    }

