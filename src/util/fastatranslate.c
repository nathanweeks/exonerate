/****************************************************************\
*                                                                *
*  fastatranslate: a utility to translate fasta format sequences *
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
#include "translate.h"

static void fasta_translate_seq(FastaDB_Seq *fdbs,
                                Translate *translate, gint frame,
                                Alphabet *protein_alphabet){
    register Sequence *rc_seq, *aa_seq;
    if(!frame){
        fasta_translate_seq(fdbs, translate, -3, protein_alphabet);
        fasta_translate_seq(fdbs, translate, -2, protein_alphabet);
        fasta_translate_seq(fdbs, translate, -1, protein_alphabet);
        fasta_translate_seq(fdbs, translate,  1, protein_alphabet);
        fasta_translate_seq(fdbs, translate,  2, protein_alphabet);
        fasta_translate_seq(fdbs, translate,  3, protein_alphabet);
        return;
        }
    if(frame < 1){
        rc_seq = Sequence_revcomp(fdbs->seq);
        aa_seq = Sequence_translate(rc_seq, translate, -frame);
        Sequence_destroy(rc_seq);
    } else {
        aa_seq = Sequence_translate(fdbs->seq, translate, frame);
        }
    if(aa_seq->len)
        Sequence_print_fasta(aa_seq, stdout, FALSE);
    Sequence_destroy(aa_seq);
    return;
    }

int Argument_main(Argument *arg){
    register FastaDB *fdb;
    register FastaDB_Seq *fdbs;
    register ArgumentSet *as
           = ArgumentSet_create("Sequence Input Options");
    register Translate *translate;
    register Alphabet *dna_alphabet
           = Alphabet_create(Alphabet_Type_DNA, FALSE);
    register Alphabet *protein_alphabet
           = Alphabet_create(Alphabet_Type_PROTEIN, FALSE);
    gchar *query_path;
    gint frame;
    ArgumentSet_add_option(as, 'f', "fasta", "path",
        "Fasta input file", NULL,
        Argument_parse_string, &query_path);
    ArgumentSet_add_option(as, 'F', "frame", NULL,
        "Reading frame to translate", "0",
        Argument_parse_int, &frame);
    Argument_absorb_ArgumentSet(arg, as);
    Translate_ArgumentSet_create(arg);
    Argument_process(arg, "fastatranslate",
        "A utility to translate fasta format sequences\n"
        "Guy St.C. Slater. guy@ebi.ac.uk. 2000-2003.\n", NULL);
    fdb = FastaDB_open(query_path, dna_alphabet);
    translate = Translate_create(FALSE);
    while((fdbs = FastaDB_next(fdb, FastaDB_Mask_ALL))){
        fasta_translate_seq(fdbs, translate, frame, protein_alphabet);
        FastaDB_Seq_destroy(fdbs);
        }
    FastaDB_close(fdb);
    Translate_destroy(translate);
    Alphabet_destroy(dna_alphabet);
    Alphabet_destroy(protein_alphabet);
    return 0;
    }

