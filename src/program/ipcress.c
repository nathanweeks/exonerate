/****************************************************************\
*                                                                *
*  ipcress : in-silico PCR experiment simulation system          *
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

#include <glib.h>
#include <stdio.h>
#include <stdlib.h> /* For atoi() */
#include <string.h> /* For strlen() */

#include "argument.h"
#include "exonerate_util.h"
#include "pcr.h"
#include "fastadb.h"
#include "lineparse.h"

/*

Input Format:

<id> <primer_A> <primer_B> <min_product_len> <max_product_len>
GO9892 CAAAGCTCTGTCCCAAAAAA TCTCTGCCCCTTTTACAGTG 90 100

Output Format:

...AAAAAAAAAAAAAAAAA.........................# forward
   |||||||||||||||||-->
   AAAAAAAAAAAAAAAAA       AAAAAAAAAAAAAAA   # primers
                        <--|||||||||||||||
...........................AAAAAAAAAAAAAAA...# revcomp

ipcress: seq_id experiment_id product_length \  3
         {primer_5:[A|B] pos_5, mismatch_5} \ # 3
         {primer_3:[A|B] pos_3, mismatch_3} \ # 3
         description: normal, revcomp, single_A, single_B  # 1

*/

typedef struct {
        FILE *input_fp;
   LineParse *input_lp;
        gint  mismatches;
        gint  seed_length;
        gint  memory_limit;
    gboolean  display_pretty;
    gboolean  display_products;
} Ipcress;

typedef enum {
    Ipcress_Type_FORWARD,
    Ipcress_Type_REVCOMP,
    Ipcress_Type_SINGLE_A,
    Ipcress_Type_SINGLE_B
} Ipcress_Type;

static gchar *Ipcress_Type_get_description(Ipcress_Type ipcress_type){
    register gchar *description = NULL;
    switch(ipcress_type){
        case Ipcress_Type_FORWARD:
            description = "forward";
            break;
        case Ipcress_Type_REVCOMP:
            description = "revcomp";
            break;
        case Ipcress_Type_SINGLE_A:
            description = "single_A";
            break;
        case Ipcress_Type_SINGLE_B:
            description = "single_B";
            break;
        default:
            g_error("Unknown ipcress type [%d]", ipcress_type);
            break;
            }
    return description;
    }

static Ipcress_Type Ipcress_Type_create(PCR_Probe *probe_a,
                                        PCR_Probe *probe_b){
    register Ipcress_Type ipcress_type;
    register PCR_Primer *primer_a = probe_a->pcr_primer,
                        *primer_b = probe_b->pcr_primer;
    register PCR_Experiment *pcr_experiment = primer_a->pcr_experiment;
    if((primer_a == pcr_experiment->primer_a)
    && (primer_b == pcr_experiment->primer_b))
        ipcress_type = Ipcress_Type_FORWARD;
    else
        ipcress_type = Ipcress_Type_REVCOMP;
    if(primer_a == primer_b){
        if(primer_a == pcr_experiment->primer_a)
            ipcress_type = Ipcress_Type_SINGLE_A;
        else
            ipcress_type = Ipcress_Type_SINGLE_B;
        }
    return ipcress_type;
    }

static gboolean Ipcress_report_product(Sequence *sequence,
                 PCR_Match *match_a, PCR_Match *match_b,
                 gint product_length, gpointer user_data){
    register PCR_Probe  *probe_a = match_a->pcr_probe,
                        *probe_b = match_b->pcr_probe;
    register PCR_Primer *primer_a = probe_a->pcr_primer,
                        *primer_b = probe_b->pcr_primer;
    register PCR_Experiment *pcr_experiment = primer_a->pcr_experiment;
    register Ipcress *ipcress = user_data;
    register gint i;
    register Ipcress_Type ipcress_type = Ipcress_Type_create(probe_a,
                                                             probe_b);
    register gchar *description, *s;
    register Sequence *subseq, *seq;
    g_assert(probe_a->strand == Sequence_Strand_FORWARD);
    g_assert(probe_b->strand == Sequence_Strand_REVCOMP);
    g_assert(pcr_experiment == primer_b->pcr_experiment);
    description = Ipcress_Type_get_description(ipcress_type);
    if(ipcress->display_pretty){
        g_print("\n");
        g_print("Ipcress result\n");
        g_print("--------------\n");
        g_print(" Experiment: %s\n", pcr_experiment->id);
        g_print("    Primers: %c %c\n",
                (primer_a == pcr_experiment->primer_a)?'A':'B',
                (primer_b == pcr_experiment->primer_a)?'A':'B');
        g_print("     Target: %s%s%s\n",
                sequence->id,
                sequence->def?" ":"",
                sequence->def?sequence->def:"");
        g_print("    Matches: %d/%d %d/%d\n",
                primer_a->length-match_a->mismatch,
                primer_a->length,
                primer_b->length-match_b->mismatch,
                primer_b->length);
        g_print("    Product: %d bp (range %d-%d)\n",
                product_length,
                pcr_experiment->min_product_len,
                pcr_experiment->max_product_len);
        g_print("Result type: %s\n", description);
        g_print("\n");
        /**/
        s = Sequence_get_substr(sequence, match_a->position,
                                primer_a->length);
        g_print("...%s.......", s); /* Start of top line */
        g_free(s);
        for(i = 0; i < primer_b->length; i++)
            g_print(".");
        g_print("... # forward\n");
        g_print("   "); /* Start of upper match line */
        for(i = 0; i < primer_a->length; i++)
            if(primer_a->forward[i]
            == Sequence_get_symbol(sequence, match_a->position+i))
                g_print("|");
            else
                g_print(" ");
        g_print("-->\n");
        g_print("5'-%.*s-3' ", /* Start of primer line */
                primer_a->length, primer_a->forward);
        Sequence_reverse_in_place(primer_b->forward, primer_b->length);
        g_print("3'-%.*s-5' # primers\n",
                primer_b->length, primer_b->forward);
        Sequence_reverse_in_place(primer_b->forward, primer_b->length);
        g_print("   "); /* Start of lower match line */
        for(i = 0; i < primer_a->length; i++)
            g_print(" ");
        g_print("    <--");
        for(i = 0; i < primer_b->length; i++)
            if(primer_b->revcomp[i]
            == Sequence_get_symbol(sequence, match_b->position+i))
                g_print("|");
            else
                g_print(" ");
        g_print("\n");
        g_print("..."); /* Start of base line */
        for(i = 0; i < primer_a->length; i++)
            g_print(".");
        g_print(".......");
        seq = Sequence_filter(sequence, Alphabet_Filter_Type_COMPLEMENT);
        s = Sequence_get_substr(seq, match_b->position,
                                primer_b->length);
        g_print("%s", s);
        Sequence_destroy(seq);
        g_free(s);
        g_print("... # revcomp\n");
        g_print("--\n");
        }
    g_print("ipcress: %s %s %d %c %d %d %c %d %d %s\n",
        sequence->id,
        pcr_experiment->id,
        product_length,
        /**/
        (primer_a == pcr_experiment->primer_a)?'A':'B',
        match_a->position,
        match_a->mismatch,
        /**/
        (primer_b == pcr_experiment->primer_a)?'A':'B',
        match_b->position,
        match_b->mismatch,
        /**/
        description);
    if(ipcress->display_products){
        subseq = Sequence_subseq(sequence, match_a->position,
                                 product_length);
        /* Should report product on the forward strand */
        if(ipcress_type == Ipcress_Type_REVCOMP)
            seq = Sequence_revcomp(subseq);
        else
            seq = Sequence_share(subseq);
        Sequence_destroy(subseq);
        g_print(">%s_product_%d seq %s start %d length %d\n",
                  pcr_experiment->id,
                  ++pcr_experiment->product_count,
                  sequence->id,
                  match_a->position,
                  product_length);
        Sequence_print_fasta_block(seq, stdout);
        Sequence_destroy(seq);
        }
    return FALSE;
    }

static Ipcress *Ipcress_create(gchar *ipcress_path, gint mismatches,
                               gint seed_length, gint memory_limit,
                               gboolean display_pretty,
                               gboolean display_products){
    register Ipcress *ipcress = g_new(Ipcress, 1);
    ipcress->input_fp = fopen(ipcress_path, "r");
    if(!ipcress->input_fp)
        g_error("Could not open ipcress file [%s]", ipcress_path);
    ipcress->input_lp = LineParse_create(ipcress->input_fp);
    ipcress->mismatches = mismatches;
    ipcress->seed_length = seed_length;
    ipcress->memory_limit = memory_limit;
    ipcress->display_pretty = display_pretty;
    ipcress->display_products = display_products;
    return ipcress;
    }

static PCR *Ipcress_prepare_PCR(Ipcress *ipcress){
    register PCR *pcr = PCR_create(Ipcress_report_product, ipcress,
                                   ipcress->mismatches,
                                   ipcress->seed_length);
    register gsize memory_usage;
    register gint experiment_count = 0, word_count;
    register gchar *id, *primer_a, *primer_b;
    register gint min_product_len, max_product_len;
    /* Read as many experiments as possible into memory */
    while((word_count = LineParse_word(ipcress->input_lp)) != EOF){
        if(word_count == 0)
            continue;
        if(word_count != 5)
            g_error("Bad line in ipcress file with [%d] fields",
                    word_count);
        id = ipcress->input_lp->word->pdata[0];
        primer_a = ipcress->input_lp->word->pdata[1];
        primer_b = ipcress->input_lp->word->pdata[2];
        strup(primer_a);
        strup(primer_b);
        min_product_len = atoi(ipcress->input_lp->word->pdata[3]);
        max_product_len = atoi(ipcress->input_lp->word->pdata[4]);
        memory_usage = PCR_add_experiment(pcr, id, primer_a, primer_b,
                                      min_product_len, max_product_len);
        experiment_count++;
        if(memory_usage > (ipcress->memory_limit << 20)){
            break;
            }
        }
    if(!experiment_count){
        PCR_destroy(pcr);
        return NULL;
        }
    g_message("Loaded [%d] experiments", experiment_count);
    PCR_prepare(pcr);
    return pcr;
    }

static void Ipcress_analyse(Ipcress *ipcress,
                            GPtrArray *sequence_path_list){
    register PCR *pcr;
    register Alphabet *alphabet = Alphabet_create(Alphabet_Type_DNA,
                                                  FALSE);
    register FastaDB *fdb = FastaDB_open_list(sequence_path_list,
                                              alphabet);
    register FastaDB_Seq *fdbs;
    register Sequence *seq;
    /* While not at the end of the input file */
    while((pcr = Ipcress_prepare_PCR(ipcress))){
        /* Analyse each sequence in DB */
        while((fdbs = FastaDB_next(fdb, FastaDB_Mask_ALL))){
            seq = Sequence_filter(fdbs->seq, Alphabet_Filter_Type_UNMASKED);
            PCR_simulate(pcr, seq);
            Sequence_destroy(seq);
            FastaDB_Seq_destroy(fdbs);
            }
        FastaDB_rewind(fdb);
        PCR_destroy(pcr);
        }
    FastaDB_close(fdb);
    Alphabet_destroy(alphabet);
    return;
    }

static void Ipcress_destroy(Ipcress *ipcress){
    LineParse_destroy(ipcress->input_lp);
    fclose(ipcress->input_fp);
    g_free(ipcress);
    return;
    }

int Argument_main(Argument *arg){
    gchar *ipcress_path;
    GPtrArray *sequence_path_list;
    gint mismatches, memory_limit, seed_length;
    gboolean display_pretty, display_products;
    /**/
    register ArgumentSet *as_input
           = ArgumentSet_create("File input options");
    register ArgumentSet *as_parameter
           = ArgumentSet_create("Ipcress parameters");
    register Ipcress *ipcress;
    /**/
    ArgumentSet_add_option(as_input, 'i', "input", "ipcress file",
          "Primer data in IPCRESS file format", NULL,
          Argument_parse_string, &ipcress_path);
    ArgumentSet_add_option(as_input, 's', "sequence", "fasta paths",
          "Fasta format sequence database", NULL,
          NULL, &sequence_path_list);
    ArgumentSet_add_option(as_parameter, 'm', "mismatch",
          "mismatches",
          "number of mismatches allowed per primer", "0",
          Argument_parse_int, &mismatches);
    ArgumentSet_add_option(as_parameter, 'M', "memory", "Mb",
          "Memory limit for FSM data", "32",
          Argument_parse_int, &memory_limit);
    ArgumentSet_add_option(as_parameter, 'p', "pretty", NULL,
          "Include \'pretty\' output", "TRUE",
          Argument_parse_boolean, &display_pretty);
    ArgumentSet_add_option(as_parameter, 'S', "seed", NULL,
          "Seed length (use zero for full length)", "12",
          Argument_parse_int, &seed_length);
    ArgumentSet_add_option(as_parameter, 'P', "products", NULL,
          "Report PCR products", "FALSE",
          Argument_parse_boolean, &display_products);
    Argument_absorb_ArgumentSet(arg, as_input);
    Argument_absorb_ArgumentSet(arg, as_parameter);
    Argument_process(arg, "ipcress",
        "In-silico PCR experiment simulation system\n"
        "Guy St.C. Slater.  guy@ebi.ac.uk  June 2001.\n", NULL);
    ipcress = Ipcress_create(ipcress_path, mismatches, seed_length,
                             memory_limit, display_pretty,
                             display_products);
    Ipcress_analyse(ipcress, sequence_path_list);
    Ipcress_destroy(ipcress);
    g_print("-- completed ipcress analysis\n");
    return 0;
    }

