/****************************************************************\
*                                                                *
*  esd2esi - generate an exonerate index file from a dataset     *
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
#include "dataset.h"
#include "index.h"

int Argument_main(Argument *arg){
    register ArgumentSet *as
           = ArgumentSet_create("Input and Output Options");
    gchar *dataset_path, *index_path;
    register Dataset *dataset;
    register Index *index;
    gboolean is_translated = FALSE;
    register gint word_length;
    gint dna_word_length, protein_word_length,
         word_jump, word_ambiguity,
         saturate_threshold, memory_limit;
    /**/
    ArgumentSet_add_option(as, 'd', "dataset", "path",
        "Exonerate dataset file", NULL,
        Argument_parse_string, &dataset_path);
    ArgumentSet_add_option(as, 'o', "output", "path",
        "Output path for .esi file", NULL,
        Argument_parse_string, &index_path);
    ArgumentSet_add_option(as, '\0', "translate", NULL,
        "Translate the dataset for comparison with protein", "FALSE",
        Argument_parse_boolean, &is_translated);
    /**/
    /* FIXME: tidy duplicated arguments also in libraries */
    ArgumentSet_add_option(as, 0, "dnawordlen", "bp",
        "Wordlength for DNA words", "12",
        Argument_parse_int, &dna_word_length);
    ArgumentSet_add_option(as, 0, "proteinwordlen", "aa",
        "Wordlength for protein words", "5",
        Argument_parse_int, &protein_word_length);
    ArgumentSet_add_option(as, 0, "wordjump", NULL,
        "Jump between database words", "1",
        Argument_parse_int, &word_jump);
    ArgumentSet_add_option(as, 0, "wordambiguity", NULL,
        "Number of ambiguous words to index", "1",
        Argument_parse_int, &word_ambiguity);
    ArgumentSet_add_option(as, 0, "saturatethreshold", NULL,
        "Word saturation threshold", "10",
        Argument_parse_int, &saturate_threshold);
    /**/
    ArgumentSet_add_option(as, 0, "memorylimit", NULL,
        "Memory limit for database indexing", "1024",
        Argument_parse_int, &memory_limit);
    /**/
    Argument_absorb_ArgumentSet(arg, as);
    Argument_process(arg, "esd2esi",
        "generate an exonerate sequence index file\n"
        "Guy St.C. Slater. guy@ebi.ac.uk. 2007.\n", NULL);
    dataset = Dataset_read(dataset_path);
    word_length = (is_translated
                  || (dataset->alphabet->type == Alphabet_Type_PROTEIN))
                ? protein_word_length
                : dna_word_length;
    if(word_ambiguity < 1)
        g_error("Word ambiguity cannot be less than one.");
    if((word_ambiguity > 1)
    && (dataset->alphabet->type == Alphabet_Type_PROTEIN))
        g_error("Protein ambuigity symbols not implemented");
    g_message("Building index");
    index = Index_create(dataset, is_translated, word_length,
                         word_jump, word_ambiguity,
                         saturate_threshold, index_path, dataset_path,
                         memory_limit);
    Index_destroy(index);
    Dataset_destroy(dataset);
    g_message("-- completed");
    return 0;
    }

/**/

