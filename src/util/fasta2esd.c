/****************************************************************\
*                                                                *
*  fasta2esd - generate an exonerate sequence database file      *
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
#include "alphabet.h"
#include "fastadb.h"

int Argument_main(Argument *arg){
    register ArgumentSet *as
           = ArgumentSet_create("Sequence Input Options");
    GPtrArray *input_path_list;
    gchar *output_path;
    gboolean softmask_input;
    Alphabet_Type alphabet_type;
    register Dataset *dataset;
    register FILE *fp;
    ArgumentSet_add_option(as, 'f', "fasta", "path",
        "Fasta format input sequences", NULL,
        NULL, &input_path_list);
    ArgumentSet_add_option(as, 'o', "output", "path",
        "Output path for .esd file", NULL,
        Argument_parse_string, &output_path);
    ArgumentSet_add_option(as, '\0', "alphabet", NULL,
        "Specify sequence alphabet type", "unknown",
        Alphabet_Argument_parse_alphabet_type, &alphabet_type);
    ArgumentSet_add_option(as, 's', "softmask", NULL,
        "Treat input sequences as softmasked", "TRUE",
        Argument_parse_boolean, &softmask_input);
    Argument_absorb_ArgumentSet(arg, as);
    Argument_process(arg, "fasta2esd",
        "generate an exonerate sequence database file\n"
        "Guy St.C. Slater. guy@ebi.ac.uk. 2006.\n", NULL);
    /* Check output absent before creating dataset */
    fp = fopen(output_path, "r");
    if(fp)
        g_error("Output file [%s] already exists", output_path);
    /**/
    g_message("Creating Dataset (alphabet:[%s] softmasked:[%s])",
            Alphabet_Type_get_name(alphabet_type),
            softmask_input?"yes":"no");
    dataset = Dataset_create(input_path_list, alphabet_type, softmask_input);
    g_message("Writing Dataset");
    Dataset_write(dataset, output_path);
    Dataset_destroy(dataset);
    g_message("-- completed");
    return 0;
    }

