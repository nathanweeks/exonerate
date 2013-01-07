/****************************************************************\
*                                                                *
*  fastacomposition : a utility to dump sequence composition      *
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

#include <string.h> /* For memset() */
#include <ctype.h>  /* For toupper() */

#include "argument.h"
#include "fastadb.h"

static void fasta_composition_report(gchar *name, glong *count,
                                     gboolean ignore_case){
    register gint i;
    g_print("%s", name);
    if(ignore_case){
        for(i = 0; i < 256; i++)
            if(count[i]){
                if(islower(i)){
                    g_print(" %c %ld", i, count[tolower(i)]
                                        + count[toupper(i)]);
                } else {
                    if(!isupper(i))
                        g_print(" %c %ld", i, count[i]);
                    }
                }
    } else {
        for(i = 0; i < 256; i++)
            if(count[i])
                g_print(" %c %ld", i, count[i]);
        }
    g_print("\n");
    return;
    }

typedef struct {
       glong count[256];
    gboolean report_separately;
    gboolean ignore_case;
} FastaComposition_Data;

static gboolean fasta_composition_traverse_func(FastaDB_Seq *fdbs,
                                                gpointer user_data){
    register gint i;
    register FastaComposition_Data *fcd = user_data;
    register gchar *str = Sequence_get_str(fdbs->seq);
    if(fcd->report_separately)
        memset(fcd->count, 0, sizeof(gint)*256);
    for(i = 0; i < fdbs->seq->len; i++)
        fcd->count[(guchar)str[i]]++;
    g_free(str);
    if(fcd->report_separately)
        fasta_composition_report(fdbs->seq->id, fcd->count, fcd->ignore_case);
    return FALSE;
    }

int Argument_main(Argument *arg){
    register FastaDB *fdb;
    register ArgumentSet *as
           = ArgumentSet_create("Sequence Input Options");
    gchar *query_path;
    FastaComposition_Data fcd = {{0},FALSE, FALSE};
    ArgumentSet_add_option(as, 'f', "fasta", "path",
        "Fasta input file", NULL,
        Argument_parse_string, &query_path);
    ArgumentSet_add_option(as, 'i', "ignorecase", NULL,
        "Ignore sequence case", "FALSE",
        Argument_parse_boolean, &fcd.ignore_case);
    ArgumentSet_add_option(as, 's', "separate", NULL,
        "Report composition for each sequence separately", "FALSE",
        Argument_parse_boolean, &fcd.report_separately);
    Argument_absorb_ArgumentSet(arg, as);
    Argument_process(arg, "fastacomposition",
        "A utility to report sequence composition\n"
        "Guy St.C. Slater. guy@ebi.ac.uk. 2000-2003.\n", NULL);
    fdb = FastaDB_open(query_path, NULL);
    FastaDB_traverse(fdb, FastaDB_Mask_ID|FastaDB_Mask_SEQ,
                     fasta_composition_traverse_func, &fcd);
    if(!fcd.report_separately)
        fasta_composition_report(query_path, fcd.count, fcd.ignore_case);
    FastaDB_close(fdb);
    return 0;
    }

