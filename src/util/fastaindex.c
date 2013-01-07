/****************************************************************\
*                                                                *
*  fastaindex : index a fasta format file                        *
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "argument.h"
#include "fastadb.h"

static int Fasta_Index_compare(const void *a, const void *b){
    register FastaDB_Seq **fdbs_a = (FastaDB_Seq**)a,
                         **fdbs_b = (FastaDB_Seq**)b;
    register gint result = strcmp((*fdbs_a)->seq->id,
                                  (*fdbs_b)->seq->id);
    if(!result)
        g_error("Duplicate fasta identifier [%s]", (*fdbs_a)->seq->id);
    return result;
    }

static void fasta_index_build(gchar *fasta_path, gchar *index_path){
    register FILE *fp;
    register FastaDB_Seq **fdbs_list;
    register gint i;
    guint total = 0;
    fp = fopen(index_path, "r");
    if(fp){
        fclose(fp);
        g_error("Index path [%s] already exists", index_path);
        }
    fp = fopen(index_path, "w");
    if(!fp)
        g_error("Could not write fasta index to [%s]", index_path);
    fdbs_list = FastaDB_all(fasta_path, NULL, FastaDB_Mask_ID, &total);
    qsort(fdbs_list, total, sizeof(FastaDB_Seq*), Fasta_Index_compare);
    for(i = 0; i < total; i++){
        fprintf(fp, "%s ", fdbs_list[i]->seq->id);
        CompoundFile_Pos_print(fp, fdbs_list[i]->location->pos);
        fprintf(fp, "\n");
        }
    FastaDB_Seq_all_destroy(fdbs_list);
    fclose(fp);
    return;
    }

int Argument_main(Argument *arg){
    register ArgumentSet *as
           = ArgumentSet_create("Sequence Input Options");
    gchar *fasta_path, *index_path;
    ArgumentSet_add_option(as, 'f', "fasta", "path",
        "Fasta input file", NULL,
        Argument_parse_string, &fasta_path);
    ArgumentSet_add_option(as, 'i', "index", "path",
        "Index file path", NULL,
        Argument_parse_string, &index_path);
    Argument_absorb_ArgumentSet(arg, as);
    Argument_process(arg, "fastaindex",
        "Index a fasta file\n"
        "Guy St.C. Slater. guy@ebi.ac.uk. 2000-2003.\n", NULL);
    fasta_index_build(fasta_path, index_path);
    return 0;
    }

