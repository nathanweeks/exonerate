/****************************************************************\
*                                                                *
*  fastanrdb : create non-redundant versions of fasta database   *
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

#include <stdlib.h>
#include <string.h> /* For strcmp() */
#include <strings.h> /* For strcasecmp() */
#include <ctype.h> /* For toupper() */

#include "argument.h"
#include "fastadb.h"

/**/

typedef struct {
           gint  checksum;
           gint  length;
    FastaDB_Key *key;
    FastaDB_Seq *seq;
} NRDB_Data;

static NRDB_Data *NRDB_Data_create(FastaDB_Seq *fdbs){
    register NRDB_Data *nd = g_new(NRDB_Data, 1);
    nd->length = fdbs->seq->len;
    nd->checksum = Sequence_checksum(fdbs->seq);
    nd->key = FastaDB_Seq_get_key(fdbs);
    nd->seq = NULL;
    return nd;
    }

static void NRDB_Data_destroy(NRDB_Data *nd){
    if(nd->key)
        FastaDB_Key_destroy(nd->key);
    if(nd->seq)
        FastaDB_Seq_destroy(nd->seq);
    g_free(nd);
    return;
    }

static int NRDB_Data_sort_checksum_function(const void *a,
                                            const void *b){
    register NRDB_Data **nd_a = (NRDB_Data**)a,
                       **nd_b = (NRDB_Data**)b;
    return (*nd_a)->checksum-(*nd_b)->checksum;
    }

/**/

typedef int (*NRDB_CompFunc)(const char *s1, const char *s2);

typedef struct {
         GPtrArray *nrdb_data_list;
          gboolean  ignore_case;
          gboolean  check_revcomp;
     NRDB_CompFunc  comp;
} NRDB_Info;

static gboolean fasta_nrdb_traverse_func(FastaDB_Seq *fdbs,
                                         gpointer user_data){
    register NRDB_Info *ni = user_data;
    register NRDB_Data *nd = NRDB_Data_create(fdbs);
    register FastaDB_Seq *revcomp_fdbs;
    register NRDB_Data *revcomp_nd;
    register gchar *seq, *rseq;
    g_ptr_array_add(ni->nrdb_data_list, nd);
    if(ni->check_revcomp){
        revcomp_fdbs = FastaDB_Seq_revcomp(fdbs);
        /* Check against palindromes */
        seq = Sequence_get_str(fdbs->seq);
        rseq = Sequence_get_str(revcomp_fdbs->seq);
        if(ni->comp(seq, rseq)){
            revcomp_nd = NRDB_Data_create(revcomp_fdbs);
            g_ptr_array_add(ni->nrdb_data_list, revcomp_nd);
            }
        FastaDB_Seq_destroy(revcomp_fdbs);
        g_free(seq);
        g_free(rseq);
        }
    return FALSE;
    }

/**/

static void NRDB_Data_report_redundant_set(GPtrArray *redundant_set){
    register gint i;
    register GPtrArray *forward_list = g_ptr_array_new(),
                       *reverse_list = g_ptr_array_new();
    register FastaDB_Seq *fdbs, *first_forward = NULL;
    register GString *merge_def;
    register gchar *curr_def;
    register FastaDB_Mask mask = FastaDB_Mask_ID
                               | FastaDB_Mask_DEF
                               | FastaDB_Mask_SEQ;
    for(i = 0; i < redundant_set->len; i++){
        fdbs = redundant_set->pdata[i];
        if(fdbs->seq->strand == Sequence_Strand_REVCOMP){
            g_ptr_array_add(reverse_list, fdbs);
        } else {
            if(first_forward){
                g_ptr_array_add(forward_list, fdbs);
            } else {
                first_forward = fdbs;
                }
            }
        }
    if(first_forward && ((forward_list->len+1) >= reverse_list->len)){
        curr_def = first_forward->seq->def;
        merge_def = g_string_sized_new(64);
        for(i = 0; i < forward_list->len; i++){
            fdbs = forward_list->pdata[i];
            g_string_append_c(merge_def, ' ');
            g_string_append(merge_def, fdbs->seq->id);
            }
        for(i = 0; i < reverse_list->len; i++){
            fdbs = reverse_list->pdata[i];
            g_string_append_c(merge_def, ' ');
            g_string_append(merge_def, fdbs->seq->id);
            g_string_append(merge_def, ".revcomp");
            }
        first_forward->seq->def = merge_def->str;
        FastaDB_Seq_print(first_forward, stdout, mask);
        g_string_free(merge_def, TRUE);
        first_forward->seq->def = curr_def;
        }
    g_ptr_array_free(forward_list, TRUE);
    g_ptr_array_free(reverse_list, TRUE);
    return;
    }

static void NRDB_Data_merge_and_print(NRDB_Info *ni, FastaDB *fdb){
    register gint i, j;
    register NRDB_Data *nd_a, *nd_b;
    register FastaDB_Seq *fdbs_a, *fdbs_b;
    register GPtrArray *redundant_set = g_ptr_array_new();
    register gchar *seq_a, *seq_b;
    for(i = 0; i < ni->nrdb_data_list->len; i++){
        nd_a = ni->nrdb_data_list->pdata[i];
        if(!nd_a) /* Already used */
            continue;
        fdbs_a = FastaDB_Key_get_seq(nd_a->key, FastaDB_Mask_ALL);
        g_ptr_array_add(redundant_set, fdbs_a);
        for(j = i+1; j < ni->nrdb_data_list->len; j++){
            nd_b = ni->nrdb_data_list->pdata[j];
            if(!nd_b) /* Already used */
                continue;
            if(nd_a->checksum != nd_b->checksum)
                break;
            if(nd_a->length != nd_b->length)
                continue;
            if(!nd_b->seq){
                /* Not available from previous checksum collision */
                nd_b->seq = FastaDB_Key_get_seq(nd_b->key,
                                                FastaDB_Mask_ALL);
                }
            fdbs_b = nd_b->seq;
            seq_a = Sequence_get_str(fdbs_a->seq);
            seq_b = Sequence_get_str(fdbs_b->seq);
            if(!ni->comp(seq_a, seq_b)){
                g_ptr_array_add(redundant_set, FastaDB_Seq_share(fdbs_b));
                NRDB_Data_destroy(nd_b);
                ni->nrdb_data_list->pdata[j] = NULL;
                }
            g_free(seq_a);
            g_free(seq_b);
            }
        NRDB_Data_destroy(nd_a);
        ni->nrdb_data_list->pdata[i] = NULL;
        NRDB_Data_report_redundant_set(redundant_set);
        for(j = 0; j < redundant_set->len; j++)
            FastaDB_Seq_destroy(redundant_set->pdata[j]);
        g_ptr_array_set_size(redundant_set, 0);
        }
    return;
    }

/**/

int Argument_main(Argument *arg){
    register FastaDB *fdb;
    register ArgumentSet *as
           = ArgumentSet_create("Sequence Input Options");
    register Alphabet *alphabet;
    NRDB_Info ni;
    gchar *query_path;
    ArgumentSet_add_option(as, 'f', "fasta", "path",
        "Fasta input file", NULL,
        Argument_parse_string, &query_path);
    ArgumentSet_add_option(as, 'i', "ignorecase", NULL,
        "Ignore sequence case", "FALSE",
        Argument_parse_boolean, &ni.ignore_case);
    ArgumentSet_add_option(as, 'r', "revcomp", NULL,
        "Check for revcomp duplicates", "FALSE",
        Argument_parse_boolean, &ni.check_revcomp);
    Argument_absorb_ArgumentSet(arg, as);
    Argument_process(arg, "fastanrdb",
        "A utility create non-redundant fasta sequence database\n"
        "Guy St.C. Slater. guy@ebi.ac.uk. 2000-2003.\n", NULL);
    ni.nrdb_data_list = g_ptr_array_new();
    if(ni.ignore_case)
        ni.comp = strcasecmp;
    else
        ni.comp = strcmp;
    if(ni.check_revcomp)
        alphabet = Alphabet_create(Alphabet_Type_DNA, FALSE);
    else
        alphabet = Alphabet_create(Alphabet_Type_UNKNOWN, FALSE);
    fdb = FastaDB_open(query_path, alphabet);
    FastaDB_traverse(fdb, FastaDB_Mask_ALL,
        fasta_nrdb_traverse_func, &ni);
    qsort(ni.nrdb_data_list->pdata, ni.nrdb_data_list->len,
          sizeof(gpointer), NRDB_Data_sort_checksum_function);
    NRDB_Data_merge_and_print(&ni, fdb);
    /* nrdb_data_list will be already cleared */
    g_ptr_array_free(ni.nrdb_data_list, TRUE);
    FastaDB_close(fdb);
    Alphabet_destroy(alphabet);
    return 0;
    }

