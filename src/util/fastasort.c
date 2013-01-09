/****************************************************************\
*                                                                *
*  fastasort : a utility for sorting fasta format databases      *
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
#include <strings.h> /* For strcasemp() */
#include <ctype.h> /* For toupper() */

#include "argument.h"
#include "fastadb.h"

/**/

typedef enum {
    Sort_KeyMask_ID      = (1<<0),
    Sort_KeyMask_LEN     = (1<<1),
    Sort_KeyMask_SEQ     = (1<<2),
    Sort_KeyMask_FORWARD = (1<<3),
    Sort_KeyMask_REVERSE = (1<<4)
} Sort_KeyMask;

typedef struct {
    FastaDB_Key *key;
    union {
        gchar *id;
         gint  len;
        gchar *seq_prefix;
    } key_data;
} Sort_Data;

static Sort_Data *Sort_Data_create(FastaDB_Seq *fdbs,
                                   Sort_KeyMask sort_key_mask){
    register Sort_Data *sd = g_new(Sort_Data, 1);
    if(sort_key_mask & Sort_KeyMask_ID){
        sd->key_data.id = g_strdup(fdbs->seq->id);
    } else if(sort_key_mask & Sort_KeyMask_LEN){
        sd->key_data.len = fdbs->seq->len;
    } else if(sort_key_mask & Sort_KeyMask_SEQ){
        sd->key_data.seq_prefix = Sequence_get_substr(fdbs->seq, 0, 16);
        }
    sd->key = FastaDB_Seq_get_key(fdbs);
    return sd;
    }

static void Sort_Data_destroy(Sort_Data *sd,
                              Sort_KeyMask sort_key_mask){
    if(sort_key_mask & Sort_KeyMask_ID){
        g_free(sd->key_data.id);
    } else if(sort_key_mask & Sort_KeyMask_SEQ){
        g_free(sd->key_data.seq_prefix);
        }
    FastaDB_Key_destroy(sd->key);
    g_free(sd);
    return;
    }

typedef int (*FastaSort_CompFunc)(const void *a, const void *b);

static int FastaSort_compare_id_forward(const void *a, const void *b){
    register Sort_Data **sd_a = (Sort_Data**)a,
                       **sd_b = (Sort_Data**)b;
    register gint ret_val = strcmp((*sd_a)->key_data.id,
                                   (*sd_b)->key_data.id);
    if(!ret_val)
        g_error("Duplicate key [%s]", (*sd_a)->key_data.id);
    return ret_val;
    }

static int FastaSort_compare_id_reverse(const void *a, const void *b){
    return FastaSort_compare_id_forward(b, a);
    }

static int FastaSort_compare_len_forward(const void *a, const void *b){
    register Sort_Data **sd_a = (Sort_Data**)a,
                       **sd_b = (Sort_Data**)b;
    return (*sd_a)->key_data.len - (*sd_b)->key_data.len;
    }

static int FastaSort_compare_len_reverse(const void *a, const void *b){
    return FastaSort_compare_len_forward(b, a);
    }

static int FastaSort_compare_seq_forward(const void *a, const void *b){
    register Sort_Data **sd_a = (Sort_Data**)a,
                       **sd_b = (Sort_Data**)b;
    register gint ret_val = strcmp((*sd_a)->key_data.seq_prefix,
                                   (*sd_b)->key_data.seq_prefix);
    register FastaDB_Seq *fdbs_a, *fdbs_b;
    register gchar *seq_a, *seq_b;
    if(!ret_val){
        fdbs_a = FastaDB_Key_get_seq((*sd_a)->key, FastaDB_Mask_SEQ);
        fdbs_b = FastaDB_Key_get_seq((*sd_b)->key, FastaDB_Mask_SEQ);
        seq_a = Sequence_get_str(fdbs_a->seq);
        seq_b = Sequence_get_str(fdbs_b->seq);
        FastaDB_Seq_destroy(fdbs_a);
        FastaDB_Seq_destroy(fdbs_b);
        ret_val = strcmp(seq_a, seq_b);
        g_free(seq_a);
        g_free(seq_b);
        }
    return ret_val;
    }

static int FastaSort_compare_seq_reverse(const void *a, const void *b){
    return FastaSort_compare_seq_forward(b, a);
    }

/* FIXME: optimisation
 *      : use a sequence cache for faster sequence sorting
 */

static Sort_KeyMask SortKey_Mask_create(gchar *sort_key,
                                        gboolean reverse_order){
    register Sort_KeyMask mask = 0;
    if(reverse_order)
        mask |= Sort_KeyMask_REVERSE;
    else
        mask |= Sort_KeyMask_FORWARD;
    if(!strcasecmp(sort_key, "id"))
        mask |= Sort_KeyMask_ID;
    else if(!strcasecmp(sort_key, "len"))
        mask |= Sort_KeyMask_LEN;
    else if(!strcasecmp(sort_key, "seq"))
        mask |= Sort_KeyMask_SEQ;
    else
        g_error("Unknown sort key [%s]", sort_key);
    return mask;
    }

static FastaSort_CompFunc fasta_sort_get_compare_func(
                          Sort_KeyMask sort_key_mask){
    switch((gint)sort_key_mask){
        case Sort_KeyMask_ID
            |Sort_KeyMask_FORWARD:
            return FastaSort_compare_id_forward;
            break;
        case Sort_KeyMask_ID
            |Sort_KeyMask_REVERSE:
            return FastaSort_compare_id_reverse;
            break;
        case Sort_KeyMask_LEN
            |Sort_KeyMask_FORWARD:
            return FastaSort_compare_len_forward;
            break;
        case Sort_KeyMask_LEN
            |Sort_KeyMask_REVERSE:
            return FastaSort_compare_len_reverse;
            break;
        case Sort_KeyMask_SEQ
            |Sort_KeyMask_FORWARD:
            return FastaSort_compare_seq_forward;
            break;
        case Sort_KeyMask_SEQ
            |Sort_KeyMask_REVERSE:
            return FastaSort_compare_seq_reverse;
            break;
        }
    g_error("Unknown sort key mask [%d]", sort_key_mask);
    return NULL;
    }

typedef struct {
             Sort_Data *prev_sd;
    FastaSort_CompFunc  comp_func;
          Sort_KeyMask  sort_key_mask;
} Sort_Check_Info;

static gboolean fasta_sort_check_traverse_func(FastaDB_Seq *fdbs,
                                               gpointer user_data){
    register Sort_Check_Info *sci = user_data;
    register gint ret_val;
    Sort_Data *sd = Sort_Data_create(fdbs, sci->sort_key_mask);
    if(sci->prev_sd){
        ret_val = sci->comp_func(&sci->prev_sd, &sd);
        if(ret_val > 1){
            g_print("File is not sorted: ");
            if(sci->sort_key_mask & Sort_KeyMask_ID){
                g_print("id [%s] followed by [%s]",
                    sci->prev_sd->key_data.id, sd->key_data.id);
            } else if(sci->sort_key_mask & Sort_KeyMask_LEN){
                g_print("len [%d] followed by [%d]",
                    sci->prev_sd->key_data.len, sd->key_data.len);
            } else if(sci->sort_key_mask & Sort_KeyMask_SEQ){
                g_print("seq [%s...] followed by [%s...]",
                    sci->prev_sd->key_data.seq_prefix,
                    sd->key_data.seq_prefix);
                }
            g_print("\n");
            exit(1);
            return TRUE;
            }
        Sort_Data_destroy(sci->prev_sd, sci->sort_key_mask);
        }
    sci->prev_sd = sd;
    return FALSE;
    }

static void fasta_sort_check_order(FastaDB *fdb,
                                   Sort_KeyMask sort_key_mask){
    Sort_Check_Info sci;
    sci.prev_sd = NULL;
    sci.comp_func = fasta_sort_get_compare_func(sort_key_mask);
    sci.sort_key_mask = sort_key_mask;
    FastaDB_traverse(fdb, FastaDB_Mask_ALL,
                     fasta_sort_check_traverse_func, &sci);
    if(sci.prev_sd)
        Sort_Data_destroy(sci.prev_sd, sci.sort_key_mask);
    return;
    }

/**/

typedef struct {
       GPtrArray *sort_data_list;
    Sort_KeyMask  sort_key_mask;
} Sort_Info;

static gboolean fasta_sort_traverse_func(FastaDB_Seq *fdbs,
                                         gpointer user_data){
    register Sort_Info *si = user_data;
    register Sort_Data *sd = Sort_Data_create(fdbs, si->sort_key_mask);
    g_ptr_array_add(si->sort_data_list, sd);
    return FALSE;
    }

static void fasta_sort_sort_data(FastaDB *fdb,
                                 Sort_KeyMask sort_key_mask){
    register gint i;
    register FastaDB_Seq *fdbs;
    register Sort_Data *sd;
    register FastaSort_CompFunc comp_func
             = fasta_sort_get_compare_func(sort_key_mask);
    Sort_Info si;
    si.sort_data_list = g_ptr_array_new();
    si.sort_key_mask = sort_key_mask;
    FastaDB_traverse(fdb, FastaDB_Mask_ALL,
                     fasta_sort_traverse_func, &si);
    qsort(si.sort_data_list->pdata, si.sort_data_list->len,
          sizeof(gpointer), comp_func);
    for(i = 0; i < si.sort_data_list->len; i++){
        sd = si.sort_data_list->pdata[i];
        fdbs = FastaDB_Key_get_seq(sd->key, FastaDB_Mask_ALL);
        FastaDB_Seq_print(fdbs, stdout, FastaDB_Mask_ID
                                       |FastaDB_Mask_DEF
                                       |FastaDB_Mask_SEQ);
        FastaDB_Seq_destroy(fdbs);
        }
    for(i = 0; i < si.sort_data_list->len; i++){
        sd = si.sort_data_list->pdata[i];
        Sort_Data_destroy(sd, sort_key_mask);
        }
    g_ptr_array_free(si.sort_data_list, TRUE);
    return;
    }

int Argument_main(Argument *arg){
    register FastaDB *fdb;
    register ArgumentSet *as
           = ArgumentSet_create("Sequence Input Options");
    gchar *query_path, *sort_key;
    gboolean check_order, reverse_order;
    register Sort_KeyMask sort_key_mask;
    ArgumentSet_add_option(as, 'f', "fasta", "path",
        "Fasta input file", NULL,
        Argument_parse_string, &query_path);
    ArgumentSet_add_option(as, 'c', "check", NULL,
        "Just check if file is already sorted (does not sort)", "FALSE",
        Argument_parse_boolean, &check_order);
    ArgumentSet_add_option(as, 'k', "key", "id | len | seq",
        "Sort key to apply to sequences", "id",
        Argument_parse_string, &sort_key);
    ArgumentSet_add_option(as, 'r', "reverse", NULL,
        "Reverse sort order", "FALSE",
        Argument_parse_boolean, &reverse_order);
    Argument_absorb_ArgumentSet(arg, as);
    Argument_process(arg, "fastasort",
        "A utility for sorting fasta sequence databases\n"
        "Guy St.C. Slater. guy@ebi.ac.uk. 2000-2003.\n", NULL);
    /**/
    fdb = FastaDB_open(query_path, NULL);
    sort_key_mask = SortKey_Mask_create(sort_key, reverse_order);
    if(check_order)
        fasta_sort_check_order(fdb, sort_key_mask);
    else
        fasta_sort_sort_data(fdb, sort_key_mask);
    FastaDB_close(fdb);
    return 0;
    }

