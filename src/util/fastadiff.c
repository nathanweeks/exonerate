/****************************************************************\
*                                                                *
*  fastadiff : a utility to compare fasta sequence files         *
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

#include <string.h> /* For strcmp() */
#include <strings.h> /* For strcasecmp() */
#include <ctype.h>  /* For toupper() */

#include "argument.h"
#include "fastadb.h"

static gint fasta_diff_ambig_strcmp(const gchar *s1,
                                    const gchar *s2){
    register gchar c1, c2;
    while(*s1 && *s2){
        c1 = toupper(*s1);
        c2 = toupper(*s2);
        if((c1 != 'N') && (c1 != 'X') && (c2 != 'N') && (c2 != 'X')){
            if(*s1 != *s2)
                return *s1 - *s2;
            }
        s1++;
        s2++;
        }
    return *s1 - *s2;
    }

static gint fasta_diff_ambig_strcasecmp(const gchar *s1,
                                        const gchar *s2){
    register gchar c1, c2;
    while(*s1 && *s2){
        c1 = toupper(*s1);
        c2 = toupper(*s2);
        if((c1 != 'N') && (c1 != 'X') && (c2 != 'N') && (c2 != 'X')){
            if(c1 != c2)
                return c1 - c2;
            }
        s1++;
        s2++;
        }
    return *s1 - *s2;
    }

typedef gint (*fasta_diff_compare_func)(const gchar *s1,
                                        const gchar *s2);

static gint fasta_diff(gchar *first_path, gchar *second_path,
                   gboolean allow_ambiguity, gboolean ignore_case,
                   gboolean check_ids){
    register FastaDB *fdb_a, *fdb_b;
    register FastaDB_Seq *fdbs_a = NULL, *fdbs_b = NULL;
    register gint ret_val = 0;
    register FastaDB_Mask mask = FastaDB_Mask_ID|FastaDB_Mask_SEQ;
    register fasta_diff_compare_func comp = NULL;
    register gchar *seq_a, *seq_b;
    if(ignore_case){
        if(allow_ambiguity){
            comp = fasta_diff_ambig_strcasecmp;
        } else {
            comp = strcasecmp;
            }
    } else {
        if(allow_ambiguity){
            comp = fasta_diff_ambig_strcmp;
        } else {
            comp = strcmp;
            }
        }
    fdb_a = FastaDB_open(first_path, NULL);
    fdb_b = FastaDB_open(second_path, NULL);
    while((fdbs_a = FastaDB_next(fdb_a, mask))){
        fdbs_b = FastaDB_next(fdb_b, mask);
        if(!fdbs_b){
            g_print("fastadiff: %s: id %s absent\n",
                    second_path, fdbs_a->seq->id);
            ret_val = 1;
            break;
            }
        /* Check ids */
        if((check_ids)
        && (strcmp(fdbs_a->seq->id, fdbs_b->seq->id))){
            g_print("fastadiff: id mismatch: %s %s\n",
                    fdbs_a->seq->id, fdbs_b->seq->id);
            ret_val = 1;
            break;
            }
        /* Check length */
        if(fdbs_a->seq->len != fdbs_b->seq->len){
            g_print("fastadiff: length mismatch: %s(%d) %s(%d)\n",
                    fdbs_a->seq->id, fdbs_a->seq->len,
                    fdbs_b->seq->id, fdbs_b->seq->len);
            ret_val = 1;
            break;
            }
        /* Check seqs */
        seq_a = Sequence_get_str(fdbs_a->seq);
        seq_b = Sequence_get_str(fdbs_b->seq);
        if(comp(seq_a, seq_b)){
            g_print("fastadiff: sequence mismatch: %s %s\n",
                    fdbs_a->seq->id, fdbs_b->seq->id);
            ret_val = 1;
            break;
            }
        g_free(seq_a);
        g_free(seq_b);
        FastaDB_Seq_destroy(fdbs_a);
        FastaDB_Seq_destroy(fdbs_b);
        fdbs_a = fdbs_b = NULL;
        }
    if(fdbs_a)
        FastaDB_Seq_destroy(fdbs_a);
    if(fdbs_b)
        FastaDB_Seq_destroy(fdbs_b);
    if(ret_val == 0){
        fdbs_b = FastaDB_next(fdb_b, mask);
        if(fdbs_b){
            g_print("fastadiff: %s: id %s absent\n",
                    first_path, fdbs_b->seq->id);
            ret_val = 1;
            FastaDB_Seq_destroy(fdbs_b);
            }
        }
    FastaDB_close(fdb_a);
    FastaDB_close(fdb_b);
    return ret_val;
    }

int Argument_main(Argument *arg){
    register ArgumentSet *as
           = ArgumentSet_create("Sequence Input Options");
    gchar *first_path, *second_path;
    gboolean allow_ambiguity, ignore_case, check_ids;
    ArgumentSet_add_option(as, '1', "first", "path",
        "First fasta input file", NULL,
        Argument_parse_string, &first_path);
    ArgumentSet_add_option(as, '2', "second", "path",
        "Second fasta input file", NULL,
        Argument_parse_string, &second_path);
    ArgumentSet_add_option(as, 'i', "ignorecase", NULL,
        "Ignore sequence case", "FALSE",
        Argument_parse_boolean, &ignore_case);
    ArgumentSet_add_option(as, 'a', "ambiguity", NULL,
        "Allow sequence ambiguity with N or X", "FALSE",
        Argument_parse_boolean, &allow_ambiguity);
    ArgumentSet_add_option(as, 'c', "checkids", NULL,
        "Check ids for consistency", "TRUE",
        Argument_parse_boolean, &check_ids);
    Argument_absorb_ArgumentSet(arg, as);
    Argument_process(arg, "fastadiff",
        "Compare fasta format sequences (assumes same order)\n"
        "Guy St.C. Slater. guy@ebi.ac.uk. 2000-2003.\n", NULL);
    return fasta_diff(first_path, second_path,
               allow_ambiguity, ignore_case, check_ids);
    }

