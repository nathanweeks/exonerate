/****************************************************************\
*                                                                *
*  fastafetch : fetch a sequence from a fasta format file        *
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
#include <strings.h>
#include <errno.h>

#include "argument.h"
#include "fastadb.h"

static gint fasta_index_seek_to_char(FILE *fp, gchar c){
    register gint ch;
    while((ch = getc(fp)) != EOF){
        if(ch == c)
            break;
        }
    return ch;
    }

static gint fasta_index_read_to_char(FILE *fp, gchar c, GString *str){
    register gint ch;
    g_string_truncate(str, 0);
    while((ch = getc(fp)) != EOF){
        if(ch == c)
            break;
        g_string_append_c(str, ch);
        }
    return ch;
    }

static CompoundFile_Pos fasta_index_lookup(FastaDB *fdb, FILE *index_fp,
                                           gchar *id){
    register long left, right, pos;
    CompoundFile_Pos_scan_type result = -1;
    register GString *curr = g_string_sized_new(100);
    register gint cond;
    left = 0;
    fseek(index_fp, 0, SEEK_END);
    right = ftell(index_fp);
    if(right == -1) {
        fprintf(stderr, "Could not seek to index end [%s]", strerror(errno));
        exit(EXIT_FAILURE);
    }
    while(left < right){ /* Binary search */
        pos = (left+right)>>1;
        fseek(index_fp, pos, SEEK_SET); /* Move between l&r */
        /* Move to next line (unless at beginning of index) */
        if(pos && (fasta_index_seek_to_char(index_fp, '\n') == EOF)){
            right = pos;
            continue;
            }
        if(fasta_index_read_to_char(index_fp, ' ', curr) == EOF){
            right = pos;
            continue;
            }
        cond = strcmp(id, curr->str);
        if(cond == 0){ /* Disco */
            fasta_index_read_to_char(index_fp, '\n', curr);
            /* Convert to position in result */
            CompoundFile_Pos_scan(curr->str, &result);
            break;
            }
        if(cond > 0){
            if((!pos) || (left == (pos-1)))
                break;
            left = pos - 1; /* Higher */
        } else {
            right = pos; /* Lower */
            }
        }
    g_string_free(curr, TRUE);
    return (CompoundFile_Pos)result;
    }

static FastaDB_Seq *fasta_index_fetch(FastaDB *fdb, FILE *index_fp,
                               FastaDB_Mask mask, gchar *id,
                               gboolean be_silent){
    register CompoundFile_Pos pos = fasta_index_lookup(fdb,
                                                       index_fp, id);
    register FastaDB_Seq *fdbs;
    if(pos == -1){
        if(be_silent)
            return NULL;
        fprintf(stderr, "Could not find identifier [%s] (missing -F ?)", id);
        exit(EXIT_FAILURE);
        }
    fdbs = FastaDB_fetch(fdb, mask, pos);
    g_assert(!strcmp(fdbs->seq->id, id));
    return fdbs;
    }

/**/

static gboolean read_next_line(FILE *fp, GString *s){
    register gint ch;
    g_string_truncate(s, 0);
    while((ch = getc(fp)) != EOF){
        if(ch == '\n'){
            return TRUE;
        } else {
            g_string_append_c(s, ch);
            }
        }
    return s->len?TRUE:FALSE;
    }

static GPtrArray *get_query_list(gchar *query, gboolean use_fosn){
    register GPtrArray *query_list = g_ptr_array_new();
    register FILE *fp = NULL;
    register GString *id;
    if(!strcasecmp(query, "stdin")){
        fp = stdin;
    } else {
        if(use_fosn)
             fp = fopen(query, "r");
        }
    if(fp){
        id = g_string_sized_new(64);
        while(read_next_line(fp, id))
            g_ptr_array_add(query_list, g_strdup(id->str));
        fclose(fp);
        g_string_free(id, TRUE);
    } else {
        g_ptr_array_add(query_list, g_strdup(query));
        }
    return query_list;
    }

static void fetch_sequences(gchar *fasta_path, gchar *index_path,
                            gboolean use_fosn, gboolean be_silent,
                            gchar *query){
    register GPtrArray *query_list = get_query_list(query, use_fosn);
    register gint i;
    register FastaDB *fdb = FastaDB_open(fasta_path, NULL);
    register FILE *index_fp = fopen(index_path, "r");
    register FastaDB_Seq *fdbs;
    register FastaDB_Mask mask = FastaDB_Mask_ID
                               | FastaDB_Mask_DEF
                               | FastaDB_Mask_SEQ;
    if(!index_fp) {
        fprintf(stderr,"Could not open fasta index [%s]", index_path);
        exit(EXIT_FAILURE);
    }
    for(i = 0; i < query_list->len; i++){
        fdbs = fasta_index_fetch(fdb, index_fp, mask,
                                 query_list->pdata[i], be_silent);
        if(fdbs){
            FastaDB_Seq_print(fdbs, stdout, mask);
            FastaDB_Seq_destroy(fdbs);
            }
        }
    for(i = 0; i < query_list->len; i++)
        g_free(query_list->pdata[i]);
    g_ptr_array_free(query_list, TRUE);
    FastaDB_close(fdb);
    fclose(index_fp);
    return;
    }

int Argument_main(Argument *arg){
    register ArgumentSet *as
           = ArgumentSet_create("Sequence Input Options");
    gchar *fasta_path, *index_path, *query;
    gboolean use_fosn, be_silent;
    ArgumentSet_add_option(as, 'f', "fasta", "path",
        "Fasta input file", NULL,
        Argument_parse_string, &fasta_path);
    ArgumentSet_add_option(as, 'i', "index", "path",
        "Index file path", NULL,
        Argument_parse_string, &index_path);
    ArgumentSet_add_option(as, 'F', "fosn", NULL,
        "Interpret query as FOSN (file of sequence names)", "FALSE",
        Argument_parse_boolean, &use_fosn);
    ArgumentSet_add_option(as, 'q', "query", "name",
        "Fetch query < id | fosn | stdin >", NULL,
        Argument_parse_string, &query);
    ArgumentSet_add_option(as, 's', "silent", NULL,
        "Silently skip ids which cannot be found", "FALSE",
        Argument_parse_boolean, &be_silent);
    Argument_absorb_ArgumentSet(arg, as);
    Argument_process(arg, "fastafetch",
        "Fetch a sequence from a fasta file\n"
        "Guy St.C. Slater. guy@ebi.ac.uk. 2000-2003.\n", NULL);
    fetch_sequences(fasta_path, index_path, use_fosn, be_silent, query);
    return 0;
    }

