/****************************************************************\
*                                                                *
*  The exonerate server                                          *
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
#include <glib.h>
#include <ctype.h> /* For isspace() */

#include "argument.h"
#include "socket.h"
#include "dataset.h"
#include "index.h"
#include "hspset.h"

typedef struct {
    Dataset *dataset;
      Index *index;
       gint  verbosity;
} Exonerate_Server;

static void Exonerate_Server_memory_usage(Exonerate_Server *exonerate_server){
    register guint64 dataset_memory, index_memory;
    dataset_memory = Dataset_memory_usage(exonerate_server->dataset),
    index_memory = exonerate_server->index
                 ? (Index_memory_usage(exonerate_server->index) - dataset_memory)
                 : 0;
    if(exonerate_server->verbosity > 0)
        g_message("Memory usage: dataset: %d Mb, index: %d Mb, Total %d Mb",
              (gint)(dataset_memory >> 20),
              (gint)(index_memory >> 20),
              (gint)((dataset_memory+index_memory) >> 20));
    return;
    }

static Exonerate_Server *Exonerate_Server_create(gchar *input_path,
                                                 gboolean preload,
                                                 gint verbosity){
    register Exonerate_Server *exonerate_server
     = g_new0(Exonerate_Server, 1);
    if(verbosity > 0)
        g_message("Starting server ...");
    if(Dataset_check_filetype(input_path)){
        exonerate_server->dataset = Dataset_read(input_path);
    } else if(Index_check_filetype(input_path)){
        exonerate_server->index = Index_open(input_path);
        exonerate_server->dataset
            = Dataset_share(exonerate_server->index->dataset);
        if(preload)
            Index_preload_index(exonerate_server->index);
    } else {
        g_error("Unknown filetype for input file [%s]", input_path);
        }
    if(preload)
        Dataset_preload_seqs(exonerate_server->dataset);
    exonerate_server->verbosity = verbosity;
    if(verbosity >= 1){
        Dataset_info(exonerate_server->dataset);
        Index_info(exonerate_server->index);
        }
    return exonerate_server;
    }

static void Exonerate_Server_destroy(Exonerate_Server *exonerate_server){
    Dataset_destroy(exonerate_server->dataset);
    if(exonerate_server->index)
        Index_destroy(exonerate_server->index);
    Dataset_destroy(exonerate_server->dataset);
    g_free(exonerate_server);
    return;
    }

/**/

typedef struct {
           Alphabet *query_alphabet;
           Alphabet_Type query_type;
    Sequence_Strand  query_strand;
           gboolean  query_is_masked;
           gboolean  revcomp_query;
           gboolean  revcomp_target;
           Sequence *query;
          HSP_Param *hsp_param;
          /**/
               gint  seed_repeat;
               gint  dna_hsp_threshold;
               gint  protein_hsp_threshold;
               gint  codon_hsp_threshold;
               gint  dna_word_limit;
               gint  protein_word_limit;
               gint  codon_word_limit;
               gint  dna_hsp_dropoff;
               gint  protein_hsp_dropoff;
               gint  codon_hsp_dropoff;
               gint  geneseed_threshold;
               gint  geneseed_repeat;
               gint  max_query_span;
               gint  max_target_span;
} Exonerate_Server_Connection;

static Exonerate_Server_Connection *Exonerate_Server_Connection_create(void){
    register Exonerate_Server_Connection *connection
     = g_new(Exonerate_Server_Connection, 1);
    register HSPset_ArgumentSet *has = HSPset_ArgumentSet_create(NULL);
    connection->query = NULL;
    connection->hsp_param = NULL;
    connection->query_alphabet = NULL;
    connection->query_type = Alphabet_Type_UNKNOWN;
    connection->query_strand = Sequence_Strand_UNKNOWN;
    connection->query_is_masked = FALSE;
    connection->revcomp_query = FALSE;
    connection->revcomp_target = FALSE;
    /**/
    connection->seed_repeat = has->seed_repeat;
    connection->dna_hsp_threshold = has->dna_hsp_threshold;
    connection->protein_hsp_threshold = has->protein_hsp_threshold;
    connection->codon_hsp_threshold = has->codon_hsp_threshold;
    connection->dna_word_limit = has->dna_word_limit;
    connection->protein_word_limit = has->protein_word_limit;
    connection->codon_word_limit = has->codon_word_limit;
    connection->dna_hsp_dropoff = has->dna_hsp_dropoff;
    connection->protein_hsp_dropoff = has->protein_hsp_dropoff;
    connection->codon_hsp_dropoff = has->codon_hsp_dropoff;
    connection->geneseed_threshold = has->geneseed_threshold;
    connection->geneseed_repeat = has->geneseed_repeat;
    connection->max_query_span = 0;
    connection->max_target_span = 0;
    return connection;
    }

static void Exonerate_Server_Connection_destroy(
            Exonerate_Server_Connection *connection){
    if(connection->query_alphabet)
        Alphabet_destroy(connection->query_alphabet);
    if(connection->hsp_param)
        HSP_Param_destroy(connection->hsp_param);
    if(connection->query)
        Sequence_destroy(connection->query);
    g_free(connection);
    return;
    }

/**/

static gpointer Exonerate_Server_Connection_open(gpointer user_data){
    return Exonerate_Server_Connection_create();
    }

static void Exonerate_Server_Connection_close(gpointer connection_data,
                                              gpointer user_data){
    register Exonerate_Server_Connection *server_connection
        = connection_data;
    Exonerate_Server_Connection_destroy(server_connection);
    return;
    }

static void Exonerate_Server_Connection_revcomp_query(
            Exonerate_Server_Connection *connection){
    register Sequence *rc_seq;
    g_assert(connection->query);
    rc_seq = Sequence_revcomp(connection->query);
    Sequence_destroy(connection->query);
    connection->query = rc_seq;
    connection->revcomp_query = connection->revcomp_query?FALSE:TRUE;
    return;
    }

static void Exonerate_Server_Connection_revcomp_target(
            Exonerate_Server_Connection *connection){
    connection->revcomp_target = connection->revcomp_target?FALSE:TRUE;
    return;
    }

static GPtrArray *Exonerate_Server_get_word_list(gchar *msg){
    register GPtrArray *word_list = g_ptr_array_new();
    register gchar *prev, *ptr;
    for(ptr = msg; isspace(*ptr); ptr++); /* skip start */
    prev = ptr;
    while(*ptr){
        if(isspace(*ptr)){
            *ptr = '\0';
            do {
                ptr++;
            } while(isspace(*ptr));
            if(!*ptr)
                break;
            g_ptr_array_add(word_list, prev); /* add a word */
            prev = ptr;
            }
        ptr++;
        }
    if(prev != ptr)
        g_ptr_array_add(word_list, prev); /* add final word */
    return word_list;
    }

static gchar *Exonerate_Server_help(void){
    return g_strdup_printf(
        "exonerate-server commands:\n"
        "    help    : print this message\n"
        "    version : show version information\n"
        "    exit    : disconnect from server\n"
        "    dbinfo  : show database info\n"
        "            : <type> <masked> <num_seqs> <max_seq_len> <total_seq_len>\n"
        "\n"
        "    lookup <eid> : get internal from external identifier\n"
        "    get info <iid> : get sequence info \n"
        "                   : <len> <checksum> <id> [<def>]\n"
        "    get seq <iid> : get sequence\n"
        "    get subseq <iid> <start> <len> : get subsequence\n"
        "\n"
        "    set query <seq> : set query sequence\n"
        "    get hsps : get hsps against current query\n"
        "             : <target_id> { <query_pos> <target_pos> <length> } \n"
        "\n"
        "    revcomp <query | target>\n"
        "    set param <name> <value>\n"
        "\n"
        "\n"
        "    valid parameters:\n"
        "        querytype\n"
        "        seedrepeat\n"
        "\n"
        "        dnahspthreshold\n"
        "        proteinhspthreshold\n"
        "        codonhspthreshold\n"
        "\n"
        "        dnawordlimit\n"
        "        proteinwordlimit\n"
        "        codonwordlimit\n"
        "\n"
        "        geneseedthreshold\n"
        "        geneseedrepeat\n"
        "        maxqueryspan\n"
        "        maxtargetspan\n"
        "--\n");
    }

static gchar *Exonerate_Server_get_info(Dataset *dataset, gint num){
    register Dataset_Sequence *ds;
    register Sequence *seq;
    register gchar *reply;
    if((num >= 0) && (num < dataset->seq_list->len)){
        ds = dataset->seq_list->pdata[num];
        seq = Dataset_get_sequence(dataset, num);
        reply = g_strdup_printf("seqinfo: %d %d %s%s%s\n",
                            (gint)ds->key->length,
                            (gint)ds->gcg_checksum,
                            ds->id,
                            seq->def?" ":"",
                            seq->def?seq->def:"");
        Sequence_destroy(seq);
    } else {
        reply = g_strdup_printf("error: sequence num out of range [%d]\n", num);
        }
    return reply;
    }

static gchar *Exonerate_Server_get_seq(Dataset *dataset, gint num){
    register Dataset_Sequence *ds;
    register Sequence *seq;
    register gchar *str, *reply;
    if((num >= 0) && (num < dataset->seq_list->len)){
        ds = dataset->seq_list->pdata[num];
        seq = Dataset_get_sequence(dataset, num);
        str = Sequence_get_str(seq);
        Sequence_destroy(seq);
        reply = g_strdup_printf("seq: %s\n", str);
        g_free(str);
    } else {
        reply = g_strdup_printf("error: sequence num out of range [%d]\n", num);
        }
    return reply;
    }

static gchar *Exonerate_Server_get_subseq(Dataset *dataset, gint num,
                                          gint start, gint len){
    register gchar *reply, *str;
    register Dataset_Sequence *ds;
    register Sequence *seq, *subseq;
    if((num >= 0) && (num < dataset->seq_list->len)){
        ds = dataset->seq_list->pdata[num];
        if(len <= 0){
            reply = g_strdup_printf("error: subseq len (%d) must be >= 0\n", len);
        } else if((start >= 0) && ((start+len) <= ds->key->length)){
            seq = Dataset_get_sequence(dataset, num);
            subseq = Sequence_subseq(seq, start, len);
            Sequence_destroy(seq);
            str = Sequence_get_str(subseq);
            Sequence_destroy(subseq);
            reply = g_strdup_printf("subseq: %s\n", str);
            g_free(str);
        } else {
            reply = g_strdup_printf("error: subsequence beyond seq len [%d]\n",
                    ds->key->length);
            }
    } else {
        reply = g_strdup_printf("error: sequence num out of range [%d]\n", num);
        }
    return reply;
    }

static gchar *Exonerate_Server_get_hsps(Exonerate_Server *exonerate_server,
                                        Exonerate_Server_Connection *connection){
    register gint i, j, hsp_total = 0;
    register GPtrArray *index_hsp_set_list;
    register Index_HSPset *index_hsp_set;
    register gchar *reply;
    register GString *str;
    register HSP *hsp;
    g_assert(connection->hsp_param);
    g_assert(connection->query);
    if(connection->revcomp_target
    && (connection->hsp_param->match->type != Match_Type_PROTEIN2DNA))
        return g_strdup_printf(
                "error: revcomp target only available for protein2dna matches");
    if(connection->geneseed_threshold > 0){
        if(connection->geneseed_threshold < connection->hsp_param->threshold)
            return g_strdup_printf(
                    "error: geneseed threshold must be >= hsp threshold");
        index_hsp_set_list = Index_get_HSPsets_geneseed(exonerate_server->index,
                                               connection->hsp_param,
                                               connection->query,
                                               connection->revcomp_target,
                                               connection->geneseed_threshold,
                                               connection->geneseed_repeat,
                                               connection->max_query_span,
                                               connection->max_target_span);
    } else {
        index_hsp_set_list = Index_get_HSPsets(exonerate_server->index,
                                               connection->hsp_param,
                                               connection->query,
                                               connection->revcomp_target);
        }
    if(index_hsp_set_list){
        str = g_string_sized_new(1024);
        for(i = 0; i < index_hsp_set_list->len; i++){
            index_hsp_set = index_hsp_set_list->pdata[i];
            g_assert(index_hsp_set->hsp_set->is_finalised);
            g_string_sprintfa(str, "hspset: %d", index_hsp_set->target_id);
            /**/
            hsp_total += index_hsp_set->hsp_set->hsp_list->len;
            for(j = 0; j < index_hsp_set->hsp_set->hsp_list->len; j++){
                hsp = index_hsp_set->hsp_set->hsp_list->pdata[j];
                g_string_sprintfa(str, " %d %d %d",
                        hsp->query_start, hsp->target_start, hsp->length);
                }
            g_string_sprintfa(str, "\n");
            Index_HSPset_destroy(index_hsp_set);
            }
        if(exonerate_server->verbosity > 1)
            g_message("served [%d] HSPsets containing [%d] hsps",
                      index_hsp_set_list->len, hsp_total);
        g_ptr_array_free(index_hsp_set_list, TRUE);
        reply = str->str;
        g_assert(reply);
        g_string_free(str, FALSE);
    } else {
        reply = g_strdup_printf("hspset: empty\n");
        }
    return reply;
    }

static Sequence *Exonerate_Server_get_query(Index *index,
                 Exonerate_Server_Connection *connection, gchar *query){
    register Alphabet_Type alphabet_type;
    register Match_Type match_type;
    register Match *match;
    if(!connection->query_alphabet){
        if (connection->query_type != Alphabet_Type_UNKNOWN)
            alphabet_type = connection->query_type; /* client-specified type */
        else
            alphabet_type = Alphabet_Type_guess(query);

        if(alphabet_type == Alphabet_Type_DNA){
            connection->query_strand = Sequence_Strand_FORWARD;
        } else {
            g_assert(alphabet_type == Alphabet_Type_PROTEIN);
            connection->query_strand = Sequence_Strand_UNKNOWN;
            if((index->dataset->alphabet->type == Alphabet_Type_DNA)
            && (!(index->header->type & 1))){
                g_message("Cannot use protein query with untranslated DNA index");
                return NULL;
                }
            }
        connection->query_alphabet = Alphabet_create(alphabet_type,
                                                     connection->query_is_masked);
        match_type = Match_Type_find(alphabet_type,
                                     index->dataset->alphabet->type, FALSE);
        /* FIXME: use Match_Type_find with translate_both for codon alignments */
        match = Match_find(match_type);
        g_assert(match);
        connection->hsp_param = HSP_Param_create(match, FALSE);
        connection->hsp_param->seed_repeat = connection->seed_repeat;
        /**/
        HSP_Param_set_dna_hsp_threshold(connection->hsp_param,
                                        connection->dna_hsp_threshold);
        HSP_Param_set_protein_hsp_threshold(connection->hsp_param,
                                            connection->protein_hsp_threshold);
        HSP_Param_set_codon_hsp_threshold(connection->hsp_param,
                                          connection->codon_hsp_threshold);
        /**/
        HSP_Param_set_dna_word_limit(connection->hsp_param,
                                     connection->dna_word_limit);
        HSP_Param_set_protein_word_limit(connection->hsp_param,
                                         connection->protein_word_limit);
        HSP_Param_set_codon_word_limit(connection->hsp_param,
                                       connection->codon_word_limit);
        /**/
        HSP_Param_set_dna_hsp_dropoff(connection->hsp_param,
                                      connection->dna_hsp_dropoff);
        HSP_Param_set_protein_hsp_dropoff(connection->hsp_param,
                                          connection->protein_hsp_dropoff);
        HSP_Param_set_codon_hsp_dropoff(connection->hsp_param,
                                        connection->codon_hsp_dropoff);
        /**/
        g_assert(connection->hsp_param);
        HSP_Param_set_wordlen(connection->hsp_param, index->header->word_length);
        }
    g_assert(connection->query_alphabet);
    return Sequence_create("query", NULL, query, 0,
                               connection->query_strand,
                               connection->query_alphabet);
    }

static gchar *Exonerate_Server_set_param_querytype(
              Exonerate_Server_Connection *connection, GPtrArray *word_list){
    const char *querytype = word_list->pdata[3];
    if (!strcmp(querytype, "dna"))
        connection->query_type = Alphabet_Type_DNA;
    else if (!strcmp(querytype, "protein"))
        connection->query_type = Alphabet_Type_PROTEIN;
    else
        return g_strdup_printf(
               "error: querytype must be \"dna\" or \"protein\"\n");
    return g_strdup_printf("ok: set\n");
    }
static gchar *Exonerate_Server_set_param_seedrepeat(
              Exonerate_Server_Connection *connection, GPtrArray *word_list){
    register gint seed_repeat = atoi(word_list->pdata[3]);
    if(seed_repeat < 1)
        return g_strdup_printf("error: seedrepeat must be > 0\n");
    connection->seed_repeat = seed_repeat;
    if(connection->hsp_param)
        connection->hsp_param->seed_repeat = seed_repeat;
    return g_strdup_printf("ok: set\n");
    }

/**/

static gchar *Exonerate_Server_set_param_dnahspthreshold(
              Exonerate_Server_Connection *connection, GPtrArray *word_list){
    register gint dnahspthreshold = atoi(word_list->pdata[3]);
    if(dnahspthreshold < 1)
        return g_strdup_printf("error: dnahspthreshold must be > 0\n");
    connection->dna_hsp_threshold = dnahspthreshold;
    if(connection->hsp_param)
        HSP_Param_set_dna_hsp_threshold(connection->hsp_param, dnahspthreshold);
    return g_strdup_printf("ok: set\n");
    }

static gchar *Exonerate_Server_set_param_proteinhspthreshold(
              Exonerate_Server_Connection *connection, GPtrArray *word_list){
    register gint proteinhspthreshold = atoi(word_list->pdata[3]);
    if(proteinhspthreshold < 1)
        return g_strdup_printf("error: proteinhspthreshold must be > 0\n");
    connection->protein_hsp_threshold = proteinhspthreshold;
    if(connection->hsp_param)
        HSP_Param_set_protein_hsp_threshold(connection->hsp_param,
                                            proteinhspthreshold);
    return g_strdup_printf("ok: set\n");
    }

static gchar *Exonerate_Server_set_param_codonhspthreshold(
              Exonerate_Server_Connection *connection, GPtrArray *word_list){
    register gint codonhspthreshold = atoi(word_list->pdata[3]);
    if(codonhspthreshold < 1)
        return g_strdup_printf("error: codonhspthreshold must be > 0\n");
    connection->codon_hsp_threshold = codonhspthreshold;
    if(connection->hsp_param)
        HSP_Param_set_codon_hsp_threshold(connection->hsp_param,
                                          codonhspthreshold);
    return g_strdup_printf("ok: set\n");
    }

/**/

static gchar *Exonerate_Server_set_param_dnawordlimit(
              Exonerate_Server_Connection *connection, GPtrArray *word_list){
    register gint dnawordlimit = atoi(word_list->pdata[3]);
    if(dnawordlimit < 0)
        return g_strdup_printf("error: dnawordlimit must be >= 0\n");
    connection->dna_word_limit = dnawordlimit;
    if(connection->hsp_param)
        HSP_Param_set_dna_word_limit(connection->hsp_param, dnawordlimit);
    return g_strdup_printf("ok: set\n");
    }

static gchar *Exonerate_Server_set_param_proteinwordlimit(
              Exonerate_Server_Connection *connection, GPtrArray *word_list){
    register gint proteinwordlimit = atoi(word_list->pdata[3]);
    if(proteinwordlimit < 0)
        return g_strdup_printf("error: proteinwordlimit must be >= 0\n");
    connection->protein_word_limit = proteinwordlimit;
    if(connection->hsp_param)
        HSP_Param_set_protein_word_limit(connection->hsp_param,
                                         proteinwordlimit);
    return g_strdup_printf("ok: set\n");
    }

static gchar *Exonerate_Server_set_param_codonwordlimit(
              Exonerate_Server_Connection *connection, GPtrArray *word_list){
    register gint codonwordlimit = atoi(word_list->pdata[3]);
    if(codonwordlimit < 0)
        return g_strdup_printf("error: codonwordlimit must be >= 0\n");
    connection->codon_word_limit = codonwordlimit;
    if(connection->hsp_param)
        HSP_Param_set_codon_word_limit(connection->hsp_param,
                                       codonwordlimit);
    return g_strdup_printf("ok: set\n");
    }

/**/

static gchar *Exonerate_Server_set_param_dnahspdropoff(
              Exonerate_Server_Connection *connection, GPtrArray *word_list){
    register gint dnahspdropoff = atoi(word_list->pdata[3]);
    if(dnahspdropoff < 0)
        return g_strdup_printf("error: dnahspdropoff must be >= 0\n");
    connection->dna_hsp_dropoff = dnahspdropoff;
    if(connection->hsp_param)
        HSP_Param_set_dna_hsp_dropoff(connection->hsp_param, dnahspdropoff);
    return g_strdup_printf("ok: set\n");
    }

static gchar *Exonerate_Server_set_param_proteinhspdropoff(
              Exonerate_Server_Connection *connection, GPtrArray *word_list){
    register gint proteinhspdropoff = atoi(word_list->pdata[3]);
    if(proteinhspdropoff < 0)
        return g_strdup_printf("error: proteinhspdropoff must be >= 0\n");
    connection->protein_hsp_dropoff = proteinhspdropoff;
    if(connection->hsp_param)
        HSP_Param_set_protein_hsp_dropoff(connection->hsp_param,
                                          proteinhspdropoff);
    return g_strdup_printf("ok: set\n");
    }

static gchar *Exonerate_Server_set_param_codonhspdropoff(
              Exonerate_Server_Connection *connection, GPtrArray *word_list){
    register gint codonhspdropoff = atoi(word_list->pdata[3]);
    if(codonhspdropoff < 0)
        return g_strdup_printf("error: codonhspdropoff must be >= 0\n");
    connection->codon_hsp_dropoff = codonhspdropoff;
    if(connection->hsp_param)
        HSP_Param_set_codon_hsp_dropoff(connection->hsp_param,
                                        codonhspdropoff);
    return g_strdup_printf("ok: set\n");
    }

/**/

static gchar *Exonerate_Server_set_param_geneseedthreshold(
              Exonerate_Server_Connection *connection, GPtrArray *word_list){
    register gint geneseed_threshold = atoi(word_list->pdata[3]);
    if(geneseed_threshold < 0)
        return g_strdup_printf("error: geneseed_threshold must be >= 0\n");
    connection->geneseed_threshold = geneseed_threshold;
    return g_strdup_printf("ok: set\n");
    }

static gchar *Exonerate_Server_set_param_geneseedrepeat(
              Exonerate_Server_Connection *connection, GPtrArray *word_list){
    register gint geneseed_repeat = atoi(word_list->pdata[3]);
    if(geneseed_repeat <= 1)
        return g_strdup_printf("error: geneseed_repeat must be > 1\n");
    connection->geneseed_repeat = geneseed_repeat;
    return g_strdup_printf("ok: set\n");
    }

static gchar *Exonerate_Server_set_param_max_query_span(
              Exonerate_Server_Connection *connection, GPtrArray *word_list){
    register gint max_query_span= atoi(word_list->pdata[3]);
    if(max_query_span < 0)
        return g_strdup_printf("error: max_query_span must be >= 0\n");
    connection->max_query_span = max_query_span;
    return g_strdup_printf("ok: set\n");
    }

static gchar *Exonerate_Server_set_param_max_target_span(
              Exonerate_Server_Connection *connection, GPtrArray *word_list){
    register gint max_target_span= atoi(word_list->pdata[3]);
    if(max_target_span < 0)
        return g_strdup_printf("error: max_target_span must be >= 0\n");
    connection->max_target_span = max_target_span;
    return g_strdup_printf("ok: set\n");
    }

static gchar *Exonerate_Server_set_param(Exonerate_Server_Connection *connection,
                                         GPtrArray *word_list){
    register gchar *reply = NULL;
    register gchar *name = word_list->pdata[2];
    if (!strcmp(name, "querytype")){
        reply = Exonerate_Server_set_param_querytype(connection, word_list);
    } else if(!strcmp(name, "seedrepeat")){
        reply = Exonerate_Server_set_param_seedrepeat(connection, word_list);
    } else if(!strcmp(name, "dnahspthreshold")){
        reply = Exonerate_Server_set_param_dnahspthreshold(connection, word_list);
    } else if(!strcmp(name, "proteinhspthreshold")){
        reply = Exonerate_Server_set_param_proteinhspthreshold(connection,
                                                               word_list);
    } else if(!strcmp(name, "codonhspthreshold")){
        reply = Exonerate_Server_set_param_codonhspthreshold(connection,
                                                             word_list);
    } else if(!strcmp(name, "dnawordlimit")){
        reply = Exonerate_Server_set_param_dnawordlimit(connection,
                                                          word_list);
    } else if(!strcmp(name, "proteinwordlimit")){
        reply = Exonerate_Server_set_param_proteinwordlimit(connection,
                                                              word_list);
    } else if(!strcmp(name, "codonwordlimit")){
        reply = Exonerate_Server_set_param_codonwordlimit(connection,
                                                            word_list);
    } else if(!strcmp(name, "dnahspdropoff")){
        reply = Exonerate_Server_set_param_dnahspdropoff(connection, word_list);
    } else if(!strcmp(name, "proteinhspdropoff")){
        reply = Exonerate_Server_set_param_proteinhspdropoff(connection,
                                                             word_list);
    } else if(!strcmp(name, "codonhspdropoff")){
        reply = Exonerate_Server_set_param_codonhspdropoff(connection,
                                                           word_list);
    } else if(!strcmp(name, "geneseedthreshold")){
        reply = Exonerate_Server_set_param_geneseedthreshold(connection,
                                                             word_list);
    } else if(!strcmp(name, "geneseedrepeat")){
        reply = Exonerate_Server_set_param_geneseedrepeat(connection,
                                                          word_list);
    } else if(!strcmp(name, "maxqueryspan")){
        reply = Exonerate_Server_set_param_max_query_span(connection,
                                                          word_list);
    } else if(!strcmp(name, "maxtargetspan")){
        reply = Exonerate_Server_set_param_max_target_span(connection,
                                                           word_list);
    } else {
        reply = g_strdup_printf(
                "warning: set param %s ignored by server\n", name);
        }
    g_assert(reply);
    return reply;
    }

// allow compilation with "-Werror"
#pragma GCC diagnostic warning "-Wformat"
static gboolean Exonerate_Server_process(gchar *msg, gchar **reply,
                                         gpointer connection_data,
                                         gpointer user_data){
    register gint msg_len = strlen(msg);
    register gboolean keep_connection = TRUE;
    register Exonerate_Server *server = user_data;
    register GPtrArray *word_list;
    register gchar *word, *id, *item, *query;
    register gint start, len, num;
    register Exonerate_Server_Connection *connection = connection_data;
    g_assert(msg);
    if(server->verbosity >= 3)
        g_print("Message: server received command [%s]\n", msg);
    if(!msg[msg_len-1] == '\n')
        msg[--msg_len] = '\0';
    word_list = Exonerate_Server_get_word_list(msg);
    (*reply) = NULL;
    if(word_list->len){
        word = word_list->pdata[0];
        if(!strcmp(word, "help")){
            (*reply) = Exonerate_Server_help();
        } else if(!strcmp(word, "version")){
            (*reply) = g_strdup_printf("version: exonerate-server %s\n", VERSION);
        } else if(!strcmp(word, "exit")){
            (*reply) = NULL;
            keep_connection = FALSE;
        } else if(!strcmp(word, "dbinfo")){
            (*reply) = g_strdup_printf("dbinfo: %s %s"
                   " %" CUSTOM_GUINT64_FORMAT
                   " %" CUSTOM_GUINT64_FORMAT
                   " %" CUSTOM_GUINT64_FORMAT "\n",
                  (server->dataset->header->type & 1)?"dna":"protein",
                  (server->dataset->header->type & (1<<1))
                      ?"softmasked":"unmasked",
                  server->dataset->header->number_of_seqs,
                  server->dataset->header->max_seq_len,
                  server->dataset->header->total_seq_len);
        } else if(!strcmp(word, "lookup")){
            if(word_list->len == 2){
                id = word_list->pdata[1];
                num = Dataset_lookup_id(server->dataset, id);
                if(num == -1)
                    (*reply) = g_strdup_printf("error: id not found\n");
                else
                    (*reply) = g_strdup_printf("lookup: %d\n", num);
            } else {
                (*reply) = g_strdup_printf("error: usage: lookup <id>\n");
                }
        } else if(!strcmp(word, "get")){
            if(word_list->len >= 2){
                item = word_list->pdata[1];
                if(!strcmp(item, "info")){
                    if(word_list->len == 3){
                        num = atoi(word_list->pdata[2]);
                        (*reply) = Exonerate_Server_get_info(server->dataset, num);
                    } else {
                        (*reply) = g_strdup_printf(
                                        "error: usage: get info <pos>\n");
                        }
                } else if(!strcmp(item, "seq")){
                    if(word_list->len == 3){
                        num = atoi(word_list->pdata[2]);
                        (*reply) = Exonerate_Server_get_seq(server->dataset, num);
                    } else {
                        (*reply) = g_strdup_printf(
                                      "error: usage: get seq <pos>\n");
                        }
                } else if(!strcmp(item, "subseq")){
                    if(word_list->len == 5){
                        num = atoi(word_list->pdata[2]);
                        start = atoi(word_list->pdata[3]);
                        len = atoi(word_list->pdata[4]);
                        (*reply) = Exonerate_Server_get_subseq(server->dataset,
                                                               num, start, len);
                    } else {
                        (*reply) = g_strdup_printf(
                            "error: usage: get subseq <pos> <start> <length>\n");
                        }
                } else if(!strcmp(item, "hsps")){
                    if(connection->query){
                        if(server->index){
                            (*reply)
                                = Exonerate_Server_get_hsps(server,
                                                            connection);
                        } else {
                            (*reply)
                                = g_strdup_printf("error: no index for hsps\n");
                            }
                    } else {
                        (*reply) = g_strdup_printf("error: query not set\n");
                        }
                } else {
                    (*reply) = g_strdup_printf("error: Unknown get command\n");
                    }
            } else {
                (*reply) = g_strdup_printf("error: get what ?\n");
                }
        } else if(!strcmp(word, "set")){
            if(word_list->len >= 2){
                item = word_list->pdata[1];
                if(!strcmp(item, "query")){
                    if(word_list->len == 3){
                        query = word_list->pdata[2];
                        if(connection->query)
                            Sequence_destroy(connection->query);
                        connection->query = Exonerate_Server_get_query(
                                                         server->index,
                                                         connection, query);
                        connection->revcomp_query = FALSE;
                        if(connection->query)
                            (*reply) = g_strdup_printf("ok: %d %d\n",
                                    connection->query->len,
                                    Sequence_checksum(connection->query));
                        else
                            (*reply) = g_strdup_printf("error: bad query\n");
                    } else {
                        (*reply) = g_strdup_printf(
                                "error: usage: set query <seq>\n");
                        }
                } else if(!strcmp(item, "param")){
                    if(word_list->len < 3){
                        (*reply) = g_strdup_printf(
                                "error: usage: set param <name> <value>\n");
                    } else {
                        (*reply) = Exonerate_Server_set_param(connection,
                                                              word_list);
                        }
                } else {
                    (*reply) = g_strdup_printf("error: Unknown set command\n");
                    }
            } else {
                (*reply) = g_strdup_printf("error: set what ?\n");
                }
        } else if(!strcmp(word, "revcomp")){
            if(word_list->len == 2){
                item = word_list->pdata[1];
                if(!strcmp(item, "query")){
                    if(!connection->query)
                        (*reply) = g_strdup_printf("error: query not set\n");
                    if(connection->query_alphabet->type == Alphabet_Type_DNA){
                        Exonerate_Server_Connection_revcomp_query(connection);
                        (*reply) = g_strdup_printf(
                                  "ok: query strand %s\n",
                                  connection->revcomp_query?"revcomp":"forward");
                    } else {
                        (*reply) = g_strdup_printf(
                                "error: cannot revcomp non-DNA query\n");
                        }
                } else if(!strcmp(item, "target")){
                    Exonerate_Server_Connection_revcomp_target(connection);
                    (*reply) = g_strdup_printf(
                                  "ok: target strand %s\n",
                                  connection->revcomp_target?"revcomp":"forward");
                } else {
                    (*reply) = g_strdup_printf("error: Unknown revcomp command\n");
                    }
            } else {
                (*reply) = g_strdup_printf("error: revcomp what ?\n");
                }
        } else {
            (*reply) = g_strdup_printf("error: Unknown command: [%s]\n", msg);
            }
        }
    g_ptr_array_free(word_list, TRUE);
    if((server->verbosity >= 3) && (*reply)){
        if((server->verbosity == 3) && (strlen(*reply) >= 80))
            g_print("Message: server returned reply [%.*s<truncated>]\n",
                    80, (*reply));
        else
            g_print("Message: server returned reply [%s]\n", (*reply));
        }
    return keep_connection;
    }

static void run_server(gint port, gchar *input_path, gboolean preload,
                       gint max_connections, gint verbosity){
    register Exonerate_Server *exonerate_server
           = Exonerate_Server_create(input_path, preload, verbosity);
    register SocketServer *ss = SocketServer_create(port, max_connections,
                       Exonerate_Server_process,
                       Exonerate_Server_Connection_open,
                       Exonerate_Server_Connection_close,
                       exonerate_server);
    Exonerate_Server_memory_usage(exonerate_server);
    if(verbosity > 0)
        g_message("listening on port [%d] ...", port);
    while(SocketServer_listen(ss));
    SocketServer_destroy(ss);
    Exonerate_Server_destroy(exonerate_server);
    return;
    }

int Argument_main(Argument *arg){
    gint port, max_connections, verbosity;
    gchar *input_path;
    gboolean preload;
    register ArgumentSet *as = ArgumentSet_create("Exonerate Server options");
    ArgumentSet_add_option(as, '\0', "port", "port",
            "Port number to run server on", "12886",
            Argument_parse_int, &port);
    ArgumentSet_add_option(as, '\0', "input", "path",
            "Path to input file (.esd or .esi)", NULL,
            Argument_parse_string, &input_path);
    ArgumentSet_add_option(as, '\0', "preload", NULL,
            "Preload index and sequence data", "TRUE",
            Argument_parse_boolean, &preload);
    ArgumentSet_add_option(as, '\0', "maxconnections", "threads",
            "Maximum concurrent server connections", "4",
            Argument_parse_int, &max_connections);
    ArgumentSet_add_option(as, 'V', "verbosity", "level",
            "Set server verbosity level", "1",
            Argument_parse_int, &verbosity);
    Argument_absorb_ArgumentSet(arg, as);
    /**/
    Match_ArgumentSet_create(arg);
    HSPset_ArgumentSet_create(arg);
    /**/
    Argument_process(arg, "exonerate-server", "Exonerate Server.\n",
                     "Guy St.C. Slater.  guy@ebi.ac.uk June 2006\n");
    run_server(port, input_path, preload, max_connections, verbosity);
    g_message("-- server exiting");
    return 0;
    }

/**/

