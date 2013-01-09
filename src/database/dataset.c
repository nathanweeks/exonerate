/****************************************************************\
*                                                                *
*  Library for manipulation of exonerate dataset files           *
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

#include <stdlib.h> /* For qsort() */
#include <string.h> /* For strcmp() */

#include "dataset.h"
#include "fastadb.h"
#include "bitarray.h"

/**/

#define DATASET_HEADER_MAGIC (('e' << 16)|('s' << 8)|('d'))
#define DATASET_HEADER_VERSION 3

static Dataset_Header *Dataset_Header_create(Alphabet_Type alphabet_type,
                                             gboolean softmask_input){
    register Dataset_Header *header = g_new0(Dataset_Header, 1);
    header->magic = DATASET_HEADER_MAGIC;
    header->version = DATASET_HEADER_VERSION;
    header->type = ((alphabet_type == Alphabet_Type_DNA)?1:0)
                 | ((softmask_input?1:0) << 1);
    /* other fields are filled during seq parsing or Dataset_finalise */
    return header;
    }

static void Dataset_Header_destroy(Dataset_Header *header){
    g_free(header);
    return;
    }

#if 0
static void Dataset_Header_info(Dataset_Header *header){
    g_message("Dataset_Header:\n"
              "     magic             [%lld]\n"
              "     version           [%lld]\n"
              "     type              [%lld]\n"
              "     line_length       [%lld]\n"
              /**/
              "     number_of_dbs     [%lld]\n"
              "     max_db_len        [%lld]\n"
              "     total_db_len      [%lld]\n"
              /**/
              "     number_of_seqs    [%lld]\n"
              "     max_seq_len       [%lld]\n"
              "     total_seq_len     [%lld]\n"
              /**/
              "     path_data_offset  [%lld]\n"
              "     seq_data_offset   [%lld]\n"
              "     seq_info_offset   [%lld]\n"
              "     total_file_length [%lld]\n",
              header->magic,
              header->version,
              header->type,
              header->line_length,
              /**/
              header->number_of_dbs,
              header->max_db_len,
              header->total_db_len,
              /**/
              header->number_of_seqs,
              header->max_seq_len,
              header->total_seq_len,
              /**/
              header->path_data_offset,
              header->seq_data_offset,
              header->seq_info_offset,
              header->total_file_length);
    return;
    }
#endif /* 0 */

static void Dataset_Header_write(Dataset_Header *header, FILE *fp){
    /* Dataset_Header_info(header); */
    /**/
    BitArray_write_int(header->magic, fp);
    BitArray_write_int(header->version, fp);
    BitArray_write_int(header->type, fp);
    BitArray_write_int(header->line_length, fp);
    /**/
    BitArray_write_int(header->number_of_dbs, fp);
    BitArray_write_int(header->max_db_len, fp);
    BitArray_write_int(header->total_db_len, fp);
    /**/
    BitArray_write_int(header->number_of_seqs, fp);
    BitArray_write_int(header->max_seq_len, fp);
    BitArray_write_int(header->total_seq_len, fp);
    /**/
    BitArray_write_int(header->path_data_offset, fp);
    BitArray_write_int(header->seq_data_offset, fp);
    BitArray_write_int(header->seq_info_offset, fp);
    BitArray_write_int(header->total_file_length, fp);
    /**/
    return;
    }

static Dataset_Header *Dataset_Header_read(FILE *fp){
    register Dataset_Header *header = g_new(Dataset_Header, 1);
    /**/
    header->magic = BitArray_read_int(fp);
    if(header->magic != DATASET_HEADER_MAGIC)
        g_error("Bad magic number in dataset file");
    header->version = BitArray_read_int(fp);
    if(header->version != DATASET_HEADER_VERSION)
        g_error("Incompatable dataset file version");
    header->type = BitArray_read_int(fp);
    header->line_length = BitArray_read_int(fp);
    /**/
    header->number_of_dbs = BitArray_read_int(fp);
    header->max_db_len = BitArray_read_int(fp);
    header->total_db_len = BitArray_read_int(fp);
    /**/
    header->number_of_seqs = BitArray_read_int(fp);
    header->max_seq_len = BitArray_read_int(fp);
    header->total_seq_len = BitArray_read_int(fp);
    /**/
    header->path_data_offset = BitArray_read_int(fp);
    header->seq_data_offset = BitArray_read_int(fp);
    header->seq_info_offset = BitArray_read_int(fp);
    header->total_file_length = BitArray_read_int(fp);
    /**/
    /* Dataset_Header_info(header); */
    return header;
    }

/**/

static Dataset_Sequence *Dataset_Sequence_create(FastaDB_Seq *fdbs){
    register Dataset_Sequence *seq = g_new0(Dataset_Sequence, 1);
    seq->key = FastaDB_Seq_get_key(fdbs);
    seq->gcg_checksum = Sequence_checksum(fdbs->seq);
    seq->id = g_strdup(fdbs->seq->id);
    seq->def = fdbs->seq->def?g_strdup(fdbs->seq->def):NULL;
    seq->cache_seq = NULL;
    return seq;
    }

static void Dataset_Sequence_destroy(Dataset_Sequence *seq){
    if(seq->cache_seq)
        Sequence_destroy(seq->cache_seq);
    FastaDB_Key_destroy(seq->key);
    g_free(seq->id);
    if(seq->def)
        g_free(seq->def);
    g_free(seq);
    return;
    }

static gsize Dataset_Sequence_memory_usage(Dataset_Sequence *ds){
    return sizeof(Dataset_Sequence)
         + sizeof(FastaDB_Key) /* FIXME: count FastaDB_Key memory properly */
         + (sizeof(gchar)*strlen(ds->id))
         + (sizeof(gchar)*(ds->def?strlen(ds->id):0))
         + (ds->cache_seq?Sequence_memory_usage(ds->cache_seq):0);
    }

static int Dataset_Sequence_compare_by_id_uniq(const void *a,
                                          const void *b){
    register Dataset_Sequence **seq_a = (Dataset_Sequence**)a,
                              **seq_b = (Dataset_Sequence**)b;
    register gint retval = strcmp((*seq_a)->id, (*seq_b)->id);
    if(!retval)
        g_error("Dataset has duplicate sequence id: [%s]",
                (*seq_a)->id);
    return retval;
    }

static int Dataset_Sequence_compare_by_id(const void *a,
                                          const void *b){
    register Dataset_Sequence **seq_a = (Dataset_Sequence**)a,
                              **seq_b = (Dataset_Sequence**)b;
    register gint retval = strcmp((*seq_a)->id, (*seq_b)->id);
    return retval;
    }

/**/

static Dataset_Width *Dataset_Width_create(Dataset_Header *header){
    register Dataset_Width *width = g_new0(Dataset_Width, 1);
    register guint64 n;
    for(n = header->number_of_dbs; n; n >>= 1)
        width->num_db_width++;
    for(n = header->max_db_len; n; n >>= 1)
        width->max_db_len_width++;
    for(n = header->max_seq_len; n; n >>= 1)
        width->max_seq_len_width++;
    n = width->num_db_width
      + width->max_db_len_width
      + width->max_seq_len_width
      + 14; /* for gcg checksum */
    width->seq_data_item_size = (n >> 3) + ((n % 8)?1:0);
    return width;
    }

static void Dataset_Width_destroy(Dataset_Width *width){
    g_free(width);
    return;
    }

/**/

Dataset *Dataset_create(GPtrArray *path_list,
                        Alphabet_Type alphabet_type, gboolean softmask_input){
    register Dataset *dataset = g_new(Dataset, 1);
    register FastaDB_Seq *fdbs;
    register Dataset_Sequence *ds;
    register gint i;
    register gsize path_data_size, seq_info_size, seq_data_size;
    register gchar *path;
    dataset->ref_count = 1;
    if(alphabet_type == Alphabet_Type_UNKNOWN){
        alphabet_type = FastaDB_guess_type((gchar*)path_list->pdata[0]);
        g_message("Guessed alphabet type as [%s]",
                  Alphabet_Type_get_name(alphabet_type));
        }
    dataset->alphabet = Alphabet_create(alphabet_type, softmask_input);
    dataset->header = Dataset_Header_create(alphabet_type, softmask_input);
    dataset->fdb = FastaDB_open_list(path_list, dataset->alphabet);
    /**/
    dataset->seq_list = g_ptr_array_new();
    while((fdbs = FastaDB_next(dataset->fdb, FastaDB_Mask_ALL))){
        ds = Dataset_Sequence_create(fdbs);
        g_ptr_array_add(dataset->seq_list, ds);
        FastaDB_Seq_destroy(fdbs);
        dataset->header->number_of_seqs++;
        if(dataset->header->max_seq_len < ds->key->length)
            dataset->header->max_seq_len = ds->key->length;
        dataset->header->total_seq_len += ds->key->length;
        }
    dataset->header->line_length = dataset->fdb->line_length;
    if(dataset->header->line_length < 1)
        g_error("Input is not a regular FASTA file, use fastareformat");
    dataset->header->number_of_dbs = path_list->len;
    dataset->header->max_db_len
        = CompoundFile_get_max_element_length(dataset->fdb->cf);
    dataset->header->total_db_len = CompoundFile_get_length(dataset->fdb->cf);
    qsort(dataset->seq_list->pdata, dataset->seq_list->len,
          sizeof(gpointer), Dataset_Sequence_compare_by_id_uniq);
    dataset->width = Dataset_Width_create(dataset->header);
    /**/
    path_data_size = 0;
    for(i = 0; i < path_list->len; i++){
        path = path_list->pdata[i];
        path_data_size += (strlen(path) + 1);
        }
    seq_data_size = 0;
    for(i = 0; i < dataset->seq_list->len; i++){
        ds = dataset->seq_list->pdata[i];
        ds->pos = i;
        seq_data_size += (strlen(ds->id) + 1);
        if(ds->def)
            seq_data_size += (strlen(ds->def) + 1);
        }
    seq_info_size = dataset->seq_list->len
                  * dataset->width->seq_data_item_size;
    /**/
    dataset->header->path_data_offset = sizeof(Dataset_Header);
    dataset->header->seq_data_offset = dataset->header->path_data_offset
                                     + path_data_size;
    dataset->header->seq_info_offset = dataset->header->seq_data_offset
                                     + seq_data_size;
    dataset->header->total_file_length = dataset->header->seq_info_offset
                                       + seq_info_size;
#ifdef USE_PTHREADS
    pthread_mutex_init(&dataset->dataset_mutex, NULL);
#endif /* USE_PTHREADS */
    return dataset;
    }

Dataset *Dataset_share(Dataset *dataset){
    dataset->ref_count++;
    return dataset;
    }

gsize Dataset_memory_usage(Dataset *dataset){
    register gint i;
    register gsize dataset_sequence_memory = 0;
    register Dataset_Sequence *ds;
    for(i = 0; i < dataset->seq_list->len; i++){
        ds = dataset->seq_list->pdata[i];
        dataset_sequence_memory += Dataset_Sequence_memory_usage(ds);
        }
    return sizeof(Dataset)
         + sizeof(Alphabet)
         + sizeof(Dataset_Header)
         + sizeof(Dataset_Width)
         + (sizeof(gpointer)*dataset->seq_list->len)
         + FastaDB_memory_usage(dataset->fdb)
         + dataset_sequence_memory;
    }

// allow compilation with "-Werror"
#pragma GCC diagnostic warning "-Wformat"
void Dataset_info(Dataset *dataset){
    g_print("Sequence Dataset:\n"
            "----------------\n"
            "             version: %d\n"
            "                type: %s, %s\n"
            "          total size: %d Mb\n"
            " number of sequences: %lld\n"
            "    longest sequence: %lld\n\n",
            (gint)dataset->header->version,
            (dataset->header->type & 1)?"DNA":"PROTEIN",
            (dataset->header->type & 2)?"masked":"unmasked",
            (gint)(dataset->header->total_db_len >> 20),
            dataset->header->number_of_seqs,
            dataset->header->max_seq_len);
    return;
    }

void Dataset_destroy(Dataset *dataset){
    register gint i;
    register Dataset_Sequence *seq;
    if(--dataset->ref_count)
        return;
    for(i = 0; i < dataset->seq_list->len; i++){
        seq = dataset->seq_list->pdata[i];
        Dataset_Sequence_destroy(seq);
        }
    FastaDB_close(dataset->fdb);
    g_ptr_array_free(dataset->seq_list, TRUE);
    Dataset_Header_destroy(dataset->header);
    if(dataset->width)
        Dataset_Width_destroy(dataset->width);
    Alphabet_destroy(dataset->alphabet);
#ifdef USE_PTHREADS
    pthread_mutex_destroy(&dataset->dataset_mutex);
#endif /* USE_PTHREADS */
    g_free(dataset);
    return;
    }

gboolean Dataset_check_filetype(gchar *path){
    register FILE *fp = fopen(path, "r");
    register guint64 magic;
    if(!fp)
        g_error("Could not open file [%s]", path);
    magic = BitArray_read_int(fp);
    fclose(fp);
    return (magic == DATASET_HEADER_MAGIC);
    }

static void Dataset_write_path_data(Dataset *dataset, FILE *fp){
    register gint i;
    register CompoundFile_Element *cfe;
    for(i = 0; i < dataset->fdb->cf->element_list->len; i++){
        cfe = dataset->fdb->cf->element_list->pdata[i];
        fprintf(fp, "%s\n", cfe->path);
        }
    return;
    }

static void Dataset_read_path_data(Dataset *dataset, FILE *fp){
    register gint i;
    gchar buf[1024];
    register GPtrArray *path_list = g_ptr_array_new();
    register gchar *path;
    for(i = 0; i < dataset->header->number_of_dbs; i++){
        if(!fgets(buf, 1024, fp))
            g_error("Problem parsing file data");
        path = g_strndup(buf, strlen(buf)-1);
        g_ptr_array_add(path_list, path);
        }
    dataset->fdb = FastaDB_open_list(path_list, dataset->alphabet);
    dataset->fdb->line_length = dataset->header->line_length;
    for(i = 0; i < path_list->len; i++){
        path = path_list->pdata[i];
        g_free(path);
        }
    g_ptr_array_free(path_list, TRUE);
    return;
    }

static void Dataset_write_seq_data(Dataset *dataset, FILE *fp){
    register gint i;
    register Dataset_Sequence *ds;
    for(i = 0; i < dataset->seq_list->len; i++){
        ds = dataset->seq_list->pdata[i];
        fprintf(fp, "%s%s%s\n",
                ds->id,
                ds->def?" ":"",
                ds->def?ds->def:"");
        }
    return;
    }

static void Dataset_read_seq_data(Dataset *dataset, FILE *fp){
    register gint64 i, j = 0, ipos, dpos;
    register Dataset_Sequence *ds;
    register guint64 len = dataset->header->seq_info_offset
                         - dataset->header->seq_data_offset;
    register gchar *buf = g_new(gchar, len+1);
    if(!fread(buf, sizeof(gchar), len, fp))
        g_error("Problem reading seq data");
    for(i = 0; i < dataset->header->number_of_seqs; i++){
        ds = g_new0(Dataset_Sequence, 1);
        g_ptr_array_add(dataset->seq_list, ds);
        ipos = j;
        dpos = -1;
        while(buf[j]){
            if((buf[j] == ' ') && (dpos == -1)){
                buf[j] = '\0';
                dpos = j+1;
                }
            if(buf[j] == '\n'){
                buf[j] = '\0';
                ds->id = g_strdup(buf+ipos);
                if(dpos != -1)
                    ds->def = g_strdup(buf+dpos);
                j++;
                break;
                }
            j++;
            }
        g_assert(ds->id);
        }
    g_free(buf);
    return;
    }

static void Dataset_write_seq_info(Dataset *dataset, FILE *fp){
    register gint i;
    register Dataset_Sequence *sequence;
    register BitArray *ba = BitArray_create();
    for(i = 0; i < dataset->seq_list->len; i++){
        sequence = dataset->seq_list->pdata[i];
        /**/
        BitArray_append(ba, sequence->key->location->element_id,
                        dataset->width->num_db_width);
        BitArray_append(ba, sequence->key->location->pos,
                        dataset->width->max_db_len_width);
        BitArray_append(ba, sequence->key->length,
                        dataset->width->max_seq_len_width);
        BitArray_append(ba, sequence->gcg_checksum, 14);
        BitArray_write(ba, fp);
        BitArray_empty(ba);
        }
    BitArray_destroy(ba);
    return;
    }

static void Dataset_read_seq_info(Dataset *dataset, FILE *fp){
    register gint i, start, element_id, length, seq_offset;
    register Dataset_Sequence *ds;
    register BitArray *ba;
    register CompoundFile_Pos offset;
    register Sequence_Strand strand
        = (dataset->alphabet->type == Alphabet_Type_DNA)
        ? Sequence_Strand_FORWARD
        : Sequence_Strand_UNKNOWN;
    register CompoundFile_Location *location;
    for(i = 0; i < dataset->seq_list->len; i++){
        ds = dataset->seq_list->pdata[i];
        ba = BitArray_read(fp, dataset->width->seq_data_item_size);
        start = 0;
        element_id = BitArray_get(ba, start, dataset->width->num_db_width);
        start += dataset->width->num_db_width;
        offset = BitArray_get(ba, start, dataset->width->max_db_len_width);
        start += dataset->width->max_db_len_width;
        length = BitArray_get(ba, start, dataset->width->max_seq_len_width);
        start += dataset->width->max_seq_len_width;
        location = CompoundFile_Location_create(dataset->fdb->cf,
                                                offset, element_id);
        seq_offset = 1 + strlen(ds->id)
                   + (ds->def?strlen(ds->def)+1:0);
        ds->key = FastaDB_Key_create(dataset->fdb,
                                     location, strand,
                                     seq_offset, length);
        CompoundFile_Location_destroy(location);
        ds->gcg_checksum = BitArray_get(ba, start, 14);
        ds->pos = i;
        BitArray_destroy(ba);
        }
    return;
    }

void Dataset_write(Dataset *dataset, gchar *path){
    register FILE *fp = fopen(path, "r");
    if(fp)
        g_error("Output file [%s] already exists", path);
    fp = fopen(path, "w");
    Dataset_Header_write(dataset->header, fp);
    Dataset_write_path_data(dataset, fp);
    Dataset_write_seq_data(dataset, fp);
    Dataset_write_seq_info(dataset, fp);
    fclose(fp);
    return;
    }

Dataset *Dataset_read(gchar *path){
    register Dataset *dataset = g_new(Dataset, 1);
    register FILE *fp = fopen(path, "r");
    if(!fp)
        g_error("Could not open esd file [%s]", path);
    dataset->ref_count = 1;
    dataset->header = Dataset_Header_read(fp);
    dataset->alphabet = Alphabet_create(
                           ((dataset->header->type&1)
                            ?Alphabet_Type_DNA
                            :Alphabet_Type_PROTEIN),
                            (dataset->header->type&2));
    dataset->width = Dataset_Width_create(dataset->header);
    dataset->seq_list = g_ptr_array_new();
    Dataset_read_path_data(dataset, fp);
    Dataset_read_seq_data(dataset, fp);
    Dataset_read_seq_info(dataset, fp);
#ifdef USE_PTHREADS
    pthread_mutex_init(&dataset->dataset_mutex, NULL);
#endif /* USE_PTHREADS */
    fclose(fp);
    return dataset;
    }

gint Dataset_lookup_id(Dataset *dataset, gchar *id){
    register Dataset_Sequence *result, **result_ptr;
    Dataset_Sequence key_seq, *key_ptr;
    key_seq.id = id;
    key_ptr = &key_seq;
    result_ptr = bsearch(&key_ptr, dataset->seq_list->pdata,
            dataset->seq_list->len,
            sizeof(gpointer), Dataset_Sequence_compare_by_id);
    if(!result_ptr)
        return -1;
    result = *result_ptr;
    return result->pos;
    }

Sequence *Dataset_get_sequence(Dataset *dataset, gint dataset_pos){
    register Sequence *seq = NULL;
    register Dataset_Sequence *ds;
    register gchar *def;
    register SparseCache *cache;
    register Sequence_Strand strand;
    g_assert(dataset_pos >= 0);
    g_assert(dataset_pos < dataset->seq_list->len);
    ds = dataset->seq_list->pdata[dataset_pos];
#ifdef USE_PTHREADS
    pthread_mutex_lock(&dataset->dataset_mutex);
#endif /* USE_PTHREADS */
    if(ds->cache_seq)
        seq = Sequence_share(ds->cache_seq);
#ifdef USE_PTHREADS
    pthread_mutex_unlock(&dataset->dataset_mutex);
#endif /* USE_PTHREADS */
    if(seq)
        return seq;
    def = FastaDB_Key_get_def(ds->key);
    strand = (dataset->alphabet->type == Alphabet_Type_DNA)
           ? Sequence_Strand_FORWARD
           : Sequence_Strand_UNKNOWN;
    cache = FastaDB_Key_get_SparseCache(ds->key);
    seq = Sequence_create_extmem(ds->id, def, ds->key->length,
                                 strand, dataset->alphabet, cache);
    g_free(def);
    SparseCache_destroy(cache);
    g_assert(Sequence_checksum(seq) == ds->gcg_checksum); /* FIXME: slow */
    return seq;
    }
/* FIXME: should keep a single compound file for whole dataset,
 *        and use that for making keys, to prevent more than
 *        one file being opened simultaneously.
 *        ?? will not work with cf sorting ... cf needs refactoring
 */

static int Dataset_Sequence_preload_compare(const void *a, const void *b){
    register Dataset_Sequence **ds_a = (Dataset_Sequence**)a,
                              **ds_b = (Dataset_Sequence**)b;
    g_assert(*ds_a);
    g_assert(*ds_b);
    g_assert((*ds_a)->key);
    g_assert((*ds_b)->key);
    if((*ds_a)->key->location->element_id
    == (*ds_b)->key->location->element_id)
        return (*ds_a)->key->location->pos
             - (*ds_b)->key->location->pos;
    return (*ds_a)->key->location->element_id
         - (*ds_b)->key->location->element_id;
    }

void Dataset_preload_seqs(Dataset *dataset){
    register gint i;
    register Dataset_Sequence *ds;
    register Dataset_Sequence **ds_list = g_new(Dataset_Sequence*,
                                                dataset->header->number_of_seqs);
    for(i = 0; i < dataset->header->number_of_seqs; i++){
        ds = dataset->seq_list->pdata[i];
        ds->cache_seq = Dataset_get_sequence(dataset, i);
        ds_list[i] = ds;
        }
    qsort(ds_list, dataset->header->number_of_seqs,
          sizeof(Dataset_Sequence*), Dataset_Sequence_preload_compare);
    /* Preload the sequence in the order they appear on disk */
    for(i = 0; i < dataset->header->number_of_seqs; i++){
        ds = ds_list[i];
        Sequence_preload_extmem(ds->cache_seq);
        }
    g_free(ds_list);
    return;
    }

/**/

