/****************************************************************\
*                                                                *
*  Library for manipulation of exonerate index files             *
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
#include <math.h>    /* For pow() */
#include <string.h>  /* For strlen() */
#include <stdlib.h>  /* For qsort() */

#include "index.h"
#include "submat.h"
#include "wordhood.h"
#include "pqueue.h"
#include "rangetree.h"
#include "noitree.h"

#define INDEX_HEADER_MAGIC (('e' << 16)|('s' << 8)|('i'))
#define INDEX_HEADER_VERSION 3

static off_t Index_ftell(FILE *fp){
    return ftello(fp);
    }
/* FIXME: need alternatives in absence of ftello()
 *        lseek(fileno(fp), 0, SEEK_CUR); ??
 */

static off_t Index_fseek(FILE *fp, off_t offset, int whence){
    return fseeko(fp, offset, whence);
    }
/* FIXME: need alternatives in absence of fseeko()
 *        lseek(fileno(fp), offset, whence); ??
 */

static guint64 Index_Strand_memory_usage(Index_Strand *strand, VFSM *vfsm){
    g_assert(strand);
    return sizeof(Index_Strand)
         + (sizeof(gint)*vfsm->lrw) /* for word_table */
         + (sizeof(Index_Word)*strand->header.word_list_length) /* word_list */
         + (strand->index_cache?BitArray_memory_usage(strand->index_cache):0);
    }

guint64 Index_memory_usage(Index *index){
    g_assert(index);
    return sizeof(Index)
         + sizeof(Index_Header)
         + (sizeof(gchar)*index->header->dataset_path_len)
         + Dataset_memory_usage(index->dataset)
         + sizeof(Index_Width)
         + sizeof(VFSM)
         + (index->forward?Index_Strand_memory_usage(index->forward,
                                                     index->vfsm):0)
         + (index->revcomp?Index_Strand_memory_usage(index->forward,
                                                     index->vfsm):0);
    }

static gsize Index_Word_get_address_list_space(Index_Word *index_word,
                                               Index *index){
    g_assert(index_word);
    if(index_word->freq_count)
        return BitArray_get_size(index_word->freq_count
                               *(index->width->number_of_seqs_width
                                +index->dataset->width->max_seq_len_width));
    return 0;
    }

/**/

static Index_Header *Index_Header_create(gint dataset_path_len,
                                         gboolean is_translated,
                                         gint word_length,
                                         gint word_jump,
                                         gint word_ambiguity,
                                         gint saturate_threshold){
    register Index_Header *index_header = g_new0(Index_Header, 1);
    index_header->magic = INDEX_HEADER_MAGIC;
    index_header->version = INDEX_HEADER_VERSION;
    index_header->type = is_translated?1:0;
    index_header->dataset_path_len = dataset_path_len;
    /**/
    index_header->word_length = word_length;
    index_header->word_jump = word_jump;
    index_header->word_ambiguity = word_ambiguity;
    index_header->saturate_threshold = saturate_threshold;
    return index_header;
    }

static void Index_Header_destroy(Index_Header *index_header){
    g_free(index_header);
    return;
    }

#if 0
static void Index_Header_info(Index_Header *index_header){
    g_message("Index_Header:\n"
              "    magic [%lld]\n"
              "    version [%lld]\n"
              "    type [%lld]\n"
              "    dataset_path_len [%lld]\n"
              "\n"
              "    word_length [%lld]\n"
              "    word_jump [%lld]\n"
              "    word_ambiguity [%lld]\n"
              "    saturate_threshold [%lld]\n",
           index_header->magic,
           index_header->version,
           index_header->type,
           index_header->dataset_path_len,
           index_header->word_length,
           index_header->word_jump,
           index_header->word_ambiguity,
           index_header->saturate_threshold);
    return;
    }
#endif /* 0 */

static void Index_Header_write(Index_Header *index_header, FILE *fp){
    /* Index_Header_info(index_header); */
    BitArray_write_int(index_header->magic, fp);
    BitArray_write_int(index_header->version, fp);
    BitArray_write_int(index_header->type, fp);
    BitArray_write_int(index_header->dataset_path_len, fp);
    /**/
    BitArray_write_int(index_header->word_length, fp);
    BitArray_write_int(index_header->word_jump, fp);
    BitArray_write_int(index_header->word_ambiguity, fp);
    BitArray_write_int(index_header->saturate_threshold, fp);
    /**/
    return;
    }

static Index_Header *Index_Header_read(FILE *fp){
    register Index_Header *index_header = g_new(Index_Header, 1);
    index_header->magic = BitArray_read_int(fp);
    if(index_header->magic != INDEX_HEADER_MAGIC)
        g_error("Bad magic number in index file");
    index_header->version = BitArray_read_int(fp);
    if(index_header->version != INDEX_HEADER_VERSION)
        g_error("Incompatible index file version");
    index_header->type = BitArray_read_int(fp);
    index_header->dataset_path_len = BitArray_read_int(fp);
    /**/
    index_header->word_length = BitArray_read_int(fp);
    index_header->word_jump  = BitArray_read_int(fp);
    index_header->word_ambiguity = BitArray_read_int(fp);
    index_header->saturate_threshold = BitArray_read_int(fp);
    /* Index_Header_info(index_header); */
    return index_header;
    }

/**/

static Index_Width *Index_Width_create(VFSM *vfsm, Dataset *dataset){
    register Index_Width *index_width = g_new0(Index_Width, 1);
    register guint64 n;
    for(n = vfsm->lrw; n; n >>= 1)
        index_width->max_word_width++;
    for(n = dataset->header->number_of_seqs; n; n >>= 1)
        index_width->number_of_seqs_width++;
    return index_width;
    }

#if 0
static void Index_Width_info(Index_Width *index_width){
    g_message("Index_Width:\n"
              "    max_word_width [%d]\n"
              "    number_of_seqs_width [%d]\n",
              index_width->max_word_width,
              index_width->number_of_seqs_width);
    return;
    }
#endif /* 0 */

static void Index_Width_destroy(Index_Width *index_width){
    g_free(index_width);
    return;
    }

/**/

typedef void (*Index_WordVisit_Func)(Index *index, Index_Strand *strand,
                                     gint seq_id,
                                     gint seq_pos, gint leaf_id,
                                     gpointer user_data);

static void Index_visit_seq_words_single(Index *index, Index_Strand *index_strand,
                                  gint seq_id, Sequence *seq,
                                  Index_WordVisit_Func iwvf, gint frame,
                                  gpointer user_data){
    register gint i, jump_ctr = 0, pos;
    register gchar *str = Sequence_get_str(seq);
    register VFSM_Int state = 0, leaf;
    for(i = 0; str[i]; i++){
        if(!index->vfsm->index[(guchar)str[i]]){
            state = 0;
            continue;
            }
        state = VFSM_change_state_M(index->vfsm, state, (guchar)str[i]);
        if(jump_ctr--)
            continue;
        jump_ctr = index->header->word_jump - 1;
        if(VFSM_state_is_leaf(index->vfsm, state)){
            leaf = VFSM_state2leaf(index->vfsm, state);
            pos = i-(index->vfsm->depth-1);
            if(frame)
                pos = (pos * 3) + frame - 1;
            iwvf(index, index_strand, seq_id, pos, leaf, user_data);
            }
        }
    g_free(str);
    return;
    }
/* FIXME: needs to work with softmasked sequences */

#ifndef Swap
#define Swap(x,y,temp) ((temp)=(x),(x)=(y),(y)=(temp))
#endif /* Swap */

static void Index_visit_seq_words_ambig(Index *index, Index_Strand *index_strand,
                                  gint seq_id, Sequence *seq,
                                  Index_WordVisit_Func iwvf, gint frame,
                                  gpointer user_data){
    register gint i, j, k, jump_ctr = 0, pos;
    register gchar *ambig, *str = Sequence_get_str(seq);
    register VFSM_Int state, next_state, leaf;
    register VFSM_Int *temp_state_list,
                      *curr_state_list = g_new0(VFSM_Int,
                                                index->header->word_ambiguity),
                      *next_state_list = g_new0(VFSM_Int,
                                                index->header->word_ambiguity);
    register gint curr_state_list_len = 1, next_state_list_len = 0;
    for(i = 0; str[i]; i++){
        if(jump_ctr--)
            continue;
        jump_ctr = index->header->word_jump - 1;
        if((!(ambig = Alphabet_nt2ambig(str[i])))
        || ((strlen(ambig) * curr_state_list_len)
             > index->header->word_ambiguity)){
            next_state_list_len = 0; /* Clear next state list */
            curr_state_list_len = 1; /* Set single null state */
            curr_state_list[0] = 0;
            continue;
            }
        for(j = 0; j < curr_state_list_len; j++){
            state = curr_state_list[j];
            for(k = 0; ambig[k]; k++){
                g_assert(index->vfsm->index[(guchar)ambig[k]]);
                next_state = VFSM_change_state_M(index->vfsm,
                                                 state, (guchar)ambig[k]);
                if(next_state_list_len
                && ((next_state == next_state_list[0])
                 || (next_state == next_state_list[next_state_list_len-1])))
                    break;
                next_state_list[next_state_list_len++] = next_state;
                if(VFSM_state_is_leaf(index->vfsm, next_state)){
                    leaf = VFSM_state2leaf(index->vfsm, next_state);
                    pos = i-(index->vfsm->depth-1);
                    if(frame)
                        pos = (pos * 3) + frame - 1;
                    iwvf(index, index_strand, seq_id, pos, leaf, user_data);
                    }
                }
            }
        curr_state_list_len = next_state_list_len;
        next_state_list_len = 0;
        Swap(next_state_list, curr_state_list, temp_state_list);
        }
    g_free(str);
    g_free(curr_state_list);
    g_free(next_state_list);
    return;
    }
/* FIXME: needs to work with softmasked sequences */
/* FIXME: optimisation: remove strlen() call */

static void Index_visit_seq_words(Index *index, Index_Strand *index_strand,
                                  gint seq_id, Sequence *seq,
                                  Index_WordVisit_Func iwvf, gint frame,
                                  gpointer user_data){
    if(index->header->word_ambiguity > 1)
        Index_visit_seq_words_ambig(index, index_strand, seq_id, seq,
                              iwvf, frame, user_data);
    else
        Index_visit_seq_words_single(index, index_strand, seq_id, seq,
                              iwvf, frame, user_data);
    return;
    }
/* FIXME: needs to work with softmasked sequences */

static void Index_visit_words(Index *index, Index_Strand *index_strand,
                              Index_WordVisit_Func iwvf,
                              gboolean is_forward, gpointer user_data){
    register gint i, j;
    register Dataset_Sequence *dataset_seq;
    register Sequence *ds_seq, *seq, *rc_seq, *aa_seq;
    register Translate *translate = Translate_create(FALSE);
    /* FIXME: move Translate_create() outside of this function */
    for(i = 0; i < index->dataset->seq_list->len; i++){
        dataset_seq = index->dataset->seq_list->pdata[i];
        ds_seq = Dataset_get_sequence(index->dataset, i);
        seq = Sequence_mask(ds_seq);
        Sequence_destroy(ds_seq);
        if(!is_forward){ /* revcomp */
            rc_seq = Sequence_revcomp(seq);
            Sequence_destroy(seq);
            seq = rc_seq;
            }
        if(index->header->type & 1){ /* is_translated */
            for(j = 0; j < 3; j++){
                aa_seq = Sequence_translate(seq, translate, j+1);
                Index_visit_seq_words(index, index_strand,
                                      i, aa_seq, iwvf, j+1, user_data);
                Sequence_destroy(aa_seq);
                }
        } else {
            Index_visit_seq_words(index, index_strand,
                                  i, seq, iwvf, 0, user_data);
            }
        Sequence_destroy(seq);
        }
    Translate_destroy(translate);
    return;
    }

/**/

static void Index_freq_count_visit(Index *index, Index_Strand *index_strand,
                                   gint seq_id, gint seq_pos,
                                   gint leaf_id, gpointer user_data){
    g_assert(leaf_id < index->vfsm->lrw);
    index_strand->word_table[leaf_id]++;
    return;
    }

static void Index_freq_count(Index *index, Index_Strand *index_strand,
                             gboolean is_forward){
    Index_visit_words(index, index_strand, Index_freq_count_visit,
                      is_forward, NULL);
    return;
    }

/**/

static gint Index_get_expect(Index *index, Index_Strand *index_strand){
    register gdouble expect = 1.0 / pow(index->vfsm->alphabet_size,
                                        index->header->word_length);
    register gint i;
    register gint64 observed_words = 0;
    for(i = 0; i < index->vfsm->lrw; i++)
        observed_words += index_strand->word_table[i];
    return (expect * observed_words) + index->header->saturate_threshold;
    }

// allow compilation with "-Werror"
#pragma GCC diagnostic warning "-Wformat"
static void Index_desaturate(Index *index, Index_Strand *index_strand){
    register VFSM_Int i;
    register gint64 removed_words = 0, removed_instances = 0;
    register gint target_expect = Index_get_expect(index, index_strand);
    if(!index->header->saturate_threshold)
        return;
    g_message("Desaturating (expecting [%d] target words)", target_expect);
    for(i = 0; i < index->vfsm->lrw; i++)
        if(index_strand->word_table[i] >= target_expect){
            removed_words++;
            removed_instances += index_strand->word_table[i];
            index_strand->word_table[i] = -1;
            }
    g_message("Removed [%" CUSTOM_GUINT64_FORMAT
              "] words, [%" CUSTOM_GUINT64_FORMAT "] instances",
              removed_words, removed_instances);
    return;
    }
/* state at function exit:
 * count for present words
 * zero for absent words
 * -1 for desaturated words
 */

static void Index_survey_word_list(Index *index,
                                   Index_Strand *index_strand){
    register gint i, pos = 0;
    register Index_Word *word;
    for(i = 0; i < index->vfsm->lrw; i++){
        if(index_strand->word_table[i] > 0)
            index_strand->header.word_list_length++;
        }
    index_strand->word_list
        = g_new(Index_Word, index_strand->header.word_list_length);
    for(i = 0; i < index->vfsm->lrw; i++){
        if(index_strand->word_table[i] > 0){
            word = &index_strand->word_list[pos];
            word->freq_count = index_strand->word_table[i];
            word->index_offset = 0;
            if(index_strand->header.max_index_length < word->freq_count)
                index_strand->header.max_index_length = word->freq_count;
            index_strand->word_table[i] = pos++;
        } else {
            if(index_strand->word_table[i] == 0)
                index_strand->word_table[i] = -2;
            }
        }
    return;
    }
/* state at function exit:
 * word_list offset for present
 * -1 for desaturated words
 * -2 for absent words
 */

gboolean Index_check_filetype(gchar *path){
    register FILE *fp = fopen(path, "r");
    register guint64 magic;
    if(!fp)
        g_error("Could not open file [%s]", path);
    magic = BitArray_read_int(fp);
    fclose(fp);
    return (magic == INDEX_HEADER_MAGIC);
    }

/**/

typedef struct {
    Index_Address *address_list;
             gint  found;
} Index_AddressList;

static Index_AddressList *Index_AddressList_create(gint freq_count){
    register Index_AddressList *address_list = g_new(Index_AddressList, 1);
    address_list->found = 0;
    address_list->address_list = g_new(Index_Address, freq_count);
    return address_list;
    }

static void Index_AddressList_destroy(Index_AddressList *address_list){
    g_free(address_list->address_list);
    g_free(address_list);
    return;
    }

static void Index_AddressList_append(Index_AddressList *address_list,
                                     gint seq_id, gint position){
    register Index_Address *address
        = &address_list->address_list[address_list->found];
    address->sequence_id = seq_id;
    address->position = position;
    address_list->found++;
    return;
    }

static void Index_AddressList_write(Index_AddressList *address_list,
                                    Index *index){
    register gint i;
    register BitArray *ba = BitArray_create();
    register Index_Address *address;
    for(i = 0; i < address_list->found; i++){
        address = &address_list->address_list[i];
        BitArray_append(ba, address->sequence_id,
                        index->width->number_of_seqs_width);
        BitArray_append(ba, address->position,
                        index->dataset->width->max_seq_len_width);
        }
    BitArray_write(ba, index->fp);
    BitArray_destroy(ba);
    return;
    }
/* FIXME: optimisation: clear and reuse BitArray */

static int Index_Address_compare(const void *a, const void *b){
    register Index_Address *addr_a = (Index_Address*)a,
                           *addr_b = (Index_Address*)b;
    if(addr_a->sequence_id == addr_b->sequence_id)
        return addr_a->position - addr_b->position;
    return addr_a->sequence_id - addr_b->sequence_id;
    }

static void Index_AddressList_sort(Index_AddressList *address_list){
    qsort(address_list->address_list, address_list->found,
          sizeof(Index_Address), Index_Address_compare);
    return;
    }
/* FIXME: optimisation: avoid this sort, by traversing all frames together,
 *                      or using more efficient 3-way merge
 */

/**/

static GArray *Index_find_pass_boundaries(Index *index,
                                          Index_Strand *index_strand,
                                          gsize available_space){
    register GArray *pass_boundary_list = g_array_new(FALSE, FALSE, sizeof(gint));
    register gsize usage = 0;
    gint i;
    register Index_Word *index_word;
    for(i = 0; i < index_strand->header.word_list_length; i++){
        index_word = &index_strand->word_list[i];
        usage += sizeof(gpointer)
               + sizeof(Index_AddressList)
               + (sizeof(Index_Address)*index_word->freq_count);
        if(usage >= available_space){
            g_array_append_val(pass_boundary_list, i);
            usage = 0;
            }
        }
    g_array_append_val(pass_boundary_list,
                       index_strand->header.word_list_length);
    return pass_boundary_list;
    }

/*
 *  0123-4567-89
 *  1234 5678 90
 */

/**/

typedef struct {
                  gint   first_word;
                  gint   interval_len;
     Index_AddressList **address_list_array;
} Index_AddressData;

static Index_AddressData *Index_AddressData_create(Index *index,
                                                   Index_Strand *index_strand,
                                                   gint first_word,
                                                   gint interval_len){
    register Index_AddressData *address_data = g_new(Index_AddressData, 1);
    register gint i;
    register Index_Word *index_word;
    address_data->first_word = first_word;
    address_data->interval_len = interval_len;
    address_data->address_list_array
                  = g_new(Index_AddressList*, interval_len);
    for(i = 0; i < interval_len; i++){
        index_word = &index_strand->word_list[first_word+i];
        address_data->address_list_array[i]
            = Index_AddressList_create(index_word->freq_count);
        }
    return address_data;
    }
/* FIXME: optimisation: use single alloc for all address_lists
 */

static void Index_AddressData_destroy(Index_AddressData *address_data){
    register gint i;
    for(i = 0; i < address_data->interval_len; i++)
        Index_AddressList_destroy(address_data->address_list_array[i]);
    g_free(address_data->address_list_array);
    g_free(address_data);
    return;
    }

static void Index_AddressData_write(Index_AddressData *address_data,
                                    Index *index, Index_Strand *index_strand){
    register gint i, word_id;
    register Index_Word *index_word;
    register Index_AddressList *address_list;
    for(i = 0; i < address_data->interval_len; i++){
        word_id = address_data->first_word+i;
        index_word = &index_strand->word_list[word_id];
        address_list = address_data->address_list_array[i];
        g_assert(index_word->freq_count == address_list->found);
        if(index->header->type & 1) /* sort address list when is_translated */
            Index_AddressList_sort(address_list);
        Index_AddressList_write(address_list, index);
        }
    return;
    }

/**/

static void Index_report_words_visit(Index *index, Index_Strand *index_strand,
                                     gint seq_id, gint seq_pos,
                                     gint leaf_id, gpointer user_data){
    register gint word_id = index_strand->word_table[leaf_id];
    register Index_Word *index_word;
    register Index_AddressData *address_data = user_data;
    register Index_AddressList *address_list;
    /* Skip word which are outside the interval for this pass */
    if((word_id < address_data->first_word)
    || (word_id >= (address_data->first_word+address_data->interval_len)))
        return;
    index_word = &index_strand->word_list[word_id];
    address_list
        = address_data->address_list_array[word_id-address_data->first_word];
    Index_AddressList_append(address_list, seq_id, seq_pos);
    return;
    }
/* FIXME: optimisation : replace allocation with RecycleBins */

static void Index_report_word_list(Index *index, Index_Strand *index_strand,
                                   gboolean is_forward, gint memory_limit){
    register Index_AddressList **address_list_array
        = g_new0(Index_AddressList*, index_strand->header.word_list_length);
    register guint64 available_memory = ((guint64)memory_limit << 20); /* Mb */
    register guint64 used_memory = Index_memory_usage(index)
                                 + Index_Strand_memory_usage(index_strand,
                                                             index->vfsm);
    register GArray *pass_boundary_list;
    register gint i, curr, prev = 0;
    register Index_AddressData *address_data;
    /**/
    if(used_memory > available_memory)
        g_error("Memory limit (%d Mb) already exceeded usage:[%d Mb]",
                memory_limit, (gint)(used_memory >> 20));
    g_message("Currently using [%d Mb] of [%d Mb]",
            (gint)(used_memory >> 20), memory_limit);
    /* Record strand offset */
    index_strand->strand_offset = Index_ftell(index->fp);
    /**/
    pass_boundary_list = Index_find_pass_boundaries(index, index_strand,
                                               (available_memory-used_memory));
    g_message("Using [%d] %s pass%s",
            pass_boundary_list->len,
            is_forward?"forward":"revcomp",
            (pass_boundary_list->len == 1)?"":"es");
    for(i = 0; i < pass_boundary_list->len; i++){
        curr = g_array_index(pass_boundary_list, gint, i);
        g_message("Sequence pass [%d] for word interval [%d] to [%d]",
                  i+1, prev, curr);
        address_data = Index_AddressData_create(index, index_strand,
                                                prev, curr-prev);
        Index_visit_words(index, index_strand,
                          Index_report_words_visit,
                          is_forward, address_data);
        Index_AddressData_write(address_data, index, index_strand);
        Index_AddressData_destroy(address_data);
        prev = curr;
        }
    g_array_free(pass_boundary_list, TRUE);
    g_free(address_list_array);
    return;
    }
/* Revised algorithm:
 *     Calculate how much space can be used
 *     Calculate how many passes are required
 *     for each pass
 *         collect words within word interval
 *         write out words
 */

static gsize Index_get_word_data_space(Index *index, Index_Strand *index_strand){
    register gint word_data_bits = index->width->max_word_width
                                 + index_strand->width.max_index_len_width
                                 + index_strand->width.total_index_len_width;
    return BitArray_get_size(word_data_bits
                           * index_strand->header.word_list_length);
    }

static void Index_find_offsets(Index *index, Index_Strand *index_strand){
    register gint i;
    register gint64 offset = 0;
    register Index_Word *index_word;
    for(i = 0; i < index_strand->header.word_list_length; i++){
        index_word = &index_strand->word_list[i];
        index_word->index_offset = offset;
        offset += Index_Word_get_address_list_space(index_word, index);
        }
    index_strand->header.total_index_length = offset;
    return;
    }

static void Index_Strand_Header_fill(Index_Strand_Header *header, FILE *fp){
    header->max_index_length = BitArray_read_int(fp);
    header->word_list_length = BitArray_read_int(fp);
    header->total_index_length = BitArray_read_int(fp);
    return;
    }

static void Index_Strand_Header_print(Index_Strand_Header *header, FILE *fp){
    BitArray_write_int(header->max_index_length, fp);
    BitArray_write_int(header->word_list_length, fp);
    BitArray_write_int(header->total_index_length, fp);
    return;
    }

static void Index_Strand_write(Index *index, Index_Strand *index_strand){
    register BitArray *word_ba = BitArray_create();
    register gint i, word_id, word_count = 0;
    register Index_Word *index_word;
    /* Create a bitarray from word_list and write it out */
    Index_Strand_Header_print(&index_strand->header, index->fp);
    for(i = 0; i < index->vfsm->lrw; i++){
        word_id = index_strand->word_table[i];
        if(word_id >= 0){
            BitArray_append(word_ba, i,
                            index->width->max_word_width);
            index_word = &index_strand->word_list[word_id];
            BitArray_append(word_ba, index_word->freq_count,
                            index_strand->width.max_index_len_width);
            BitArray_append(word_ba, index_word->index_offset,
                            index_strand->width.total_index_len_width);
            word_count++;
            }
        }
    g_assert(word_count == index_strand->header.word_list_length);
    BitArray_write(word_ba, index->fp);
    BitArray_destroy(word_ba);
    return;
    }

static void Index_Strand_Width_fill(Index_Strand *index_strand){
    register guint64 n;
    for(n = index_strand->header.max_index_length; n; n >>= 1)
        index_strand->width.max_index_len_width++;
    for(n = index_strand->header.total_index_length; n; n >>= 1)
        index_strand->width.total_index_len_width++;
    return;
    }

static Index_Strand *Index_Strand_create(Index *index, gboolean is_forward,
                                         gint memory_limit){
    register Index_Strand *index_strand = g_new0(Index_Strand, 1);
    g_assert(index);
    index_strand->word_table = g_new0(gint, index->vfsm->lrw);
    g_assert(index_strand->word_table);
    /* 1st sequence pass */
    Index_freq_count(index, index_strand, is_forward);
    /* Remove words occuring too frequently */
    Index_desaturate(index, index_strand);
    Index_survey_word_list(index, index_strand);
    Index_find_offsets(index, index_strand);
    Index_Strand_Width_fill(index_strand);
    Index_Strand_write(index, index_strand);
    /* 2nd seq pass */
    Index_report_word_list(index, index_strand, is_forward, memory_limit);
    index_strand->index_cache = NULL;
    return index_strand;
    }

static void Index_Strand_destroy(Index_Strand *index_strand){
    g_assert(index_strand);
    g_free(index_strand->word_table);
    g_free(index_strand->word_list);
    if(index_strand->index_cache)
        BitArray_destroy(index_strand->index_cache);
    g_free(index_strand);
    return;
    }

Index *Index_create(Dataset *dataset, gboolean is_translated, gint word_length,
                    gint word_jump, gint word_ambiguity, gint saturate_threshold,
                    gchar *index_path, gchar *dataset_path, gint memory_limit){
    register Index *index = g_new0(Index, 1);
    register gchar *member;
    register Alphabet *alphabet;
    index->ref_count = 1;
    /* Open index path for reading */
    index->fp = fopen(index_path, "r");
    if(index->fp){
        fclose(index->fp);
        g_error("Index file [%s] already exists", index_path);
        }
    /* Open index path for writing */
    index->fp = fopen(index_path, "w");
    if(!index->fp)
        g_error("Could not open [%s] to write index file", index_path);
    index->dataset_path = g_strdup(dataset_path);
    index->dataset = Dataset_share(dataset);
    index->header = Index_Header_create(strlen(dataset_path),
                                        is_translated, word_length,
                                        word_jump, word_ambiguity,
                                        saturate_threshold);
    if(is_translated){
        alphabet = Alphabet_create(Alphabet_Type_PROTEIN, FALSE);
        member = (gchar*)alphabet->member;
        Alphabet_destroy(alphabet);
    } else {
        member = (gchar*)dataset->alphabet->member;
        }
    index->vfsm = VFSM_create(member, word_length);
    index->width = Index_Width_create(index->vfsm, dataset);
    /* Index_Width_info(index->width); */
    /* Write out header and dataset path */
    Index_Header_write(index->header, index->fp);
    g_assert(Index_ftell(index->fp) == sizeof(Index_Header));
    fwrite(index->dataset_path, sizeof(gchar),
           index->header->dataset_path_len, index->fp);
    g_assert(Index_ftell(index->fp)
        == (sizeof(Index_Header) + index->header->dataset_path_len));
    /* Create forward and revcomp Index_Strand */
    index->forward = Index_Strand_create(index, TRUE, memory_limit);
    if(is_translated)
        index->revcomp = Index_Strand_create(index, FALSE, memory_limit);
    /* Close index->fp and reopen read-only */
    fclose(index->fp);
    index->fp = fopen(index_path, "r");
    if(!index->fp)
        g_error("Could not repoen [%s] for reading", index_path);
#ifdef USE_PTHREADS
    pthread_mutex_init(&index->index_mutex, NULL);
#endif /* USE_PTHREADS */
    return index;
    }

/* Indexing algorithm:
 *
 * word_table: +ve for word_list, -ve for absent_word_list
 *
 * data:[absent_list|neighbour_list|address_list,address_count]
 *
 * - count word frequencies in word_table (SEQ PASS 1)
 * - desaturate word_table
 * - convert word_table to pos in word_list and fill freq_count
 * - calculate neighbours (WORDHOOD TRAVERSE)
 *     fill neighbour_count and store list in build_data
 *     when freq_count == 0
 *         add AbsentWord to build_data
 * - expand word_list with absent words
 *     (update word_list ordering in word_table)
 *  (the data array is now empty)
 * - write out the data for each word
 * - write out and free neighbour lists and padding for index addresses
 *
 * - parse seqs again, on each word: (SEQ PASS 2)
 *       fill Index_Data (in data)
 *       when (index_len == freq_count)
 *           write the Index_Data and free
 */

/* Algorithm overview:
 *
 * +count word freqs
 * +desaturate words
 * +convert word freqs to pos in word_list and fill index_word->freq_count
 * +survey word neighbourhood
 * +expand word_list with absent words
 * +Index_Width_create()
 * +write out zero-padded file
 * +calc and write out neighbourhoods
 * +calc and write out word addresses
 */

Index *Index_share(Index *index){
    index->ref_count++;
    return index;
    }

void Index_destroy(Index *index){
    if(--index->ref_count)
        return;
    if(index->fp)
        fclose(index->fp);
    g_free(index->dataset_path);
    Dataset_destroy(index->dataset);
    Index_Header_destroy(index->header);
    VFSM_destroy(index->vfsm);
    Index_Width_destroy(index->width);
    if(index->forward)
        Index_Strand_destroy(index->forward);
    if(index->revcomp)
        Index_Strand_destroy(index->revcomp);
#ifdef USE_PTHREADS
    pthread_mutex_destroy(&index->index_mutex);
#endif /* USE_PTHREADS */
    g_free(index);
    return;
    }

void Index_info(Index *index){
    g_print("Sequence Index:\n"
            "--------------\n"
            "                version: %d\n"
            "                   type: %s\n"
            "            word length: %d\n"
            "              word jump: %d\n"
            "         word ambiguity: %d\n"
            "     saturate threshold: %d\n\n",
            (gint)index->header->version,
            (index->header->type & 1) ? "translated" : "normal",
            (gint)index->header->word_length,
            (gint)index->header->word_jump,
            (gint)index->header->word_ambiguity,
            (gint)index->header->saturate_threshold);
    return;
    }

static Index_Word *Index_Strand_read_word_list(Index *index,
                                               Index_Strand *index_strand){
    register Index_Word *word_list = g_new(Index_Word,
                                           index_strand->header.word_list_length);
    register BitArray *word_list_ba = BitArray_read(index->fp,
                            Index_get_word_data_space(index, index_strand));
    register gint i, start = 0;
    register guint64 leaf;
    register Index_Word *index_word;
    for(i = 0; i < index->vfsm->lrw; i++)
        index_strand->word_table[i] = -1;
    for(i = 0; i < index_strand->header.word_list_length; i++){
        leaf = BitArray_get(word_list_ba, start,
                            index->width->max_word_width);
        start += index->width->max_word_width;
        index_strand->word_table[leaf] = i;
        index_word = &word_list[i];
        g_assert(index_word);
        /**/
        index_word->freq_count = BitArray_get(word_list_ba, start,
                               index_strand->width.max_index_len_width);
        start += index_strand->width.max_index_len_width;
        /**/
        index_word->index_offset = BitArray_get(word_list_ba, start,
                               index_strand->width.total_index_len_width);
        start += index_strand->width.total_index_len_width;
        }
    BitArray_destroy(word_list_ba);
    return word_list;
    }

static Index_Strand *Index_Strand_read(Index *index, gboolean is_forward){
    register Index_Strand *index_strand = g_new0(Index_Strand, 1);
    Index_Strand_Header_fill(&index_strand->header, index->fp);
    Index_Strand_Width_fill(index_strand);
    /**/
    index_strand->word_table = g_new(gint, index->vfsm->lrw);
    index_strand->word_list = Index_Strand_read_word_list(index, index_strand);
    index_strand->index_cache = NULL;
    /**/
    index_strand->strand_offset = Index_ftell(index->fp);
    /* Seek to end of strand index */
    Index_fseek(index->fp, index_strand->header.total_index_length, SEEK_CUR);
    return index_strand;
    }

Index *Index_open(gchar *index_path){
    register Index *index = g_new(Index, 1);
    register gchar *member;
    register Alphabet *alphabet;
    index->ref_count = 1;
    index->fp = fopen(index_path, "r");
    if(!index->fp)
        g_error("Could not open index [%s]", index_path);
    index->header = Index_Header_read(index->fp);
    index->dataset_path = g_new(gchar, index->header->dataset_path_len+1);
    index->dataset_path[index->header->dataset_path_len] = '\0';
    fread(index->dataset_path, sizeof(gchar),
          index->header->dataset_path_len, index->fp);
    index->dataset = Dataset_read(index->dataset_path);
    /**/
    if(index->header->type & 1){ /* is_translated */
        alphabet = Alphabet_create(Alphabet_Type_PROTEIN, FALSE);
        member = (gchar*)alphabet->member;
        Alphabet_destroy(alphabet);
    } else {
        member = (gchar*)index->dataset->alphabet->member;
        }
    index->vfsm = VFSM_create(member, index->header->word_length);
    index->width = Index_Width_create(index->vfsm, index->dataset);
    /* Index_Width_info(index->width); */
    /**/
    index->forward = Index_Strand_read(index, TRUE);
    if(index->header->type == 1) /* is_translated */
        index->revcomp = Index_Strand_read(index, FALSE);
    else
        index->revcomp = NULL;
#ifdef USE_PTHREADS
    pthread_mutex_init(&index->index_mutex, NULL);
#endif /* USE_PTHREADS */
    return index;
    }

typedef struct {
VFSM_Int leaf;
    gint query_pos;
} Index_WordSeed;

static void Index_get_query_word_list(Index *index, Index_Strand *index_strand,
                                      Sequence *query, GArray *word_seed_list,
                                      gint frame, HSP_Param *hsp_param){
    register gint i;
    register gchar *str = Sequence_get_str(query);
    register VFSM_Int state = 0;
    Index_WordSeed seed;
    for(i = 0; str[i]; i++){
        if(!index->vfsm->index[(guchar)str[i]]){
            state = 0;
            continue;
            }
        state = VFSM_change_state_M(index->vfsm, state, (guchar)str[i]);
        if(VFSM_state_is_leaf(index->vfsm, state)){
            seed.leaf = VFSM_state2leaf(index->vfsm, state);
            if((index_strand->word_table[seed.leaf] >= 0) || hsp_param->wordhood){
                seed.query_pos = i - (index->header->word_length-1);
                if(frame)
                    seed.query_pos = (seed.query_pos * 3) + frame - 1;
                g_array_append_val(word_seed_list, seed);
                }
            }
        }
    g_free(str);
    return;
    }
/* FIXME: needs to work with softmasked sequences */

typedef struct {
    gint start;
    gint length;
    gint target_id;
} Index_Interval;

static gboolean Index_Address_lt_Interval_start(Index_Address *address,
                                                Index_Interval *interval){
    if(address->sequence_id == interval->target_id)
        return address->position < interval->start;
    return address->sequence_id < interval->target_id;
    }

static gboolean Index_Address_lt_Interval_end(Index_Address *address,
                                              Index_Interval *interval){
    if(address->sequence_id == interval->target_id)
        return address->position < (interval->start + interval->length);
    return address->sequence_id < interval->target_id;
    }

static gboolean Index_Address_gt_Interval_end(Index_Address *address,
                                              Index_Interval *interval){
    if(address->sequence_id == interval->target_id)
        return address->position > (interval->start + interval->length);
    return address->sequence_id > interval->target_id;
    }

static gint Index_Address_list_trim_start(Index_Address *address_list,
                                          gint address_list_start,
                                          gint address_list_len,
                                          Index_Interval *limit){
    register gint left = address_list_start,
                  right = address_list_start+address_list_len-1, mid;
    while(left < right){
        mid = (left+right) >> 1;
        if(Index_Address_lt_Interval_start(&address_list[mid], limit)){
            if((mid == (address_list_start+address_list_len))
            || (!Index_Address_lt_Interval_start(&address_list[mid+1], limit)))
                return mid+1-address_list_start;
            left = mid + 1;
        } else {
            right = mid;
            }
        }
    return address_list_len;
    }

static gint Index_Address_list_trim_end(Index_Address *address_list,
                                        gint address_list_start,
                                        gint address_list_len,
                                        Index_Interval *limit){
    register gint left = address_list_start,
                  right = address_list_start+address_list_len-1, mid;
    while(left < right){
        mid = (left+right) >> 1;
        if(Index_Address_lt_Interval_end(&address_list[mid], limit)){
            left = mid + 1;
        } else {
            if((mid == address_list_start)
            || (Index_Address_lt_Interval_end(&address_list[mid-1], limit)))
                return address_list_start+address_list_len-mid;
            right = mid;
            }
        }
    return address_list_len;
    }

static void Index_Address_list_refine_recur(
            Index_Address *address_list,
            gint address_list_start, gint address_list_len,
            Index_Interval *interval_list,
            gint interval_list_start, gint interval_list_len,
            Index_Address *refined_address_list, gint *refined_address_list_len){
    register gint i, move, interval_mid;
    if(Index_Address_lt_Interval_start(&address_list[address_list_start],
                                       &interval_list[interval_list_start])){
        move = Index_Address_list_trim_start(address_list, address_list_start,
                                       address_list_len,
                                       &interval_list[interval_list_start]);
        address_list_start += move;
        address_list_len -= move;
        }
    if(Index_Address_gt_Interval_end(
                &address_list[address_list_start+address_list_len-1],
                &interval_list[interval_list_start+interval_list_len-1])){
        move = Index_Address_list_trim_end(address_list,
            address_list_start, address_list_len,
            &interval_list[interval_list_start+interval_list_len-1]);
        address_list_len -= move;
        }
    if(!address_list_len)
        return;
    if(interval_list_len == 1){
        for(i = 0; i < address_list_len; i++){
            refined_address_list[(*refined_address_list_len)].sequence_id
                = address_list[address_list_start+i].sequence_id;
            refined_address_list[(*refined_address_list_len)].position
                = address_list[address_list_start+i].position;
            (*refined_address_list_len)++;
            }
    } else {
        interval_mid = interval_list_len >> 1;
        if(interval_mid)
            Index_Address_list_refine_recur(address_list, address_list_start,
                                 address_list_len, interval_list,
                                 interval_list_start, interval_mid,
                                 refined_address_list, refined_address_list_len);
        if(interval_list_len-interval_mid)
            Index_Address_list_refine_recur(address_list, address_list_start,
                                 address_list_len, interval_list,
                                 interval_list_start+interval_mid,
                                 interval_list_len-interval_mid,
                                 refined_address_list, refined_address_list_len);
        }
    return;
    }
/* Algorithm: search address_list
 *
 *      input:
 *       address_list (ordered positions)
 *       interval_list (ordered intervals)
 *
 *     output:
 *       list of regions in address_list
 *
 *  algorithm:
 *      get_addressses(address_list, region_list)
 *          if(address_list starts before region_list)
 *              trim address_list start to region list start
 *          if(address_list ends after region_list)
 *              trim address_list end to region list end
 *          if(address list is empty)
 *              return
 *          if(region_list->len == 1){
 *              add address_list to output_list
 *          } else {
 *              region_mid = region+_list_len / 2
 *              get_addressses(address_list, region_start, region_mid);
 *              get_addresses(address_list, region_mid, region_end);
 *              }
 *
 */

/**/

static gint Index_Address_list_refine(Index_Address *address_list,
                                      gint address_list_len,
                                      GArray *interval_list){
    register Index_Address *refined_address_list
        = g_new(Index_Address, address_list_len);
    register gint i;
    gint refined_address_list_len = 0;
    Index_Address_list_refine_recur(address_list, 0, address_list_len,
            (Index_Interval*)interval_list->data, 0, interval_list->len,
            refined_address_list, &refined_address_list_len);
    g_assert(refined_address_list_len <= address_list_len);
    /* Copy refined_address_list to address_list */
    /* FIXME: optimisation: replace this */
    for(i = 0; i < refined_address_list_len; i++){
        address_list[i].sequence_id = refined_address_list[i].sequence_id;
        address_list[i].position = refined_address_list[i].position;
        }
    g_free(refined_address_list);
    return refined_address_list_len;
    }
/* Index_Address_list_refine shortens address_list in place,
 * returning the new_address_list_len, which may be zero
 * if there is no intersection between the address_list and interval_list
 *
 * input: interval_list {target_id, start, length}
 *          addresslist {sequence_id, position}
 */

static gint Index_Word_get_address_list(Index *index,
                                        Index_Strand *index_strand,
                                        Index_Word *index_word,
                                        Index_Address *address_list,
                                        GArray *interval_list){
    register BitArray *ba;
    register gint i, pos = 0;
    register guint64 start = 0;
    register gint address_list_len = index_word->freq_count;
    g_assert(index_word->freq_count);
    if(index_strand->index_cache){
        ba = index_strand->index_cache;
        start = index_word->index_offset * CHAR_BIT;
    } else {
#ifdef USE_PTHREADS
        pthread_mutex_lock(&index->index_mutex);
#endif /* USE_PTHREADS */
        Index_fseek(index->fp,
                    index_strand->strand_offset+index_word->index_offset,
                    SEEK_SET);
        ba = BitArray_read(index->fp,
                      Index_Word_get_address_list_space(index_word, index));
#ifdef USE_PTHREADS
        pthread_mutex_unlock(&index->index_mutex);
#endif /* USE_PTHREADS */
        }
    for(i = 0; i < index_word->freq_count; i++){
        address_list[pos].sequence_id = BitArray_get(ba, start,
                                         index->width->number_of_seqs_width);
        start += index->width->number_of_seqs_width;
        address_list[pos++].position = BitArray_get(ba, start,
                                       index->dataset->width->max_seq_len_width);
        start += index->dataset->width->max_seq_len_width;
        }
    if(!index_strand->index_cache)
        BitArray_destroy(ba);
    if(interval_list)
        address_list_len = Index_Address_list_refine(
                                  address_list, index_word->freq_count,
                                  interval_list);
    return address_list_len;
    }
/* FIXME: optimisation: use BitArray in pq_seed, to avoid duplication */

typedef struct {
           Index  *index;
    Index_Strand  *index_strand;
  Index_WordSeed   seed;
          GArray  *word_seed_list;
} Index_Word_Collect_Traverse_Data;

static gboolean Index_WordHood_collect_traverse_func(gchar *word, gint score,
                                                     gpointer user_data){
    register VFSM_Int state;
    VFSM_Int leaf;
    register Index_Word_Collect_Traverse_Data *iwctd = user_data;
    state = VFSM_word2state(iwctd->index->vfsm, word);
    leaf = VFSM_state2leaf(iwctd->index->vfsm, state);
    iwctd->seed.leaf = leaf;
    if(iwctd->index_strand->word_table[leaf] >= 0)
        g_array_append_val(iwctd->word_seed_list, iwctd->seed);
    return FALSE;
    }
/* FIXME: optimisation: change to VFSM-based wordhoods */

static void Index_expand_word_seed_list(Index *index,
                                        Index_Strand *index_strand,
                                        WordHood *wh, GArray *word_seed_list){
    register gint i;
    register VFSM_Int state;
    /* FIXME: set threshold, use_dropoff elsewhere */
    /* FIXME: set word limit for dna/protein/codon and from default or client */
    register gchar *word = g_new(gchar, index->header->word_length+1);
    register gint init_seed_list_len = word_seed_list->len;
    register Index_WordSeed *init_seed_list = g_new(Index_WordSeed,
                                                    init_seed_list_len);
    register Index_WordSeed *orig_seed;
    for(i = 0; i < word_seed_list->len; i++){
        orig_seed = &g_array_index(word_seed_list, Index_WordSeed, i);
        init_seed_list[i].leaf = orig_seed->leaf;
        init_seed_list[i].query_pos = orig_seed->query_pos;
        }
    g_array_set_size(word_seed_list, 0);
    Index_Word_Collect_Traverse_Data iwctd;
    iwctd.index = index;
    iwctd.index_strand = index_strand;
    iwctd.word_seed_list = word_seed_list;
    /**/
    for(i = 0; i < init_seed_list_len; i++){
        state = VFSM_leaf2state(index->vfsm, init_seed_list[i].leaf);
        VFSM_state2word(index->vfsm, state, word);
        iwctd.seed.query_pos = init_seed_list[i].query_pos;
        WordHood_traverse(wh, Index_WordHood_collect_traverse_func,
                          word, index->header->word_length, &iwctd);
        }
    g_free(init_seed_list);
    g_free(word);
    return;
    }
/* FIXME: also need ability to catch redundancy over a whole set of queries ?
 *       - do this by doing Index vs Index matching
 *       - merge word lists
 *           - at each matching word pair
 *           - get query seed word
 *           - get target seed word
 *           - expand EITHER query OR target with its neighbourlist
 *           - for each query seed address
 *               - for each target seed address
 *                   - seed <query_id> <query_pos> <target_id> <target_pos>
 */
/* FIXME: skip storage requirement for intermediate results,
 *        by using PQ to reuse results.
 *        This can be done for both the single and joint algorithms.
 */

static HSPset *Index_get_HSPset(Index *index, gint target_id,
                                HSP_Param *hsp_param, Sequence *query,
                                HSPset_SList_Node *node,
                                gboolean revcomp_target){
    register Sequence *target, *tmp_seq;
    register HSPset *hsp_set;
    register HSPset_SList_Node *curr;
    register gint i;
    /* Check there are enough seeds */
    curr = node;
    for(i = 0; i < hsp_param->has->seed_repeat; i++){
        if(!curr)
            return NULL;
        curr = curr->next;
        }
    /* Make HSPset for current target */
    target = Dataset_get_sequence(index->dataset, target_id);
    /**/
    if(revcomp_target){
        tmp_seq = target;
        target = Sequence_revcomp(target);
        Sequence_destroy(tmp_seq);
        }
    hsp_set = HSPset_create(query, target, hsp_param);
    HSPset_seed_all_qy_sorted(hsp_set, node);
    Sequence_destroy(target);
    if(HSPset_is_empty(hsp_set)){
        HSPset_destroy(hsp_set);
        return NULL;
        }
    /* HSPset_print(hsp_set); */
    return hsp_set;
    }

/**/

static Index_HSPset *Index_HSPset_create(HSPset *hsp_set,
                                         gint target_id){
    register Index_HSPset *index_hsp_set = g_new(Index_HSPset, 1);
    index_hsp_set->hsp_set = HSPset_share(hsp_set);
    index_hsp_set->target_id = target_id;
    return index_hsp_set;
    }

void Index_HSPset_destroy(Index_HSPset *index_hsp_set){
    HSPset_destroy(index_hsp_set->hsp_set);
    g_free(index_hsp_set);
    return;
    }

/**/

static GPtrArray *Index_seed_HSPsets(Index *index, Index_Strand *index_strand,
                                     GArray *seed_list, HSP_Param *hsp_param,
                                     Sequence *query, gboolean revcomp_target,
                                     GArray *interval_list){
    register GPtrArray *hsp_set_list = g_ptr_array_new();
    register gint i, j, word_id, address_list_len;
    gint target_id;
    register Index_WordSeed *seed;
    register Index_Address *address_list;
    register Index_Word *index_word;
    register HSPset_SList_Node *node, **target_bin
        = g_new0(HSPset_SList_Node*, index->dataset->header->number_of_seqs);
    /* FIXME: allocate target_bin outside of this function
     *        (and clear after use) (per-thread)
     */
    register GArray *target_id_list = g_array_new(FALSE, FALSE, sizeof(gint));
    register HSPset *hsp_set;
    register Index_HSPset *index_hsp_set;
    register RecycleBin *hsp_slist_recycle = HSPset_SList_RecycleBin_create();
    /**/
    for(i = 0; i < seed_list->len; i++){
        seed = &g_array_index(seed_list, Index_WordSeed, i);
        word_id = index_strand->word_table[seed->leaf];
        index_word = &index_strand->word_list[word_id];
        g_assert(word_id >= 0);
        /* FIXME: optimisation: reuse address list to avoid allocation */
        address_list = g_new(Index_Address, index_word->freq_count);
        address_list_len = Index_Word_get_address_list(index, index_strand,
                                                   index_word, address_list,
                                                   interval_list);
        if(!address_list_len){
            g_free(address_list);
            continue; /* May be absent when using interval_list */
            }
        g_assert(index_word->freq_count);
        g_assert(address_list_len);
        for(j = 0; j < address_list_len; j++){
            target_id = address_list[j].sequence_id;
            if(!target_bin[target_id])
                g_array_append_val(target_id_list, target_id);
            target_bin[target_id]
                = HSPset_SList_append(hsp_slist_recycle, target_bin[target_id],
                       seed->query_pos, address_list[j].position);
            }
        g_free(address_list);
        }
    for(i = 0; i < target_id_list->len; i++){
        target_id = g_array_index(target_id_list, gint, i);
        node = target_bin[target_id];
        g_assert(node);
        hsp_set = Index_get_HSPset(index, target_id, hsp_param, query,
                                   node, revcomp_target);
        if(hsp_set){
            index_hsp_set = Index_HSPset_create(hsp_set, target_id);
            HSPset_destroy(hsp_set);
            /* g_message("Have [%d] seeds", hsp_set->hsp_list->len); */
            g_ptr_array_add(hsp_set_list, index_hsp_set);
            }
        }
    g_free(target_bin);
    g_array_free(target_id_list, TRUE);
    RecycleBin_destroy(hsp_slist_recycle);
    if(!hsp_set_list->len){
        g_ptr_array_free(hsp_set_list, TRUE);
        return NULL;
        }
    return hsp_set_list;
    }

static GArray *Index_get_word_seed_list(Index*index, Sequence *query,
                                        Index_Strand *index_strand,
                                        HSP_Param *hsp_param){
    register GArray *word_seed_list = g_array_new(FALSE, FALSE,
                                                  sizeof(Index_WordSeed));
    register gint i;
    register Sequence *aa_seq;
    if(hsp_param->match->query->is_translated){
        g_assert(hsp_param->match->mas->translate);
        for(i = 0; i < 3; i++){
            aa_seq = Sequence_translate(query,
                                        hsp_param->match->mas->translate, i+1);
            Index_get_query_word_list(index, index_strand,
                                      aa_seq, word_seed_list, i+1, hsp_param);
            Sequence_destroy(aa_seq);
            }
    } else {
        Index_get_query_word_list(index, index_strand,
                                  query, word_seed_list, 0, hsp_param);
        }
    /**/
    /* g_message("Using [%d] initial word seeds", word_seed_list->len); */
    if(hsp_param->wordhood){
        Index_expand_word_seed_list(index, index_strand,
                                    hsp_param->wordhood, word_seed_list);
        /*
        g_message("... expanded to [%d] words with wordhood",
                  word_seed_list->len);
        */
        }
    return word_seed_list;
    }

static Index_Strand *Index_get_index_strand(Index *index,
                                            gboolean revcomp_target){
    if(revcomp_target){
        if(!index->revcomp)
            g_error("No revcomp strand available in Index");
        return index->revcomp;
        }
    return index->forward;
    }

static GPtrArray *Index_get_HSPsets_interval(Index *index, HSP_Param *hsp_param,
                             Sequence *query, gboolean revcomp_target,
                             GArray *interval_list, GArray *word_seed_list){
    register GPtrArray *hsp_set_list;
    register Index_Strand *index_strand = Index_get_index_strand(index,
                                                                 revcomp_target);
    register gboolean free_word_seed_list = FALSE;
    /* Traverse the query using VSFM */
    g_assert(index_strand);
    if(!word_seed_list){
        word_seed_list = Index_get_word_seed_list(index, query,
                                                  index_strand, hsp_param);
        free_word_seed_list = TRUE;
        }
    hsp_set_list = Index_seed_HSPsets(index, index_strand,
                                      word_seed_list, hsp_param, query,
                                      revcomp_target, interval_list);
    if(free_word_seed_list)
        g_array_free(word_seed_list, TRUE);
    return hsp_set_list;
    }
/* FIXME: optimisation: detect repeated words with two query passes */

GPtrArray *Index_get_HSPsets(Index *index, HSP_Param *hsp_param,
                             Sequence *query, gboolean revcomp_target){
    return Index_get_HSPsets_interval(index, hsp_param, query,
                                      revcomp_target, NULL, NULL);
    }

/**/

typedef struct {
    gboolean  go_fwd;
    gboolean  go_rev;
         HSP *hsp;
} Index_Subseed;

typedef struct {
    RangeTree *keeper_hsp_tree;    /* accepted HSPs */
    RangeTree *cand_hsp_tree;      /* candidate HSPs */
          HSP *max_cobs_cand_hsp;  /* candidate HSP with largest cobs */
     NOI_Tree *coverage;           /* searched regions of target */
    GPtrArray *hspset_list;        /* contributing HSPsets */
    GPtrArray *subseed_list;
    GPtrArray *next_subseed_list;
         gint  target_id;
} Index_Geneseed;

static void Index_Subseed_add_to_subseed_list(GPtrArray *subseed_list, HSP *hsp,
                                              gboolean go_fwd, gboolean go_rev){
    register Index_Subseed *index_subseed = g_new(Index_Subseed, 1);
    index_subseed->hsp = hsp;
    index_subseed->go_fwd = go_fwd;
    index_subseed->go_rev = go_rev;
    g_ptr_array_add(subseed_list, index_subseed);
    return;
    }

static Index_Geneseed *Index_Geneseed_create(Index_HSPset *index_hspset,
                                             NOI_Tree_Set *nts){
    register Index_Geneseed *index_geneseed = g_new(Index_Geneseed, 1);
    register gint i;
    register HSP *hsp;
    index_geneseed->keeper_hsp_tree = RangeTree_create();
    index_geneseed->cand_hsp_tree = RangeTree_create();
    index_geneseed->max_cobs_cand_hsp = NULL;
    index_geneseed->coverage = NOI_Tree_create(nts);
    index_geneseed->hspset_list = g_ptr_array_new();
    index_geneseed->subseed_list = g_ptr_array_new();
    index_geneseed->next_subseed_list = g_ptr_array_new();
    index_geneseed->target_id = index_hspset->target_id;
    for(i = 0; i < index_hspset->hsp_set->hsp_list->len; i++){
        hsp = index_hspset->hsp_set->hsp_list->pdata[i];
        RangeTree_add(index_geneseed->keeper_hsp_tree,
                      HSP_query_cobs(hsp), HSP_target_cobs(hsp), hsp);
        Index_Subseed_add_to_subseed_list(index_geneseed->subseed_list, hsp,
                                          TRUE, TRUE);
        }
    g_ptr_array_add(index_geneseed->hspset_list,
                    HSPset_share(index_hspset->hsp_set));
    return index_geneseed;
    }

static void Index_Geneseed_destroy(Index_Geneseed *index_geneseed,
                                   NOI_Tree_Set *nts){
    register gint i;
    register HSPset *hspset;
    for(i = 0; i < index_geneseed->hspset_list->len; i++){
        hspset = index_geneseed->hspset_list->pdata[i];
        HSPset_destroy(hspset);
        }
    RangeTree_destroy(index_geneseed->keeper_hsp_tree, NULL, NULL);
    RangeTree_destroy(index_geneseed->cand_hsp_tree, NULL, NULL);
    NOI_Tree_destroy(index_geneseed->coverage, nts);
    g_free(index_geneseed);
    return;
    }

static gboolean Index_Geneseed_collect_hspset_Rangetree_traverse(gint x, gint y,
                gpointer info, gpointer user_data){
    register HSP *hsp = info;
    register HSPset *hsp_set = user_data;
    HSPset_add_known_hsp(hsp_set, hsp->query_start, hsp->target_start,
                         hsp->length);
    return FALSE;
    }

static Index_HSPset *Index_Geneseed_collect_hspset(
                     Index_Geneseed *index_geneseed){
    register Index_HSPset *index_hspset;
    register HSPset *hsp_set, *source_hspset;
    if(RangeTree_is_empty(index_geneseed->keeper_hsp_tree))
        return NULL;
    g_assert(index_geneseed->hspset_list->len);
    source_hspset = index_geneseed->hspset_list->pdata[0];
    hsp_set = HSPset_create(source_hspset->query, source_hspset->target,
                            source_hspset->param);
    RangeTree_traverse(index_geneseed->keeper_hsp_tree,
                       Index_Geneseed_collect_hspset_Rangetree_traverse,
                       hsp_set);
    g_assert(!HSPset_is_empty(hsp_set));
    HSPset_finalise(hsp_set);
    /* g_message("collected [%d] hsps", hsp_set->hsp_list->len); */
    index_hspset = Index_HSPset_create(hsp_set, index_geneseed->target_id);
    HSPset_destroy(hsp_set);
    return index_hspset;
    }

/**/

typedef struct {
       GPtrArray *index_geneseed_list;
    NOI_Tree_Set *nts;
    /**/
          GArray *curr_interval_list;
            gint  curr_target_id;
    /**/
            gint  max_query_span;
            gint  max_target_span;
} Index_Geneseed_List;

static void Index_Interval_NOI_Tree_traverse(gint start, gint length,
                                                      gpointer user_data){
    register Index_Geneseed_List *geneseed_list = user_data;
    Index_Interval interval;
    /**/
    interval.start = start;
    interval.length = length;
    interval.target_id = geneseed_list->curr_target_id;
    g_array_append_val(geneseed_list->curr_interval_list, interval);
    /*
    g_message("add interval [%d,%d,%d] to interval_list",
            start, length, geneseed_list->curr_target_id);
    */
    return;
    }
/* Need to have interval_list, interval, target_id */

static int Index_Geneseed_Subseed_compare(const void *a,
                                          const void *b){
    register Index_Subseed **subseed_a = (Index_Subseed**)a,
                           **subseed_b = (Index_Subseed**)b;
    return (*subseed_a)->hsp->target_start
         - (*subseed_b)->hsp->target_start;
    }

#if 0
static void Index_Geneseed_get_regions(Index_Geneseed *index_geneseed,
                                       NOI_Tree_Set *nts, gint max_target_span){
    register gint i, start, length;
    register gint target_range;
    register Index_Subseed *subseed;
    register HSP *hsp;
    NOI_Tree_delta_init(index_geneseed->coverage, nts);
    /* FIXME: optimisation replace qsort of HSPset with radix type sort */
    qsort(index_geneseed->subseed_list->pdata, index_geneseed->subseed_list->len,
          sizeof(gpointer), Index_Geneseed_Subseed_compare);
    for(i = 0; i < index_geneseed->subseed_list->len; i++){
        subseed = index_geneseed->subseed_list->pdata[i];
        hsp = subseed->hsp;
        target_range = max_target_span
                     + ((HSP_target_cobs(hsp)-hsp->target_start) * 2);
        if(subseed->go_rev){
            start = hsp->target_start - target_range;
            if(start < 0)
                start = 0;
            length = HSP_target_end(hsp) + target_range - start;
            if((start+length) >= hsp->hsp_set->target->len)
                length = hsp->hsp_set->target->len - start;
            NOI_Tree_insert(index_geneseed->coverage, nts, start, length);
            }
        if(subseed->go_fwd){
            start = HSP_target_cobs(hsp);
            length = start + target_range;
            if((start+length) >= hsp->hsp_set->target->len)
                length = hsp->hsp_set->target->len - start;
            NOI_Tree_insert(index_geneseed->coverage, nts, start, length);
            }
        }
    NOI_Tree_delta_traverse(index_geneseed->coverage, nts,
                            Index_Interval_NOI_Tree_traverse);
    return;
    }
#endif /* 0 */

static void Index_Geneseed_get_regions(Index_Geneseed *index_geneseed,
                                       NOI_Tree_Set *nts, gint max_target_span){
    register gint i, start, length;
    register gint target_range;
    register Index_Subseed *subseed;
    register HSP *hsp;
    NOI_Tree_delta_init(index_geneseed->coverage, nts);
    /* FIXME: optimisation replace qsort of HSPset with radix type sort */
    qsort(index_geneseed->subseed_list->pdata, index_geneseed->subseed_list->len,
          sizeof(gpointer), Index_Geneseed_Subseed_compare);
    for(i = 0; i < index_geneseed->subseed_list->len; i++){
        subseed = index_geneseed->subseed_list->pdata[i];
        hsp = subseed->hsp;
        target_range = max_target_span
                     + ((HSP_target_cobs(hsp)-hsp->target_start) * 2);
        if(subseed->go_rev){
            start = HSP_target_cobs(hsp) - target_range;
            if(start < 0)
                start = 0;
            length = HSP_target_cobs(hsp) - start;
            g_assert((start+length) <= hsp->hsp_set->target->len);
            NOI_Tree_insert(index_geneseed->coverage, nts, start, length);
            }
        if(subseed->go_fwd){
            start = HSP_target_cobs(hsp);
            length = start + target_range;
            if((start+length) >= hsp->hsp_set->target->len)
                length = hsp->hsp_set->target->len - start;
            NOI_Tree_insert(index_geneseed->coverage, nts, start, length);
            }
        }
    NOI_Tree_delta_traverse(index_geneseed->coverage, nts,
                            Index_Interval_NOI_Tree_traverse);
    return;
    }

/**/

static int Index_Geneseed_compare(const void *a,
                                  const void *b){
    register Index_Geneseed **geneseed_a = (Index_Geneseed**)a,
                            **geneseed_b = (Index_Geneseed**)b;
    return (*geneseed_a)->target_id
         - (*geneseed_b)->target_id;
    }

static Index_Geneseed_List *Index_Geneseed_List_create(GPtrArray *hspset_list,
                                     gint max_query_span, gint max_target_span){
    register Index_Geneseed_List *geneseed_list
     = g_new(Index_Geneseed_List, hspset_list->len);
    register Index_HSPset *index_hsp_set;
    register Index_Geneseed *index_geneseed;
    register gint i;
    geneseed_list->curr_interval_list = NULL;
    geneseed_list->curr_target_id = -1;
    /**/
    geneseed_list->max_query_span = max_query_span;
    geneseed_list->max_target_span = max_target_span;
    /**/
    geneseed_list->nts = NOI_Tree_Set_create(geneseed_list);
    geneseed_list->index_geneseed_list = g_ptr_array_new();
    for(i = 0; i < hspset_list->len; i++){
        index_hsp_set = hspset_list->pdata[i];
        index_geneseed = Index_Geneseed_create(index_hsp_set, geneseed_list->nts);
        g_ptr_array_add(geneseed_list->index_geneseed_list, index_geneseed);
        }
    g_assert(hspset_list->len == geneseed_list->index_geneseed_list->len);
    qsort(geneseed_list->index_geneseed_list->pdata,
          geneseed_list->index_geneseed_list->len,
          sizeof(gpointer), Index_Geneseed_compare);
    return geneseed_list;
    }

static void Index_Geneseed_List_destroy(Index_Geneseed_List *geneseed_list){
    register gint i;
    register Index_Geneseed *index_geneseed;
    for(i = 0; i < geneseed_list->index_geneseed_list->len; i++){
        index_geneseed = geneseed_list->index_geneseed_list->pdata[i];
        Index_Geneseed_destroy(index_geneseed, geneseed_list->nts);
        }
    g_ptr_array_free(geneseed_list->index_geneseed_list, TRUE);
    NOI_Tree_Set_destroy(geneseed_list->nts);
    g_free(geneseed_list);
    return;
    }

static GArray *Index_Geneseed_List_get_regions(
               Index_Geneseed_List *geneseed_list){
    register gint i;
    register Index_Geneseed *index_geneseed;
    register GArray *interval_list = g_array_new(FALSE, FALSE,
                                                 sizeof(Index_Interval));
    geneseed_list->curr_interval_list = interval_list;
    for(i = 0; i < geneseed_list->index_geneseed_list->len; i++){
        index_geneseed = geneseed_list->index_geneseed_list->pdata[i];
        geneseed_list->curr_target_id = index_geneseed->target_id;
        if(index_geneseed->subseed_list->len)
            Index_Geneseed_get_regions(index_geneseed, geneseed_list->nts,
                                       geneseed_list->max_target_span);
        }
    geneseed_list->curr_interval_list = NULL;
    geneseed_list->curr_target_id = -1;
    if(!interval_list->len){
        g_array_free(interval_list, TRUE);
        return NULL;
        }
    return interval_list;
    }

/**/

static GPtrArray *Index_Geneseed_List_get_subseeds(Index *index,
                                HSP_Param *hsp_param,
                                Sequence *query, gboolean revcomp_target,
                                GArray *interval_list, GArray *word_seed_list){
    return Index_get_HSPsets_interval(index, hsp_param, query, revcomp_target,
                             interval_list, word_seed_list);
    }

/**/

static gboolean Index_Geneseed_refine_RangeTree_report_fwd(gint x, gint y,
                                    gpointer info, gpointer user_data){
    /* Called for each selected new HSP */
    register HSP *hsp = info;
    register Index_Geneseed *index_geneseed = user_data;
    /* Add to keeper rangetree if not already present */
    if(!RangeTree_check_pos(index_geneseed->keeper_hsp_tree, x, y)){
        RangeTree_add(index_geneseed->keeper_hsp_tree, x, y, hsp);
        Index_Subseed_add_to_subseed_list(index_geneseed->next_subseed_list,
                                          hsp, TRUE, FALSE);
        }
    return FALSE;
    }

static gboolean Index_Geneseed_refine_RangeTree_report_rev(gint x, gint y,
                                    gpointer info, gpointer user_data){
    /* Called for each selected new HSP */
    register HSP *hsp = info;
    register Index_Geneseed *index_geneseed = user_data;
    /* Add to keeper rangetree if not already present */
    if(!RangeTree_check_pos(index_geneseed->keeper_hsp_tree, x, y)){
        RangeTree_add(index_geneseed->keeper_hsp_tree, x, y, hsp);
        Index_Subseed_add_to_subseed_list(index_geneseed->next_subseed_list,
                                          hsp, FALSE, TRUE);
        }
    return FALSE;
    }

static int Index_HSPset_compare(const void *a, const void *b){
    register Index_HSPset **index_hspset_a = (Index_HSPset**)a,
                          **index_hspset_b = (Index_HSPset**)b;
    return (*index_hspset_a)->target_id
         - (*index_hspset_b)->target_id;
    }

static void Index_Geneseed_refine_subseeds(
            Index_Geneseed_List *index_geneseed_list,
            GPtrArray *subseed_hsp_list){
    register gint i, j, ig_pos = 0, query_range, target_range;
    register Index_HSPset *index_hspset;
    register Index_Geneseed *index_geneseed;
    register HSP *hsp;
    register GPtrArray *swap_list;
    register Index_Subseed *subseed;
    /* Sort subseed_hsp_list by target_id */
    qsort(subseed_hsp_list->pdata, subseed_hsp_list->len,
          sizeof(gpointer), Index_HSPset_compare);
    /* For each subseed HSPset */
    for(i = 0; i < subseed_hsp_list->len; i++){
        index_hspset = subseed_hsp_list->pdata[i];
        /* Find the corresponding Index_Geneseed */
        index_geneseed = NULL;
        while(ig_pos < index_geneseed_list->index_geneseed_list->len){
            index_geneseed
                = index_geneseed_list->index_geneseed_list->pdata[ig_pos];
            if(index_hspset->target_id == index_geneseed->target_id)
                break;
            index_geneseed = NULL;
            ig_pos++;
            }
        g_assert(index_geneseed);
        if(!index_geneseed->subseed_list->len)
            continue;
        /* Put new HSPs into cand RangeTree */
        for(j = 0; j < index_hspset->hsp_set->hsp_list->len; j++){
            hsp = index_hspset->hsp_set->hsp_list->pdata[j];
            if(!RangeTree_check_pos(index_geneseed->cand_hsp_tree,
                              HSP_query_cobs(hsp), HSP_target_cobs(hsp))){
                RangeTree_add(index_geneseed->cand_hsp_tree,
                              HSP_query_cobs(hsp), HSP_target_cobs(hsp), hsp);
                }
            /* Store max cobs candidate hsp */
            if((!index_geneseed->max_cobs_cand_hsp)
               || (index_geneseed->max_cobs_cand_hsp->cobs < hsp->cobs))
                index_geneseed->max_cobs_cand_hsp = hsp;
            }
        /* For each old HSP */
        for(j = 0; j < index_geneseed->subseed_list->len; j++){
            subseed = index_geneseed->subseed_list->pdata[j];
            hsp = subseed->hsp;
            query_range = index_geneseed_list->max_query_span
                    + (((HSP_query_end(hsp) - HSP_query_cobs(hsp))
                      + (HSP_query_cobs(index_geneseed->max_cobs_cand_hsp)
                       - index_geneseed->max_cobs_cand_hsp->query_start)) * 2);
            target_range = index_geneseed_list->max_target_span
                    + (((HSP_target_end(hsp) - HSP_target_cobs(hsp))
                      + (HSP_target_cobs(index_geneseed->max_cobs_cand_hsp)
                       - index_geneseed->max_cobs_cand_hsp->target_start)) * 2);
            /* Select HSP within range forwards from cand RangeTree */
            if(subseed->go_fwd)
                RangeTree_find(index_geneseed->cand_hsp_tree,
                           HSP_query_cobs(hsp), query_range,
                           HSP_target_cobs(hsp), target_range,
                           Index_Geneseed_refine_RangeTree_report_fwd,
                           index_geneseed);
            /* Select HSP within range backwards from cand RangeTree */
            if(subseed->go_rev)
                RangeTree_find(index_geneseed->cand_hsp_tree,
                           HSP_query_cobs(hsp)-query_range, query_range,
                           HSP_target_cobs(hsp)-target_range, target_range,
                           Index_Geneseed_refine_RangeTree_report_rev,
                           index_geneseed);
            g_free(subseed);
            }
        /* Add new HSPset to hspset_list */
        g_ptr_array_add(index_geneseed->hspset_list,
                        HSPset_share(index_hspset->hsp_set));
        /* Empty subseed_list and swap with next_subseed_list */
        g_ptr_array_set_size(index_geneseed->subseed_list, 0);
        swap_list = index_geneseed->subseed_list;
        index_geneseed->subseed_list = index_geneseed->next_subseed_list;
        index_geneseed->next_subseed_list = swap_list;
        }
    return;
    }

static GPtrArray *Index_Geneseed_collect_hsps(
                  Index_Geneseed_List *index_geneseed_list){
    register GPtrArray *hsp_list = g_ptr_array_new();
    register gint i;
    register Index_Geneseed *index_geneseed;
    register Index_HSPset *index_hspset;
    for(i = 0; i < index_geneseed_list->index_geneseed_list->len; i++){
        index_geneseed = index_geneseed_list->index_geneseed_list->pdata[i];
        index_hspset = Index_Geneseed_collect_hspset(index_geneseed);
        if(index_hspset)
            g_ptr_array_add(hsp_list, index_hspset);
        }
    if(!hsp_list->len){
        g_ptr_array_free(hsp_list, TRUE);
        return NULL;
        }
    return hsp_list;
    }

static void Index_HSPset_List_destroy(GPtrArray *index_hspset_list){
    register gint i;
    register Index_HSPset *index_hsp_set;
    for(i = 0; i < index_hspset_list->len; i++){
        index_hsp_set = index_hspset_list->pdata[i];
        Index_HSPset_destroy(index_hsp_set);
        }
    g_ptr_array_free(index_hspset_list, TRUE);
    return;
    }

GPtrArray *Index_get_HSPsets_geneseed(Index *index, HSP_Param *hsp_param,
                                  Sequence *query, gboolean revcomp_target,
                                  gint geneseed_threshold, gint geneseed_repeat,
                                  gint max_query_span, gint max_target_span){
    register gint original_hsp_threshold = hsp_param->threshold,
                  original_seed_repeat = hsp_param->seed_repeat;
    register GPtrArray *geneseed_hsp_list, *subseed_hsp_list,
                       *final_hsp_list;
    register Index_Geneseed_List *index_geneseed_list;
    register GArray *interval_list;
    register Index_Strand *index_strand = Index_get_index_strand(index,
                                                                 revcomp_target);
    register GArray *word_seed_list = Index_get_word_seed_list(index, query,
                                                     index_strand, hsp_param);
    /* g_message("have [%d] words", word_seed_list->len); */
    /* Find geneseed HSPs */
    HSP_Param_set_hsp_threshold(hsp_param, geneseed_threshold);
    HSP_Param_set_seed_repeat(hsp_param, geneseed_repeat);
    geneseed_hsp_list = Index_get_HSPsets_interval(index, hsp_param, query,
                                                   revcomp_target, NULL,
                                                   word_seed_list);
    HSP_Param_set_hsp_threshold(hsp_param, original_hsp_threshold);
    HSP_Param_set_seed_repeat(hsp_param, original_seed_repeat);
    /**/
    if(!geneseed_hsp_list){ /* if no hsps */
        g_array_free(word_seed_list, TRUE);
        return NULL;
        }
    /* g_message("start with [%d] geneseed hspsets", geneseed_hsp_list->len); */
    index_geneseed_list = Index_Geneseed_List_create(geneseed_hsp_list,
                                       max_query_span, max_target_span);
    Index_HSPset_List_destroy(geneseed_hsp_list);
    /* Make geneseed list */
    do {
        interval_list = Index_Geneseed_List_get_regions(index_geneseed_list);
        if(!interval_list)
            break;
        subseed_hsp_list = Index_Geneseed_List_get_subseeds(index, hsp_param,
                                          query, revcomp_target,
                                          interval_list, word_seed_list);
        g_array_free(interval_list, TRUE);
        if(!subseed_hsp_list)
            break;
        Index_Geneseed_refine_subseeds(index_geneseed_list, subseed_hsp_list);
        Index_HSPset_List_destroy(subseed_hsp_list);
    } while(TRUE);
    /* Collect kept hsps to return */
    final_hsp_list = Index_Geneseed_collect_hsps(index_geneseed_list);
    Index_Geneseed_List_destroy(index_geneseed_list);
    g_array_free(word_seed_list, TRUE);
    return final_hsp_list;
    }

/**/

static void Index_preload_index_Strand(Index *index, Index_Strand *index_strand){
    if(index_strand->index_cache)
        return;
    Index_fseek(index->fp, index_strand->strand_offset, SEEK_SET);
    index_strand->index_cache = BitArray_read(index->fp,
                                  index_strand->header.total_index_length);
    return;
    }

void Index_preload_index(Index *index){
    g_message("Preloading index");
    if(index->forward)
         Index_preload_index_Strand(index, index->forward);
    if(index->revcomp)
         Index_preload_index_Strand(index, index->revcomp);
    return;
    }

/* seq vs index algorithm:
 *
 * for each word in query
 *     append(data[word]->qpos_list, qpos)
 *     for each neighbour of word
 *         append(data[neighbour]->qpos_list, qpos)
 * for each word with a qpos_list
 *    PQ_push(word->address_list, qpos_list)
 * while(PQ_pop) order: <tid>
 *    seed each query_pos against all addresss in current target
 *    if not at end of address_list
 *        push back onto address_list
 *    else
 *        free address_list
 *    if changing target
 *        HSPset_seed_all_hsps(seed_list)
 *
 * Summary: word,qpos_list:PQ_push(address_list,qpos_list):pop(tid order)
 */

/* index vs index algorithm:
 * for each query word
 *     if corresponding target word exists
 *         PQ_push word pair
 *     for each neighbour word
 *         if corresponding target word exists
 *             PQpush word_pair
 * while(PQpop) order: <qid,tid>
 *     collate all seeds of <qid,tid>
 *     if sequence pair <qid,tid>
 *         HSPset_seed_all_hsps(seed_list)
 * Summary: PQ_push(word_pairs),PQ_pop()
 */

/**/

/* FIXME: optimisation: rice coding ??
 */

