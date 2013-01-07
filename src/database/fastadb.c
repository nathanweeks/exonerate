/****************************************************************\
*                                                                *
*  Library for manipulation of FASTA format databases            *
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

#include <stdio.h> /* For BUFSIZ */
#include <errno.h>
#include <string.h> /* For strerror() */
#include <ctype.h>  /* For isspace(),isprint() */

#include <sys/types.h> /* For stat() */
#include <sys/stat.h>  /* For stat() */
#include <unistd.h>    /* For stat() */
#include <dirent.h>    /* For readdir() */

#include "fastadb.h"

#ifndef ALPHABET_SIZE
#define ALPHABET_SIZE (sizeof(gchar)<<CHAR_BIT)
#endif /* ALPHABET_SIZE */

#define FASTADB_OUT_BUFFER_CHUNK_SIZE (BUFSIZ * 64)

FastaDB_ArgumentSet *FastaDB_ArgumentSet_create(Argument *arg){
    register ArgumentSet *as;
    static FastaDB_ArgumentSet fas = {NULL};
    if(arg){
        as = ArgumentSet_create("Fasta Database Options");
        ArgumentSet_add_option(as, '\0', "fastasuffix", "suffix",
           "Fasta file suffix filter (in subdirectories)", ".fa",
           Argument_parse_string, &fas.suffix_filter);
        Argument_absorb_ArgumentSet(arg, as);
        }
    return &fas;
    }

static gboolean FastaDB_file_is_directory(gchar *path){
    struct stat buf;
    if(stat(path, &buf))
        g_error("Cannot read file [%s]", path);
    return S_ISDIR(buf.st_mode);
    }

static GPtrArray *FastaDB_get_dir_contents(gchar *path){
    register GPtrArray *list = g_ptr_array_new();
    register DIR *dir = opendir(path);
    register struct dirent *entry;
    if(!dir)
        g_error("Could not open directory [%s]", path);
    readdir(dir); /* Skip '.' */
    readdir(dir); /* Skip '..' */
    while((entry = readdir(dir)))
        g_ptr_array_add(list, g_strdup(entry->d_name));
    closedir(dir);
    return list;
    }

static void FastaDB_expand_path_list_recur(GPtrArray *path_list,
               gchar *path, gchar *suffix_filter, gboolean is_subdir){
    register GPtrArray *dir_content;
    register gchar *name, *full_path;
    register gint i;
    if(FastaDB_file_is_directory(path)){
        dir_content = FastaDB_get_dir_contents(path);
        for(i = 0; i < dir_content->len; i++){
            name = dir_content->pdata[i];
            full_path = g_strdup_printf("%s%c%s",
                                        path, G_DIR_SEPARATOR, name);
            FastaDB_expand_path_list_recur(path_list, full_path,
                                           suffix_filter, TRUE);
            g_free(full_path);
            g_free(name);
            }
        g_ptr_array_free(dir_content, TRUE);
    } else {
        if(is_subdir && suffix_filter){
            for(i = strlen(path)-1; i >= 0; i--)
                if(path[i] == '.')
                    break;
            if((!i) || strcmp(path+i+1, suffix_filter))
                return;
            }
        g_ptr_array_add(path_list, g_strdup(path));
        }
    return;
    }

static GPtrArray *FastaDB_expand_path_list(GPtrArray *path_list,
                                           gchar *suffix_filter){
    register GPtrArray *expanded_path_list = g_ptr_array_new();
    register gchar *path;
    register gint i;
    if(suffix_filter && (suffix_filter[0] == '.'))
        suffix_filter++; /* Skip any leading dot on suffix_filter */
    for(i = 0; i < path_list->len; i++){
        path = path_list->pdata[i];
        FastaDB_expand_path_list_recur(expanded_path_list,
                                       path, suffix_filter, FALSE);
        }
    if(!expanded_path_list->len){
        for(i = 0; i < path_list->len; i++)
            g_message("Empty path: [%s]", (gchar*)path_list->pdata[i]);
        g_error("No valid file paths found");
        }
    return expanded_path_list;
    }

FastaDB *FastaDB_open_list(GPtrArray *path_list, Alphabet *alphabet){
    register FastaDB *fdb = g_new(FastaDB, 1);
    register FastaDB_ArgumentSet *fas
           = FastaDB_ArgumentSet_create(NULL);
    register gint i;
    register gchar *path;
    register GPtrArray *full_path_list;
    g_assert(path_list);
    fdb->ref_count = 1;
    if(alphabet)
        fdb->alphabet = Alphabet_share(alphabet);
    else
        fdb->alphabet = Alphabet_create(Alphabet_Type_UNKNOWN, FALSE);
    full_path_list = FastaDB_expand_path_list(path_list, fas->suffix_filter);
    fdb->cf = CompoundFile_create(full_path_list, TRUE);
    for(i = 0; i < full_path_list->len; i++){
        path = full_path_list->pdata[i];
        g_free(path);
        }
    g_ptr_array_free(full_path_list, TRUE);
    FastaDB_rewind(fdb);
    fdb->out_buffer_alloc = FASTADB_OUT_BUFFER_CHUNK_SIZE;
    fdb->out_buffer = g_malloc(sizeof(gchar)*fdb->out_buffer_alloc);
    fdb->line_length = -1;
    return fdb;
    }

FastaDB *FastaDB_open_list_with_limit(GPtrArray *path_list,
             Alphabet *alphabet, gint chunk_id, gint chunk_total){
    register FastaDB *fdb = FastaDB_open_list(path_list, alphabet);
    register CompoundFile_Pos start, stop, total_length, chunk_size;
    register CompoundFile_Location *start_cfl, *stop_cfl;
    if(chunk_total){
        g_assert(chunk_id);
        if(chunk_total < 1)
            g_error("Chunk total [%d] is too small", chunk_total);
        if((chunk_id < 1) || (chunk_id > chunk_total))
            g_error("Chunk id should be between 1 and %d", chunk_total);
        total_length = CompoundFile_get_length(fdb->cf);
        chunk_size = total_length / chunk_total;
        start = (chunk_id-1) * chunk_size;
        start = FastaDB_find_next_start(fdb, start);
        if(chunk_id == chunk_total){
            stop = total_length - 1;
        } else {
            stop = chunk_id * chunk_size;
            stop = FastaDB_find_next_start(fdb, stop);
            }
        start_cfl = CompoundFile_Location_from_pos(fdb->cf, start);
        stop_cfl = CompoundFile_Location_from_pos(fdb->cf, stop);
        CompoundFile_set_limits(fdb->cf, start_cfl, stop_cfl);
        CompoundFile_Location_destroy(start_cfl);
        CompoundFile_Location_destroy(stop_cfl);
        FastaDB_rewind(fdb);
        }
    return fdb;
    }

FastaDB *FastaDB_open(gchar *path, Alphabet *alphabet){
    register FastaDB *fdb;
    register GPtrArray *path_list = g_ptr_array_new();
    g_ptr_array_add(path_list, path);
    fdb = FastaDB_open_list(path_list, alphabet);
    g_ptr_array_free(path_list, TRUE);
    return fdb;
    }

#define FastaDB_putc(fdb, ch) \
    if(fdb->out_buffer_pos == fdb->out_buffer_alloc) \
        FastaDB_putc_realloc(fdb); \
    fdb->out_buffer[fdb->out_buffer_pos++] = (ch)

static void FastaDB_putc_realloc(FastaDB *fdb){
    fdb->out_buffer_alloc += FASTADB_OUT_BUFFER_CHUNK_SIZE;
    fdb->out_buffer = g_realloc(fdb->out_buffer,
                                fdb->out_buffer_alloc);
    return;
    }

/**/

FastaDB *FastaDB_share(FastaDB *fdb){
    g_assert(fdb);
    fdb->ref_count++;
    return fdb;
    }

FastaDB *FastaDB_dup(FastaDB *fdb){
    register FastaDB *nfdb = g_new(FastaDB, 1);
    nfdb->ref_count = 1;
    nfdb->alphabet = Alphabet_create(fdb->alphabet->type,
                                     fdb->alphabet->is_soft_masked);
    nfdb->cf = CompoundFile_dup(fdb->cf);
    nfdb->out_buffer_alloc = FASTADB_OUT_BUFFER_CHUNK_SIZE;
    nfdb->out_buffer = g_malloc(sizeof(gchar)*nfdb->out_buffer_alloc);
    nfdb->line_length = -1;
    FastaDB_rewind(nfdb);
    return nfdb;
    }

void FastaDB_close(FastaDB *fdb){
    g_assert(fdb);
    if(--fdb->ref_count)
        return;
    CompoundFile_destroy(fdb->cf);
    Alphabet_destroy(fdb->alphabet);
    g_free(fdb->out_buffer);
    g_free(fdb);
    return;
    }

void FastaDB_rewind(FastaDB *fdb){
    register gint ch, prev = '\n';
    g_assert(fdb);
    CompoundFile_rewind(fdb->cf);
    while((ch = CompoundFile_getc(fdb->cf)) != EOF){
        if((ch == '>') && (prev == '\n'))
            break;
        prev = ch;
        }
    return;
    }

CompoundFile_Pos FastaDB_find_next_start(FastaDB *fdb,
                                         CompoundFile_Pos pos){
    register gint ch, prev = '\n';
    g_assert(fdb->cf->element_list->len == 1);
    CompoundFile_seek(fdb->cf, pos);
    while((ch = CompoundFile_getc(fdb->cf)) != EOF){
        if((ch == '>') && (prev == '\n'))
            break;
        prev = ch;
        }
    return CompoundFile_ftell(fdb->cf)-1;
    }

gboolean FastaDB_file_is_fasta(gchar *path){
    register FILE *fp = fopen(path, "r");
    register gint ch;
    register gboolean result = FALSE;
    if(!fp)
        g_error("Could not open file [%s]", path);
    while(((ch = getc(fp)) != EOF)){
        if(isspace(ch))
            continue;
        else
            if(ch == '>'){ /* Assume is fasta format */
                result = TRUE;
                break;
                }
        else
            break;
        }
    fclose(fp);
    return result;
    }
/* Assumes a file is fasta format
 * if the first non-whitespace character is '>'
 */

gboolean FastaDB_is_finished(FastaDB *fdb){
    g_assert(fdb);
    return CompoundFile_is_finished(fdb->cf);
    }

void FastaDB_traverse(FastaDB *fdb, FastaDB_Mask mask,
                      FastaDB_TraverseFunc fdtf, gpointer user_data){
    register FastaDB_Seq *fdbs;
    g_assert(fdb);
    g_assert(fdtf);
    while((fdbs = FastaDB_next(fdb, mask))){
        if(fdtf(fdbs, user_data)){
            FastaDB_Seq_destroy(fdbs);
            break;
            }
        FastaDB_Seq_destroy(fdbs);
        }
    return;
    }

gsize FastaDB_memory_usage(FastaDB *fdb){
    return sizeof(FastaDB)
         + sizeof(Alphabet)
         + sizeof(CompoundFile)
         + (sizeof(gchar)*fdb->out_buffer_alloc);
    }

/**/

static FastaDB_Seq *FastaDB_Seq_create(FastaDB *fdb, Sequence *seq,
                                       CompoundFile_Location *cfl){
    register FastaDB_Seq *fdbs = g_new0(FastaDB_Seq, 1);
    fdbs->ref_count = 1;
    fdbs->source = FastaDB_share(fdb);
    fdbs->location = CompoundFile_Location_share(cfl);
    fdbs->seq = Sequence_share(seq);
    return fdbs;
    }

FastaDB_Seq *FastaDB_Seq_share(FastaDB_Seq *fdbs){
    g_assert(fdbs);
    fdbs->ref_count++;
    return fdbs;
    }

void FastaDB_Seq_destroy(FastaDB_Seq *fdbs){
    g_assert(fdbs);
    if(--fdbs->ref_count)
        return;
    CompoundFile_Location_destroy(fdbs->location);
    FastaDB_close(fdbs->source);
    Sequence_destroy(fdbs->seq);
    g_free(fdbs);
    return;
    }

FastaDB_Seq *FastaDB_Seq_revcomp(FastaDB_Seq *fdbs){
    register FastaDB_Seq *revcomp_fdbs;
    g_assert(fdbs);
    revcomp_fdbs = FastaDB_Seq_create(fdbs->source,
                                      fdbs->seq, fdbs->location);
    Sequence_destroy(revcomp_fdbs->seq);
    revcomp_fdbs->seq = Sequence_revcomp(fdbs->seq);
    return revcomp_fdbs;
    }

/**/

FastaDB_Seq *FastaDB_next(FastaDB *fdb, FastaDB_Mask mask){
    register FastaDB_Seq *fdbs;
    register gint ch;
    register gint id_pos = -1, def_pos = -1, seq_pos = -1,
                  seq_len = -1;
    register gint prev_len = -1, curr_length, line_start_seq_pos;
    register Sequence *seq;
    register CompoundFile_Location *location
           = CompoundFile_Location_tell(fdb->cf);
    fdb->out_buffer_pos = 0; /* Clear output buffer */
    if(mask & FastaDB_Mask_ID){ /* Record ID */
        id_pos = 0;
        while((ch = CompoundFile_getc(fdb->cf)) != EOF){
            /* if((ch == '\n') || (ch == ' ') || (ch == '\t')) */
            if(isspace(ch))
                break;
            FastaDB_putc(fdb, ch);
            }
        FastaDB_putc(fdb, '\0');
    } else { /* Skip ID */
        while((ch = CompoundFile_getc(fdb->cf)) != EOF){
            /* if((ch == '\n') || (ch == ' ') || (ch == '\t')) */
            if(isspace(ch))
                break;
            }
        }
    if(ch == EOF){
        CompoundFile_Location_destroy(location);
        return NULL;
        }
    if(ch != '\n'){ /* If DEF is present */
        if(mask & FastaDB_Mask_DEF){ /* Record DEF */
            def_pos = fdb->out_buffer_pos;
            while((ch = CompoundFile_getc(fdb->cf)) != EOF){
                if(ch == '\n')
                    break;
                FastaDB_putc(fdb, ch);
                }
            FastaDB_putc(fdb, '\0');
        } else { /* Skip DEF */
            while((ch = CompoundFile_getc(fdb->cf)) != EOF){
                if(ch == '\n')
                    break;
                }
            }
        }
    if(ch == EOF){
        CompoundFile_Location_destroy(location);
        return NULL;
        }
    if(mask & FastaDB_Mask_SEQ){ /* Record SEQ and LEN */
        seq_pos = fdb->out_buffer_pos;
        line_start_seq_pos = seq_pos;
        while(((ch = CompoundFile_getc(fdb->cf)) != EOF)
            && (ch != '>')){
            if(Alphabet_symbol_is_valid(fdb->alphabet, ch)){
                FastaDB_putc(fdb, ch);
            } else if(isspace(ch)){
                if(ch == '\n'){
                    curr_length = fdb->out_buffer_pos
                                - line_start_seq_pos;
                    line_start_seq_pos = fdb->out_buffer_pos;
                    if(prev_len != -1){
                        if(fdb->line_length == -1){
                            fdb->line_length = prev_len;
                        } else {
                            if(prev_len != fdb->line_length)
                                fdb->line_length = 0;
                            }
                        }
                    prev_len = curr_length;
                    /**/
                    }
            } else {
                g_error("Unrecognised symbol \'%c\' (ascii:%d)"
                        " file:[%s] seq:[%s] pos:[%d]",
                     isprint(ch)?ch:' ', ch,
                     CompoundFile_current_path(fdb->cf),
                     (id_pos == -1)?"unknown":&fdb->out_buffer[id_pos],
                     (fdb->out_buffer_pos - seq_pos - 1));
                }
            }
        FastaDB_putc(fdb, '\0');
        seq_len = fdb->out_buffer_pos - seq_pos - 1;
        if(fdb->line_length >= 0)
            if(fdb->line_length < prev_len)
                fdb->line_length = 0;
    } else if(mask & FastaDB_Mask_LEN){ /* Skip SEQ and record LEN */
        seq_len = 0;
        while(((ch = CompoundFile_getc(fdb->cf)) != EOF)
            && (ch != '>')){
            if( ((ch >= 'A') && (ch <= 'Z'))
             || ((ch >= 'a') && (ch <= 'z'))){
                seq_len++;
                }
           }
    } else { /* Just skip SEQ */
        while((ch = CompoundFile_getc(fdb->cf)) != EOF){
            if(ch == '>')
                break;
            }
        }
    seq = Sequence_create(
            ( id_pos == -1)?NULL:&fdb->out_buffer[id_pos],
            (def_pos == -1)?NULL:&fdb->out_buffer[def_pos],
            (seq_pos == -1)?NULL:&fdb->out_buffer[seq_pos],
            (seq_len == -1)?0:seq_len,
            (fdb->alphabet->type == Alphabet_Type_DNA)
                ? Sequence_Strand_FORWARD
                : Sequence_Strand_UNKNOWN,
            fdb->alphabet);
    fdbs = FastaDB_Seq_create(fdb, seq, location);
    CompoundFile_Location_destroy(location);
    Sequence_destroy(seq);
    if((mask & FastaDB_Mask_LEN) &&(!seq->len))
        g_warning("Warning zero length sequence [%s]", seq->id);
    return fdbs;
    }

/**/

FastaDB_Key *FastaDB_Key_create(FastaDB *source,
                                CompoundFile_Location *location,
                                Sequence_Strand strand,
                                gint seq_offset, gint length){
    register FastaDB_Key *fdbk = g_new(FastaDB_Key, 1);
    g_assert(source);
    fdbk->source = FastaDB_share(source);
    fdbk->location = CompoundFile_Location_share(location);
    fdbk->strand = strand;
    fdbk->seq_offset = seq_offset;
    fdbk->length = length;
    return fdbk;
    }

void FastaDB_Key_destroy(FastaDB_Key *fdbk){
    g_assert(fdbk);
    FastaDB_close(fdbk->source);
    CompoundFile_Location_destroy(fdbk->location);
    g_free(fdbk);
    return;
    }

FastaDB_Key *FastaDB_Seq_get_key(FastaDB_Seq *fdbs){
    register FastaDB_Key *fdbk;
    register gint seq_offset = strlen(fdbs->seq->id) + 2
                             + (fdbs->seq->def ? (strlen(fdbs->seq->def)+1) : 0);
    g_assert(fdbs);
    fdbk = FastaDB_Key_create(fdbs->source, fdbs->location,
                              fdbs->seq->strand, seq_offset,
                              fdbs->seq->len);
    return fdbk;
    }

FastaDB_Seq *FastaDB_fetch(FastaDB *fdb, FastaDB_Mask mask,
                           CompoundFile_Pos pos){
    register FastaDB_Seq *fdbs = NULL;
    register CompoundFile_Pos orig_pos = CompoundFile_ftell(fdb->cf);
    g_assert(fdb);
    g_assert(fdb->cf->element_list->len == 1);
    CompoundFile_seek(fdb->cf, pos);
    fdbs = FastaDB_next(fdb, mask);
    CompoundFile_seek(fdb->cf, orig_pos);
    return fdbs;
    }

FastaDB_Seq *FastaDB_Key_get_seq(FastaDB_Key *fdbk, FastaDB_Mask mask){
    register FastaDB_Seq *fdbs;
    register CompoundFile_Location *orig;
    register Sequence *revcomp_sequence;
    g_assert(fdbk);
    orig = CompoundFile_Location_tell(fdbk->source->cf);
    CompoundFile_Location_seek(fdbk->location);
    fdbs = FastaDB_next(fdbk->source, mask);
    g_assert(fdbs);
    CompoundFile_Location_seek(orig);
    CompoundFile_Location_destroy(orig);
    if(fdbk->strand == Sequence_Strand_REVCOMP){
        g_assert(fdbs->seq->strand == Sequence_Strand_FORWARD);
        fdbs->seq->strand = Sequence_Strand_REVCOMP;
        revcomp_sequence = Sequence_revcomp(fdbs->seq);
        Sequence_destroy(fdbs->seq);
        fdbs->seq = revcomp_sequence;
        }
    return fdbs;
    }

gchar *FastaDB_Key_get_def(FastaDB_Key *fdbk){
    register gint ch;
    register gchar *def = NULL;
    register GString *s;
#ifdef USE_PTHREADS
    pthread_mutex_lock(&fdbk->source->cf->compoundfile_mutex);
#endif /* USE_PTHREADS */
    CompoundFile_Location_seek(fdbk->location);
    while((ch = CompoundFile_getc(fdbk->source->cf)) != EOF)
        if((ch == '\n') || (ch == ' ') || (ch == '\t'))
            break;
    if(ch == '\n'){
#ifdef USE_PTHREADS
        pthread_mutex_unlock(&fdbk->source->cf->compoundfile_mutex);
#endif /* USE_PTHREADS */
        return NULL;
        }
    s = g_string_sized_new(64);
    while((ch = CompoundFile_getc(fdbk->source->cf)) != EOF){
        if(ch == '\n')
            break;
        s = g_string_append_c(s, ch);
        }
#ifdef USE_PTHREADS
    pthread_mutex_unlock(&fdbk->source->cf->compoundfile_mutex);
#endif /* USE_PTHREADS */
    def = s->str;
    g_string_free(s, FALSE);
    return def;
    }

/**/

static gpointer FastaDB_SparseCache_get_func_8_BIT(gint pos, gpointer page_data,
                                             gpointer user_data){
    return GINT_TO_POINTER((gint)((gchar*)page_data)[pos]);
    }

static gint FastaDB_SparseCache_copy_func_8_BIT(gint start, gint length,
                                                gchar *dst, gpointer page_data,
                                                gpointer user_data){
    register gint i;
    register gchar *str = page_data;
    for(i = 0; i < length; i++)
        dst[i] = str[start+i];
    return length;
    }

static gpointer FastaDB_SparseCache_get_func_0_BIT_LC(gint pos,
                                                      gpointer page_data,
                                                      gpointer user_data){
    return GINT_TO_POINTER((gint)'n');
    }

static gint FastaDB_SparseCache_copy_func_0_BIT_LC(gint start, gint length,
                                                   gchar *dst, gpointer page_data,
                                                   gpointer user_data){
    register gint i;
    for(i = 0; i < length; i++)
        dst[i] = 'n';
    return length;
    }

static gpointer FastaDB_SparseCache_get_func_0_BIT_UC(gint pos,
                                                      gpointer page_data,
                                                      gpointer user_data){
    return GINT_TO_POINTER((gint)'N');
    }

static gint FastaDB_SparseCache_copy_func_0_BIT_UC(gint start, gint length,
                                                   gchar *dst, gpointer page_data,
                                                   gpointer user_data){
    register gint i;
    for(i = 0; i < length; i++)
        dst[i] = 'N';
    return length;
    }

typedef enum {
    FastaDB_SparseCache_Mask_0_BIT_LC, /* [n] */
    FastaDB_SparseCache_Mask_0_BIT_UC, /* [N] */
    FastaDB_SparseCache_Mask_2_BIT_LC, /* [acgt] */
    FastaDB_SparseCache_Mask_2_BIT_UC, /* [ACGT] */
    FastaDB_SparseCache_Mask_4_BIT,    /* [ACGTNacgtn] */
    FastaDB_SparseCache_Mask_8_BIT
} FastaDB_SparseCache_Mask;

static void FastaDB_SparseCache_fill_mask(FastaDB_SparseCache_Mask *mask_index){
    register gint i;
    for(i = 0; i < ALPHABET_SIZE; i++)
        mask_index[i] = (1 << FastaDB_SparseCache_Mask_8_BIT);
    mask_index['N'] = (1 << FastaDB_SparseCache_Mask_0_BIT_UC);
    mask_index['n'] = (1 << FastaDB_SparseCache_Mask_0_BIT_LC);
    mask_index['A'] =
    mask_index['C'] =
    mask_index['G'] =
    mask_index['T'] = (1 << FastaDB_SparseCache_Mask_2_BIT_UC);
    mask_index['a'] =
    mask_index['c'] =
    mask_index['g'] =
    mask_index['t'] = (1 << FastaDB_SparseCache_Mask_2_BIT_LC);
    return;
    }

/**/

static void FastaDB_SparseCache_fill_index(gint *index, gchar *alphabet){
    register gint i;
    for(i = 0; i < ALPHABET_SIZE; i++)
        index[i] = -1;
    for(i = 0; alphabet[i]; i++)
        index[(guchar)alphabet[i]] = i;
    return;
    }

static gpointer FastaDB_SparseCache_get_func_4_BIT(gint pos,
                                                   gpointer page_data,
                                                   gpointer user_data){
    register gint ch = ((gchar*)page_data)[(pos >> 1)],
                  num = (ch >> ((pos & 1) << 2)) & 15;
    static gchar *alphabet = "ACGTNacgtn------";
    return GINT_TO_POINTER((gint)alphabet[num]);
    }

static gint FastaDB_SparseCache_copy_func_4_BIT(gint start, gint length,
                                                gchar *dst, gpointer page_data,
                                                gpointer user_data){
#ifndef G_DISABLE_ASSERT
    register gint i;
#endif /* G_DISABLE_ASSERT */
    register gint pos = 0, end = length;
    register gchar *str = page_data, *word;
    static gchar *table[256] = {
        "AA", "CA", "GA", "TA", "NA", "aA", "cA", "gA",
        "tA", "nA", "-A", "-A", "-A", "-A", "-A", "-A",
        "AC", "CC", "GC", "TC", "NC", "aC", "cC", "gC",
        "tC", "nC", "-C", "-C", "-C", "-C", "-C", "-C",
        "AG", "CG", "GG", "TG", "NG", "aG", "cG", "gG",
        "tG", "nG", "-G", "-G", "-G", "-G", "-G", "-G",
        "AT", "CT", "GT", "TT", "NT", "aT", "cT", "gT",
        "tT", "nT", "-T", "-T", "-T", "-T", "-T", "-T",
        "AN", "CN", "GN", "TN", "NN", "aN", "cN", "gN",
        "tN", "nN", "-N", "-N", "-N", "-N", "-N", "-N",
        "Aa", "Ca", "Ga", "Ta", "Na", "aa", "ca", "ga",
        "ta", "na", "-a", "-a", "-a", "-a", "-a", "-a",
        "Ac", "Cc", "Gc", "Tc", "Nc", "ac", "cc", "gc",
        "tc", "nc", "-c", "-c", "-c", "-c", "-c", "-c",
        "Ag", "Cg", "Gg", "Tg", "Ng", "ag", "cg", "gg",
        "tg", "ng", "-g", "-g", "-g", "-g", "-g", "-g",
        "At", "Ct", "Gt", "Tt", "Nt", "at", "ct", "gt",
        "tt", "nt", "-t", "-t", "-t", "-t", "-t", "-t",
        "An", "Cn", "Gn", "Tn", "Nn", "an", "cn", "gn",
        "tn", "nn", "-n", "-n", "-n", "-n", "-n", "-n",
        "A-", "C-", "G-", "T-", "N-", "a-", "c-", "g-",
        "t-", "n-", "--", "--", "--", "--", "--", "--",
        "A-", "C-", "G-", "T-", "N-", "a-", "c-", "g-",
        "t-", "n-", "--", "--", "--", "--", "--", "--",
        "A-", "C-", "G-", "T-", "N-", "a-", "c-", "g-",
        "t-", "n-", "--", "--", "--", "--", "--", "--",
        "A-", "C-", "G-", "T-", "N-", "a-", "c-", "g-",
        "t-", "n-", "--", "--", "--", "--", "--", "--",
        "A-", "C-", "G-", "T-", "N-", "a-", "c-", "g-",
        "t-", "n-", "--", "--", "--", "--", "--", "--",
        "A-", "C-", "G-", "T-", "N-", "a-", "c-", "g-",
        "t-", "n-", "--", "--", "--", "--", "--", "--"
        };
    g_assert(length);
    if(start&1){ /* If odd start, get first base */
        dst[0] = GPOINTER_TO_INT(FastaDB_SparseCache_get_func_4_BIT(start,
                                                     page_data, user_data));
        pos = 1;
        }
    if((length > 2) && ((start+length)&1)){ /* If odd end, get last base */
        dst[length-1] = GPOINTER_TO_INT(FastaDB_SparseCache_get_func_4_BIT(
                                       start+length-1, page_data, user_data));
        end--;
        }
    while(pos < end){
        word = table[(guchar)str[(start+pos)>>1]];
        g_assert(word[0] != '-');
        g_assert(word[1] != '-');
        dst[pos++] = word[0];
        dst[pos++] = word[1];
        }
#ifndef G_DISABLE_ASSERT
    for(i = 0; i < length; i++){
        g_assert(dst[i]
                 == GPOINTER_TO_INT(FastaDB_SparseCache_get_func_4_BIT(
                                    start+i, page_data, user_data)));
        }
#endif /* G_DISABLE_ASSERT */
    return length;
    }

static void FastaDB_SparseCache_compress_4bit(SparseCache_Page *page,
                                              gint len){
    register gint i;
    register gchar *seq =  page->data;
    static gint index[ALPHABET_SIZE];
    static gboolean index_is_filled = FALSE;
    if(!index_is_filled){
        FastaDB_SparseCache_fill_index(index, "ACGTNacgtn");
        index_is_filled = TRUE;
        }
    page->data = g_new0(gchar, (len >> 1)+1);
    page->data_size = sizeof(gchar)*((len >> 1)+1);
    for(i = 0; i < len; i++)
        ((gchar*)page->data)[(i >> 1)]
        |= (index[(guchar)seq[i]] << ((i&1) << 2));
    g_free(seq);
    page->get_func = FastaDB_SparseCache_get_func_4_BIT;
    page->copy_func = FastaDB_SparseCache_copy_func_4_BIT;
    return;
    }

/**/

static gpointer FastaDB_SparseCache_get_func_2_BIT_LC(gint pos,
                                                      gpointer page_data,
                                                      gpointer user_data){
    register gint ch = ((gchar*)page_data)[(pos >> 2)],
                  num = (ch >> ((pos & 3) << 1)) & 3;
    static gchar *alphabet = "acgt";
    return GINT_TO_POINTER((gint)alphabet[num]);
    }

static gint FastaDB_SparseCache_copy_func_2_BIT_LC(gint start, gint length,
                                                   gchar *dst, gpointer page_data,
                                                   gpointer user_data){
#ifndef G_DISABLE_ASSERT
    register gint i;
#endif /* G_DISABLE_ASSERT */
    register gint pos = 0, end = length;
    register gchar *str = page_data, *word;
    gchar *table[256] = {
        "aaaa", "caaa", "gaaa", "taaa", "acaa", "ccaa", "gcaa", "tcaa",
        "agaa", "cgaa", "ggaa", "tgaa", "ataa", "ctaa", "gtaa", "ttaa",
        "aaca", "caca", "gaca", "taca", "acca", "ccca", "gcca", "tcca",
        "agca", "cgca", "ggca", "tgca", "atca", "ctca", "gtca", "ttca",
        "aaga", "caga", "gaga", "taga", "acga", "ccga", "gcga", "tcga",
        "agga", "cgga", "ggga", "tgga", "atga", "ctga", "gtga", "ttga",
        "aata", "cata", "gata", "tata", "acta", "ccta", "gcta", "tcta",
        "agta", "cgta", "ggta", "tgta", "atta", "ctta", "gtta", "ttta",
        "aaac", "caac", "gaac", "taac", "acac", "ccac", "gcac", "tcac",
        "agac", "cgac", "ggac", "tgac", "atac", "ctac", "gtac", "ttac",
        "aacc", "cacc", "gacc", "tacc", "accc", "cccc", "gccc", "tccc",
        "agcc", "cgcc", "ggcc", "tgcc", "atcc", "ctcc", "gtcc", "ttcc",
        "aagc", "cagc", "gagc", "tagc", "acgc", "ccgc", "gcgc", "tcgc",
        "aggc", "cggc", "gggc", "tggc", "atgc", "ctgc", "gtgc", "ttgc",
        "aatc", "catc", "gatc", "tatc", "actc", "cctc", "gctc", "tctc",
        "agtc", "cgtc", "ggtc", "tgtc", "attc", "cttc", "gttc", "tttc",
        "aaag", "caag", "gaag", "taag", "acag", "ccag", "gcag", "tcag",
        "agag", "cgag", "ggag", "tgag", "atag", "ctag", "gtag", "ttag",
        "aacg", "cacg", "gacg", "tacg", "accg", "cccg", "gccg", "tccg",
        "agcg", "cgcg", "ggcg", "tgcg", "atcg", "ctcg", "gtcg", "ttcg",
        "aagg", "cagg", "gagg", "tagg", "acgg", "ccgg", "gcgg", "tcgg",
        "aggg", "cggg", "gggg", "tggg", "atgg", "ctgg", "gtgg", "ttgg",
        "aatg", "catg", "gatg", "tatg", "actg", "cctg", "gctg", "tctg",
        "agtg", "cgtg", "ggtg", "tgtg", "attg", "cttg", "gttg", "tttg",
        "aaat", "caat", "gaat", "taat", "acat", "ccat", "gcat", "tcat",
        "agat", "cgat", "ggat", "tgat", "atat", "ctat", "gtat", "ttat",
        "aact", "cact", "gact", "tact", "acct", "ccct", "gcct", "tcct",
        "agct", "cgct", "ggct", "tgct", "atct", "ctct", "gtct", "ttct",
        "aagt", "cagt", "gagt", "tagt", "acgt", "ccgt", "gcgt", "tcgt",
        "aggt", "cggt", "gggt", "tggt", "atgt", "ctgt", "gtgt", "ttgt",
        "aatt", "catt", "gatt", "tatt", "actt", "cctt", "gctt", "tctt",
        "agtt", "cgtt", "ggtt", "tgtt", "attt", "cttt", "gttt", "tttt",
        };
    while((start+pos)&3){
        dst[pos] = GPOINTER_TO_INT(FastaDB_SparseCache_get_func_2_BIT_LC(
                                   start+pos, page_data, user_data));
        pos++;
        }
    if(length > 4){
        while(end&3){
            dst[end-start-1]
                = GPOINTER_TO_INT(FastaDB_SparseCache_get_func_2_BIT_LC(
                                  end-1, page_data, user_data));
            end--;
            }
        }
    while(pos < end){
        word = table[(guchar)str[(start+pos)>>2]];
        dst[pos++] = word[0];
        dst[pos++] = word[1];
        dst[pos++] = word[2];
        dst[pos++] = word[3];
        }
#ifndef G_DISABLE_ASSERT
    for(i = 0; i < length; i++){
        g_assert(dst[i]
          == GPOINTER_TO_INT(FastaDB_SparseCache_get_func_2_BIT_LC(
                             start+i, page_data, user_data)));
        }
#endif /* G_DISABLE_ASSERT */
    return length;
    }

static void FastaDB_SparseCache_compress_2bit_lc(SparseCache_Page *page,
                                                 gint len){
    register gint i;
    register gchar *seq =  page->data;
    static gint index[ALPHABET_SIZE];
    static gboolean index_is_filled = FALSE;
    if(!index_is_filled){
        FastaDB_SparseCache_fill_index(index, "acgt");
        index_is_filled = TRUE;
        }
    page->data = g_new0(gchar, (len >> 2)+1);
    page->data_size = sizeof(gchar)*((len >> 2)+1);
    for(i = 0; i < len; i++)
        ((gchar*)page->data)[(i >> 2)] |= (index[(guchar)seq[i]] << ((i&3) << 1));
    g_free(seq);
    page->get_func = FastaDB_SparseCache_get_func_2_BIT_LC;
    page->copy_func = FastaDB_SparseCache_copy_func_2_BIT_LC;
    return;
    }

static gpointer FastaDB_SparseCache_get_func_2_BIT_UC(gint pos,
                                                      gpointer page_data,
                                                      gpointer user_data){
    register gint ch = ((gchar*)page_data)[(pos >> 2)],
                  num = (ch >> ((pos & 3) << 1)) & 3;
    static gchar *alphabet = "ACGT";
    return GINT_TO_POINTER((gint)alphabet[num]);
    }

static gint FastaDB_SparseCache_copy_func_2_BIT_UC(gint start, gint length,
                                                   gchar *dst, gpointer page_data,
                                                   gpointer user_data){
#ifndef G_DISABLE_ASSERT
    register gint i;
#endif /* G_DISABLE_ASSERT */
    register gint pos = 0, end = length;
    register gchar *str = page_data, *word;
    gchar *table[256] = {
        "AAAA", "CAAA", "GAAA", "TAAA", "ACAA", "CCAA", "GCAA", "TCAA",
        "AGAA", "CGAA", "GGAA", "TGAA", "ATAA", "CTAA", "GTAA", "TTAA",
        "AACA", "CACA", "GACA", "TACA", "ACCA", "CCCA", "GCCA", "TCCA",
        "AGCA", "CGCA", "GGCA", "TGCA", "ATCA", "CTCA", "GTCA", "TTCA",
        "AAGA", "CAGA", "GAGA", "TAGA", "ACGA", "CCGA", "GCGA", "TCGA",
        "AGGA", "CGGA", "GGGA", "TGGA", "ATGA", "CTGA", "GTGA", "TTGA",
        "AATA", "CATA", "GATA", "TATA", "ACTA", "CCTA", "GCTA", "TCTA",
        "AGTA", "CGTA", "GGTA", "TGTA", "ATTA", "CTTA", "GTTA", "TTTA",
        "AAAC", "CAAC", "GAAC", "TAAC", "ACAC", "CCAC", "GCAC", "TCAC",
        "AGAC", "CGAC", "GGAC", "TGAC", "ATAC", "CTAC", "GTAC", "TTAC",
        "AACC", "CACC", "GACC", "TACC", "ACCC", "CCCC", "GCCC", "TCCC",
        "AGCC", "CGCC", "GGCC", "TGCC", "ATCC", "CTCC", "GTCC", "TTCC",
        "AAGC", "CAGC", "GAGC", "TAGC", "ACGC", "CCGC", "GCGC", "TCGC",
        "AGGC", "CGGC", "GGGC", "TGGC", "ATGC", "CTGC", "GTGC", "TTGC",
        "AATC", "CATC", "GATC", "TATC", "ACTC", "CCTC", "GCTC", "TCTC",
        "AGTC", "CGTC", "GGTC", "TGTC", "ATTC", "CTTC", "GTTC", "TTTC",
        "AAAG", "CAAG", "GAAG", "TAAG", "ACAG", "CCAG", "GCAG", "TCAG",
        "AGAG", "CGAG", "GGAG", "TGAG", "ATAG", "CTAG", "GTAG", "TTAG",
        "AACG", "CACG", "GACG", "TACG", "ACCG", "CCCG", "GCCG", "TCCG",
        "AGCG", "CGCG", "GGCG", "TGCG", "ATCG", "CTCG", "GTCG", "TTCG",
        "AAGG", "CAGG", "GAGG", "TAGG", "ACGG", "CCGG", "GCGG", "TCGG",
        "AGGG", "CGGG", "GGGG", "TGGG", "ATGG", "CTGG", "GTGG", "TTGG",
        "AATG", "CATG", "GATG", "TATG", "ACTG", "CCTG", "GCTG", "TCTG",
        "AGTG", "CGTG", "GGTG", "TGTG", "ATTG", "CTTG", "GTTG", "TTTG",
        "AAAT", "CAAT", "GAAT", "TAAT", "ACAT", "CCAT", "GCAT", "TCAT",
        "AGAT", "CGAT", "GGAT", "TGAT", "ATAT", "CTAT", "GTAT", "TTAT",
        "AACT", "CACT", "GACT", "TACT", "ACCT", "CCCT", "GCCT", "TCCT",
        "AGCT", "CGCT", "GGCT", "TGCT", "ATCT", "CTCT", "GTCT", "TTCT",
        "AAGT", "CAGT", "GAGT", "TAGT", "ACGT", "CCGT", "GCGT", "TCGT",
        "AGGT", "CGGT", "GGGT", "TGGT", "ATGT", "CTGT", "GTGT", "TTGT",
        "AATT", "CATT", "GATT", "TATT", "ACTT", "CCTT", "GCTT", "TCTT",
        "AGTT", "CGTT", "GGTT", "TGTT", "ATTT", "CTTT", "GTTT", "TTTT",
        };
    while((start+pos)&3){
        dst[pos] = GPOINTER_TO_INT(FastaDB_SparseCache_get_func_2_BIT_UC(
                                   start+pos, page_data, user_data));
        pos++;
        }
    if(length > 4){
        while(end&3){
            dst[end-start-1]
                = GPOINTER_TO_INT(FastaDB_SparseCache_get_func_2_BIT_UC(end-1,
                                              page_data, user_data));
            end--;
            }
        }
    while(pos < end){
        word = table[(guchar)str[(start+pos)>>2]];
        dst[pos++] = word[0];
        dst[pos++] = word[1];
        dst[pos++] = word[2];
        dst[pos++] = word[3];
        }
#ifndef G_DISABLE_ASSERT
    for(i = 0; i < length; i++){
        g_assert(dst[i]
          == GPOINTER_TO_INT(FastaDB_SparseCache_get_func_2_BIT_UC(
                             start+i, page_data, user_data)));
        }
#endif /* G_DISABLE_ASSERT */
    return length;
    }

static void FastaDB_SparseCache_compress_2bit_uc(SparseCache_Page *page,
                                                 gint len){
    register gint i;
    register gchar *seq =  page->data;
    static gint index[ALPHABET_SIZE];
    static gboolean index_is_filled = FALSE;
    if(!index_is_filled){
        FastaDB_SparseCache_fill_index(index, "ACGT");
        index_is_filled = TRUE;
        }
    page->data = g_new0(gchar, (len >> 2)+1);
    page->data_size = sizeof(gchar)*((len >> 2)+1);
    for(i = 0; i < len; i++)
        ((gchar*)page->data)[(i >> 2)] |= (index[(guchar)seq[i]] << ((i&3) << 1));
    g_free(seq);
    page->get_func = FastaDB_SparseCache_get_func_2_BIT_UC;
    page->copy_func = FastaDB_SparseCache_copy_func_2_BIT_UC;
    return;
    }

/**/

void FastaDB_SparseCache_compress(SparseCache_Page *page, gint len){
    register gint i, mask = 0;
    register gchar *seq =  page->data;
    static FastaDB_SparseCache_Mask mask_index[ALPHABET_SIZE];
    static gboolean mask_is_filled = FALSE;
    if(!mask_is_filled){
        FastaDB_SparseCache_fill_mask(mask_index);
        mask_is_filled = TRUE;
        }
    for(i = 0; i < len; i++)
        mask |= mask_index[(guchar)seq[i]];
    if(mask & (1 << FastaDB_SparseCache_Mask_8_BIT)){
        page->copy_func = FastaDB_SparseCache_copy_func_8_BIT;
        return; /* Keep current 8 bit encoding */
        }
    if(mask & (mask-1)){ /* more than one mask bit is set */
        FastaDB_SparseCache_compress_4bit(page, len);
        return;
        }
    if(mask & (1 << FastaDB_SparseCache_Mask_2_BIT_LC)){
        FastaDB_SparseCache_compress_2bit_lc(page, len);
        return;
        }
    if(mask & (1 << FastaDB_SparseCache_Mask_2_BIT_UC)){
        FastaDB_SparseCache_compress_2bit_uc(page, len);
        return;
        }
    if(mask & (1 << FastaDB_SparseCache_Mask_0_BIT_LC)){
        g_free(page->data);
        page->data = NULL;
        page->data_size = 0;
        page->get_func = FastaDB_SparseCache_get_func_0_BIT_LC;
        page->copy_func = FastaDB_SparseCache_copy_func_0_BIT_LC;
        return;
        }
    if(mask & (1 << FastaDB_SparseCache_Mask_0_BIT_UC)){
        g_free(page->data);
        page->data = NULL;
        page->data_size = 0;
        page->get_func = FastaDB_SparseCache_get_func_0_BIT_UC;
        page->copy_func = FastaDB_SparseCache_copy_func_0_BIT_UC;
        return;
        }
    g_error("Bad FastaDB_SparseCache compression mask");
    return;
    }

static SparseCache_Page *FastaDB_SparseCache_fill_func(gint start,
                                                       gpointer user_data){
    register SparseCache_Page *page = g_new(SparseCache_Page, 1);
    register FastaDB_Key *fdbk = user_data;
    register gint ch, pos = 0, len;
    register gint db_start, db_end;
    if(fdbk->source->line_length < 1)
        g_error("Unknown or irregular fasta line length");
    db_start = fdbk->seq_offset
             + start
             + (start / (fdbk->source->line_length));
    len = MIN(SparseCache_PAGE_SIZE, fdbk->length-start);
    db_end = db_start
           + len
           + (len / (fdbk->source->line_length));
#ifdef USE_PTHREADS
    pthread_mutex_lock(&fdbk->source->cf->compoundfile_mutex);
#endif /* USE_PTHREADS */
    fdbk->location->pos += db_start;
    CompoundFile_Location_seek(fdbk->location);
    fdbk->location->pos -= db_start;
    page->data = g_new(gchar, len);
    page->get_func = FastaDB_SparseCache_get_func_8_BIT;
    page->copy_func = NULL; /* FIXME: temp */
    page->data_size = sizeof(gchar)*len;
    do {
        ch = CompoundFile_getc(fdbk->source->cf);
        g_assert(ch != EOF);
        g_assert(pos < len);
        if(!isspace(ch))
            ((gchar*)page->data)[pos++] = ch;
    } while(pos < len);
#ifdef USE_PTHREADS
    pthread_mutex_unlock(&fdbk->source->cf->compoundfile_mutex);
#endif /* USE_PTHREADS */
    FastaDB_SparseCache_compress(page, len);
    return page;
    }
/* FIXME: this is really inefficient: should use unbuffered read/write
 *        do this by adding CompoundFile_read()
 */
/* FIXME add suppport for compressed pages etc */

SparseCache *FastaDB_Key_get_SparseCache(FastaDB_Key *fdbk){
    register SparseCache *cache = SparseCache_create(fdbk->length,
                  FastaDB_SparseCache_fill_func, NULL, NULL, fdbk);
    return cache;
    }

/**/

FastaDB_Seq **FastaDB_all(gchar *path, Alphabet *alphabet,
                          FastaDB_Mask mask, guint *total){
    register FastaDB *fdb = FastaDB_open(path, alphabet);
    register GPtrArray *ptra = g_ptr_array_new();
    register FastaDB_Seq *fdbs, **fdbsa;
    register guint numseqs = 0;
    while((fdbs = FastaDB_next(fdb, mask))){
        g_ptr_array_add(ptra, fdbs);
        numseqs++;
        }
    g_ptr_array_add(ptra, NULL); /* NULL delimit */
    if(total)
        *total = numseqs; /* Record number of sequences */
    fdbsa = (FastaDB_Seq**)ptra->pdata;
    g_ptr_array_free(ptra, FALSE); /* Leave data intact */
    FastaDB_close(fdb);
    return fdbsa;
    }

void FastaDB_Seq_all_destroy(FastaDB_Seq **fdbs){
    register gint i;
    g_assert(fdbs);
    for(i = 0; fdbs[i]; i++)
        FastaDB_Seq_destroy(fdbs[i]);
    g_free(fdbs);
    return;
    }

gint FastaDB_Seq_print(FastaDB_Seq *fdbs, FILE *fp, FastaDB_Mask mask){
    register gint written = 0;
    if(mask & (FastaDB_Mask_ID|FastaDB_Mask_DEF|FastaDB_Mask_LEN)){
        written += fprintf(fp, ">%s", (mask & FastaDB_Mask_ID)
                             ?fdbs->seq->id:"[unknown]");
        if(fdbs->seq->def && (mask & FastaDB_Mask_DEF))
            written += fprintf(fp, " %s", fdbs->seq->def);
        if(mask & FastaDB_Mask_LEN)
            written += fprintf(fp, " [len:%d]", fdbs->seq->len);
        written += fprintf(fp, "\n");
        }
    if(mask & FastaDB_Mask_SEQ) /* to print seq, length must be set */
        written += Sequence_print_fasta_block(fdbs->seq, fp);
    return written;
    }

gint FastaDB_Seq_all_print(FastaDB_Seq **fdbs, FILE *fp,
                          FastaDB_Mask mask){
    register int i, written = 0;
    for(i = 0; fdbs[i]; i++)
        written += FastaDB_Seq_print(fdbs[i], fp, mask);
    return written;
    }

FastaDB_Seq *FastaDB_get_single(gchar *path, Alphabet *alphabet){
    register FastaDB *fdb = FastaDB_open(path, alphabet);
    register FastaDB_Seq *fdbs = FastaDB_next(fdb, FastaDB_Mask_ALL);
    if(!fdbs)
        g_error("No sequences found in [%s]\n", path);
    FastaDB_close(fdb);
    return fdbs;
    }

Alphabet_Type FastaDB_guess_type(gchar *path){
    register Alphabet *alphabet = Alphabet_create(Alphabet_Type_UNKNOWN,
                                                  FALSE);
    register FastaDB_Seq *fdbs = FastaDB_get_single(path, alphabet);
    register gchar *str = Sequence_get_str(fdbs->seq);
    register Alphabet_Type type = Alphabet_Type_guess(str);
    g_free(str);
    FastaDB_Seq_destroy(fdbs);
    Alphabet_destroy(alphabet);
    return type;
    }

