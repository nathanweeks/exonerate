/****************************************************************\
*                                                                *
*  Simple Sequence Object                                        *
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

#include <stdio.h>  /* For fopen() */
#include <string.h> /* For strlen() */
#include <ctype.h>  /* For toupper() */
#include <search.h> /* For tdelete(), tfind(), tsearch(), and twalk() */
#include <stdlib.h> /* For atoi() */

#include "sequence.h"
#include "lineparse.h"
#include "translate.h"

static Sequence_Annotation *Sequence_Annotation_create(gchar *id,
        Sequence_Strand strand, gint cds_start, gint cds_length){
    register Sequence_Annotation *annotation
     = g_new(Sequence_Annotation, 1);
    annotation->id = g_strdup(id);
    annotation->strand = strand;
    annotation->cds_start = cds_start;
    annotation->cds_length = cds_length;
    return annotation;
    }

static void Sequence_Annotation_destroy(
            Sequence_Annotation *annotation){
    g_free(annotation->id);
    g_free(annotation);
    return;
    }

static gint Sequence_Annotation_compare(gconstpointer a,
                                        gconstpointer b){
    register gchar *id_a = (gchar*)a,
                   *id_b = (gchar*)b;
    return strcmp(id_a, id_b);
    }

static GTree *Sequence_create_annotation_tree(gchar *path){
    void *tree = NULL;
    register FILE *fp = fopen(path, "r");
    register LineParse *lp;
    register gchar *id, strand_char;
    register gint cds_start = 0, cds_length = 0;
    register Sequence_Strand strand;
    register Sequence_Annotation *annotation;
    /**/
    if(!fp)
        g_error("Could not open annotation file [%s]", path);
    lp = LineParse_create(fp);
    while(LineParse_word(lp) != EOF){
        if((lp->word->len == 2) || (lp->word->len == 4)){
            id = lp->word->pdata[0];
            strand_char = *((gchar*)(lp->word->pdata[1]));
            strand = (strand_char == '+')
                   ? Sequence_Strand_FORWARD
                   : ((strand_char == '-')
                     ? Sequence_Strand_REVCOMP
                     : Sequence_Strand_UNKNOWN);
            if(strand == Sequence_Strand_UNKNOWN)
                g_error("Strand unknown for [%s]", id);
            if(lp->word->len == 4){
                cds_start = atoi(lp->word->pdata[2]) - 1;
                cds_length = atoi(lp->word->pdata[3]);
                }
            annotation = Sequence_Annotation_create(id, strand,
                                             cds_start, cds_length);
            g_assert(!tfind((void*)annotation->id, &tree, 
                          Sequence_Annotation_compare));
            tsearch((void*)annotation, &tree, Sequence_Annotation_compare);
            }
        }
    LineParse_destroy(lp);
    fclose(fp);
    return tree;
    }

static void Sequence_destroy_annotation_tree(void *tree){
    while (tree) {
        Sequence_Annotation* tmp = *(Sequence_Annotation **)tree; /*root node*/
        tdelete((void*)tmp, &tree, Sequence_Annotation_compare);
        Sequence_Annotation_destroy((Sequence_Annotation*)tmp);
        }
    }

static void Sequence_Argument_cleanup(gpointer user_data){
    register Sequence_ArgumentSet *sas = user_data;
    if(sas->annotation_tree){
        Sequence_destroy_annotation_tree(sas->annotation_tree);
        sas->annotation_tree = NULL;
        }
    return;
    }

Sequence_ArgumentSet *Sequence_ArgumentSet_create(Argument *arg){
    register ArgumentSet *as;
    static Sequence_ArgumentSet sas = {NULL};
    if(arg){
        as = ArgumentSet_create("Sequence Options");
        ArgumentSet_add_option(as, 'A', "annotation", "path",
           "Path to sequence annotation file", "none",
           Argument_parse_string, &sas.annotation_path);
        Argument_absorb_ArgumentSet(arg, as);
        Argument_add_cleanup(arg, Sequence_Argument_cleanup, &sas);
    } else {
        if(sas.annotation_path
        && strcmp(sas.annotation_path, "none")
        && (!sas.annotation_tree))
            sas.annotation_tree
                = Sequence_create_annotation_tree(sas.annotation_path);
        }
    return &sas;
    }

/**/

static gboolean Sequence_is_valid(Sequence *seq){
    register gint i, ch;
    g_assert(seq);
    g_assert(seq->id);
    if(seq->alphabet->type == Alphabet_Type_PROTEIN){
        g_assert(seq->strand == Sequence_Strand_UNKNOWN);
    } else if(seq->alphabet->type == Alphabet_Type_UNKNOWN){
        g_assert(seq->strand != Sequence_Strand_FORWARD);
        g_assert(seq->strand != Sequence_Strand_REVCOMP);
        }
    if(seq->data)
        for(i = 0; i < seq->len; i++){
            ch = Sequence_get_symbol(seq, i);
            if(!Alphabet_symbol_is_valid(seq->alphabet, ch))
                g_warning("Invalid symbol [%c](%d) at [%s][%d]",
                          ch, ch, seq->id, i);
            g_assert(Alphabet_symbol_is_valid(seq->alphabet, ch));
            }
    return TRUE;
    }

static void Sequence_IntMemory_data_destroy(gpointer data){
    g_free(data);
    return;
    }

static gint Sequence_IntMemory_get_symbol(gpointer data, gint pos){
    register gchar *seq = data;
    return seq[pos];
    }

static Sequence *Sequence_create_internal(gchar *id, gchar *def, guint len,
                          Sequence_Strand strand, Alphabet *alphabet){
    register Sequence *s = g_new0(Sequence, 1);
    register Sequence_ArgumentSet *sas
           = Sequence_ArgumentSet_create(NULL);
    void *tree_node;
    s->ref_count = 1;
    if(id)
        s->id  = g_strdup(id);
    if(def)
        s->def = g_strdup(def);
    s->strand = strand;
    if(alphabet)
        s->alphabet = Alphabet_share(alphabet);
    else
        s->alphabet = Alphabet_create(Alphabet_Type_UNKNOWN, FALSE);
    tree_node = tfind((void*)s->id, &sas->annotation_tree,
                       Sequence_Annotation_compare);
    s->annotation = tree_node ? *(Sequence_Annotation **)tree_node : NULL;
    s->len = len;
#ifdef USE_PTHREADS
    pthread_mutex_init(&s->seq_lock, NULL);
#endif /* USE_PTHREADS */
    return s;
    }

static void Sequence_ExtMemory_data_destroy(gpointer data){
    register SparseCache *cache = data;
    SparseCache_destroy(cache);
    return;
    }

static gint Sequence_ExtMemory_get_symbol(gpointer data, gint pos){
    register SparseCache *cache = data;
    return GPOINTER_TO_INT(SparseCache_get(cache, pos));
    }

Sequence *Sequence_create(gchar *id, gchar *def, gchar *seq, guint len,
                          Sequence_Strand strand, Alphabet *alphabet){
    register Sequence *s = Sequence_create_internal(id, def, len,
                                                    strand, alphabet);
    s->type = Sequence_Type_INTMEM;
    s->get_symbol = Sequence_IntMemory_get_symbol;
    if(seq){
        s->len = len ? len : strlen(seq);
        g_assert((!len) || (s->len == strlen(seq)));
        s->data = g_strndup(seq, s->len);
    } else {
        s->len = len;
        s->data = NULL;
        }
    g_assert(Sequence_is_valid(s));
    return s;
    }

Sequence *Sequence_create_extmem(gchar *id, gchar *def, guint len,
                          Sequence_Strand strand, Alphabet *alphabet,
                          SparseCache *cache){
    register Sequence *s = Sequence_create_internal(id, def, len,
                                                    strand, alphabet);
    s->type = Sequence_Type_EXTMEM;
    s->get_symbol = Sequence_ExtMemory_get_symbol;
    s->ref_count = 1;
    s->data = SparseCache_share(cache);
    return s;
    }

void Sequence_preload_extmem(Sequence *s){
    register gint i;
    if(s->type != Sequence_Type_EXTMEM)
        return;
    for(i = 0; i < s->len; i += SparseCache_PAGE_SIZE)
        Sequence_get_symbol(s, i);
    return;
    }

Sequence *Sequence_share(Sequence *s){
    g_assert(s);
    Sequence_lock(s);
    g_assert(s->ref_count);
    s->ref_count++;
    Sequence_unlock(s);
    return s;
    }

typedef struct {
    Sequence *sequence;
       guint  start;
} Sequence_Subseq;

static void Sequence_Subseq_data_destroy(gpointer data){
    register Sequence_Subseq *subseq = data;
    Sequence_destroy(subseq->sequence);
    g_free(subseq);
    return;
    }

static gint Sequence_Subseq_get_symbol(gpointer data, gint pos){
    register Sequence_Subseq *subseq = data;
    return Sequence_get_symbol(subseq->sequence, subseq->start+pos);
    }

Sequence *Sequence_subseq(Sequence *s, guint start, guint length){
    register Sequence *ns;
    register Sequence_Subseq *subseq;
    g_assert(s);
    g_assert(s->data);
    g_assert(s->len);
    g_assert(start < s->len);
    g_assert((start+length) <= s->len);
    if((start == 0) && (s->len == length))
        return Sequence_share(s);
    ns = Sequence_create(s->id, s->def, NULL, 0, s->strand, s->alphabet);
    g_free(ns->id);
    ns->id = g_strdup_printf("%s:subseq(%d,%d)", s->id, start, length);
    subseq = g_new(Sequence_Subseq, 1);
    subseq->sequence = Sequence_share(s);
    subseq->start = start;
    ns->type = Sequence_Type_SUBSEQ;
    ns->get_symbol = Sequence_Subseq_get_symbol;
    ns->data = subseq;
    ns->len = length;
    return ns;
    }

/**/

gint Sequence_print_fasta_block(Sequence *s, FILE *fp){
    register gint i, pos = 0, pause, width = 70, total = s->len;
    if(s->len){
        pos = 0;
        pause = pos+s->len-width;
        while(pos < pause){
            for(i = 0; i < width; i++)
                fputc(Sequence_get_symbol(s, pos+i), fp);
            fputc('\n', fp);
            pos += width;
            total++;
            }
        for(i = pos; i < s->len; i++)
            fputc(Sequence_get_symbol(s, i), fp);
        fputc('\n', fp);
        total++;
        }
#if 0
    /* FIXME: optimisation: use this for Sequence_Mode_IntMemory */
    register gchar *ptr, *pause;
    /**/
    if(s->seq){
        ptr = s->seq;
        pause = ptr+s->len-width;
        while(ptr < pause){
            ptr += fwrite(ptr, sizeof(gchar), width, fp);
            fputc('\n', fp);
            }
        ptr += fwrite(ptr, sizeof(gchar), s->seq+s->len-ptr, fp);
        fputc('\n', fp);
        }
#endif /* 0 */
    return total;
    }
/* FIXME: use fast version where possible
 *        (Sequence_strncpy each buffer, then write()
 *
 */

void Sequence_print_fasta(Sequence *s, FILE *fp, gboolean show_info){
    fprintf(fp, ">%s", s->id?s->id:"[unknown]");
    if(s->def)
        fprintf(fp, " %s", s->def);
    if(show_info){
        if(s->strand != Sequence_Strand_UNKNOWN)
            fprintf(fp, " [%s]", (s->strand == Sequence_Strand_FORWARD)
                               ?"forward":"revcomp");
        if(s->alphabet->type != Alphabet_Type_UNKNOWN)
            fprintf(fp, " [%s]",
                    Alphabet_Type_get_name(s->alphabet->type));
        if(s->len)
            fprintf(fp, " [length %d]", s->len);
        }
    fprintf(fp, "\n");
    Sequence_print_fasta_block(s, fp);
    return;
    }

static void Sequence_revcomp_data_destroy(gpointer data){
    register Sequence *sequence = data;
    Sequence_destroy(sequence);
    return;
    }

static gint Sequence_revcomp_get_symbol(gpointer data, gint pos){
    register Sequence *sequence = data;
    register gint ch = Sequence_get_symbol(sequence, sequence->len-pos-1);
    return sequence->alphabet->complement[ch];
    }

Sequence_Strand Sequence_Strand_revcomp(Sequence_Strand strand){
    g_assert((strand == Sequence_Strand_FORWARD)
           ||(strand == Sequence_Strand_REVCOMP));
    return (strand == Sequence_Strand_FORWARD)
           ? Sequence_Strand_REVCOMP
           : Sequence_Strand_FORWARD;
    }

void Sequence_revcomp_in_place(gchar *seq, guint length){
    register guchar *a, *z, swap;
    register gint pos;
    register Alphabet *alphabet = Alphabet_create(Alphabet_Type_DNA,
                                                  FALSE);
    register guchar *complement
        = Alphabet_get_filter_by_type(alphabet,
                                      Alphabet_Filter_Type_COMPLEMENT);
    for(a = (guchar*)seq, z = (guchar*)seq+length-1; a < z; a++, z--){
        swap = complement[*a];
        *a = complement[*z];
        *z = swap;
        }
    if(length & 1){ /* If odd length, complement the central base */
        pos = length >> 1;
        seq[pos] = complement[(guchar)seq[pos]];
        }
    Alphabet_destroy(alphabet);
    return;
    }
/* FIXME: optimisation: should avoid repeated creation of alphabet. */

void Sequence_reverse_in_place(gchar *seq, guint length){
    register guchar *a, *z, swap;
    for(a = (guchar*)seq, z = (guchar*)seq+length-1; a < z; a++, z--){
        swap = *a;
        *a = *z;
        *z = swap;
        }
    return;
    }

Sequence *Sequence_revcomp(Sequence *s){
    register Sequence *ns;
    register Sequence_Strand strand;
    /* Prevent creation of revcomp(revcomp(seq)) */
    if(s->type == Sequence_Type_REVCOMP)
        return Sequence_share((Sequence*)s->data);
    strand = Sequence_Strand_revcomp(s->strand);
    ns = Sequence_create(s->id, s->def, NULL, s->len, strand, s->alphabet);
    if(ns->def){
        g_free(ns->def);
        ns->def = g_strdup_printf("%s:[revcomp]", s->def);
    } else {
        ns->def = g_strdup("[revcomp]");
        }
    ns->data = Sequence_share(s);
    ns->type = Sequence_Type_REVCOMP;
    ns->get_symbol = Sequence_revcomp_get_symbol;
    return ns;
    }

/**/

void Sequence_filter_in_place(gchar *seq, guint length,
                              Alphabet *alphabet,
                              Alphabet_Filter_Type filter_type){
    register gint i;
    register const guchar *filter
                         = Alphabet_get_filter_by_type(alphabet,
                                                       filter_type);
    g_assert(filter);
    for(i = 0; i < length; i++){
        g_assert(filter[(guchar)seq[i]] != '-');
        seq[i] = filter[(guchar)seq[i]];
        }
    return;
    }

typedef struct {
                  Sequence *sequence;
                    guchar *filter;
      Alphabet_Filter_Type  filter_type;
} Sequence_Filter;

static void Sequence_Filter_data_destroy(gpointer data){
    register Sequence_Filter *sequence_filter = data;
    Sequence_destroy(sequence_filter->sequence);
    g_free(sequence_filter);
    return;
    }

static gint Sequence_Filter_get_symbol(gpointer data, gint pos){
    register Sequence_Filter *sequence_filter = data;
    return sequence_filter->filter
          [Sequence_get_symbol(sequence_filter->sequence, pos)];
    }

Sequence *Sequence_filter(Sequence *s,
                          Alphabet_Filter_Type filter_type){
    register Sequence *ns = Sequence_create(s->id, s->def, NULL,
                                        s->len, s->strand, s->alphabet);
    register Sequence_Filter *sequence_filter = g_new(Sequence_Filter, 1);
    g_assert(Alphabet_get_filter_by_type(s->alphabet, filter_type));
    g_free(ns->id);
    ns->id = g_strdup_printf("%s:filter(%s)", s->id,
                  Alphabet_Filter_Type_get_name(filter_type));
    ns->data = sequence_filter;
    ns->type = Sequence_Type_FILTER;
    ns->get_symbol = Sequence_Filter_get_symbol;
    sequence_filter->sequence = Sequence_share(s);
    sequence_filter->filter = Alphabet_get_filter_by_type(s->alphabet,
                                                          filter_type);
    sequence_filter->filter_type = filter_type;
    return ns;
    }
/* FIXME: optimisaton : if sequence is already filtered,
 *                      use merged filter with original seq
 */

gint Sequence_checksum(Sequence *s){
    register guint64 check = 0;
    register gint i, ch;
    register gchar *index =
    "--------------------------------------&---*---.-----------------"
    "@ABCDEFGHIJKLMNOPQRSTUVWXYZ------ABCDEFGHIJKLMNOPQRSTUVWXYZ---~-"
    "----------------------------------------------------------------"
    "----------------------------------------------------------------";
    for(i = 0; i < s->len; i++){
        ch = index[Sequence_get_symbol(s, i)];
        if(ch != '-')
            check += ((i % 57) + 1) * ch;
        }
    return check % 10000; /* gcg checksum */
    }

/**/

typedef struct {
    Sequence *sequence;
        gint  frame;
   Translate *translate;
} Sequence_Translation;

static void Sequence_Translation_data_destroy(gpointer data){
    register Sequence_Translation *translation = data;
    Sequence_destroy(translation->sequence);
    Translate_destroy(translation->translate);
    g_free(translation);
    return;
    }

static gint Sequence_translate_get_symbol(gpointer data, gint pos){
    register Sequence_Translation *translation = data;
    register gint p = (pos*3)+(translation->frame-1);
    return Translate_base(translation->translate,
                     Sequence_get_symbol(translation->sequence, p),
                     Sequence_get_symbol(translation->sequence, p+1),
                     Sequence_get_symbol(translation->sequence, p+2));
    }

Sequence *Sequence_translate(Sequence *s, Translate *translate, gint frame){
    register Alphabet *protein_alphabet
           = Alphabet_create(Alphabet_Type_PROTEIN, FALSE);
    register Sequence *ts = Sequence_create(s->id, s->def, NULL, 0,
                                            Sequence_Strand_UNKNOWN,
                                            protein_alphabet);
    register Sequence_Translation *translation
           = g_new(Sequence_Translation, 1);
    g_assert((frame >= 1) && (frame <= 3));
    if(ts->def){
        g_free(ts->def);
        ts->def = g_strdup_printf("%s:[translate(%d)]", s->def, frame);
    } else {
        ts->def = g_strdup_printf("[translate(%d)]", frame);
        }
    translation->sequence = Sequence_share(s);
    translation->frame = frame;
    translation->translate = Translate_share(translate);
    ts->len = (s->len-(frame-1))/3;
    ts->data = translation;
    ts->type = Sequence_Type_TRANSLATE;
    ts->get_symbol = Sequence_translate_get_symbol;
    Alphabet_destroy(protein_alphabet);
    return ts;
    }

/**/

#if 0
static void Sequence_print_type(Sequence *s){
    register Sequence_Subseq *subseq;
    register Sequence *revcomp;
    register Sequence_Filter *filter;
    register Sequence_Translation *translation;
    switch(s->type){
        case Sequence_Type_INTMEM:
            g_print("intmem");
            /* s->data is seq */
            break;
        case Sequence_Type_EXTMEM:
            g_print("extmem");
            /* s->data is SparseCache */
            break;
        case Sequence_Type_SUBSEQ:
            g_print("subseq:");
            subseq = s->data;
            Sequence_print_type(subseq->sequence);
            break;
        case Sequence_Type_REVCOMP:
            g_print("revcomp:");
            revcomp = s->data;
            Sequence_print_type(revcomp);
            break;
        case Sequence_Type_FILTER:
            g_print("filter:");
            filter = s->data;
            Sequence_print_type(filter->sequence);
            break;
        case Sequence_Type_TRANSLATE:
            g_print("translate:");
            translation = s->data;
            Sequence_print_type(translation->sequence);
            break;
        default:
            g_error("unknown sequence type [%d]", s->type);
            break;
        }
    return;
    }
#endif /* 0 */


void Sequence_strncpy(Sequence *s, gint start, gint length, gchar *dst){
    register gint i;
    register gchar *str;
    register Sequence_Subseq *subseq;
    register SparseCache *cache;
    g_assert(start >= 0);
    g_assert(s->len > 0);
    g_assert(start < s->len);
    g_assert((start+length) <= s->len);
    /* FIXME: g_warning("using slow implementation of [%s]", __FUNCTION__); */
#if 0
    g_print("Sequence_strncpy [");
    Sequence_print_type(s);
    g_print("] [%d] [%d] [%s]\n", start, length, s->id);
#endif /* 0 */
    switch(s->type){
        case Sequence_Type_INTMEM:
            str = s->data;
            for(i = 0; i < length; i++)
                dst[i] = str[start+i];
            return;
        case Sequence_Type_EXTMEM:
            cache = s->data;
            SparseCache_copy(cache, start, length, dst);
            break;
        case Sequence_Type_SUBSEQ:
            subseq = s->data;
            Sequence_strncpy(subseq->sequence, start+subseq->start, length, dst);
            return;
        default:
            for(i = 0; i < length; i++)
                dst[i] = Sequence_get_symbol(s, start+i);
            break;
        }
    return;
    }
/* FIXME: optimisation : implement fast versions
 *                       for specific sequence types
 *                       target types:
 *                           subseq:extmem
 *                           intmem
 */

void Sequence_strcpy(Sequence *s, gchar *dst){
    Sequence_strncpy(s, 0, s->len, dst);
    return;
    }

gchar *Sequence_get_substr(Sequence *s, gint start, gint length){
    register gchar *str = g_new(gchar, length+1);
    Sequence_strncpy(s, start, length, str);
    str[length] = '\0';
    return str;
    }

gchar *Sequence_get_str(Sequence *s){
    return Sequence_get_substr(s, 0, s->len);
    }
/* FIXME: optimisation: implement fast version
 *                      returning shared String for simple sequences
 */

/**/

void Sequence_destroy(Sequence *s){
    g_assert(s);
#ifdef USE_PTHREADS
    pthread_mutex_lock(&s->seq_lock);
#endif /* USE_PTHREADS */
    if(--s->ref_count){
#ifdef USE_PTHREADS
        pthread_mutex_unlock(&s->seq_lock);
#endif /* USE_PTHREADS */
        return;
        }
#ifdef USE_PTHREADS
    pthread_mutex_unlock(&s->seq_lock);
#endif /* USE_PTHREADS */
    if(s->id)
        g_free(s->id);
    if(s->def)
        g_free(s->def);
    if(s->data){
        switch(s->type){
            case Sequence_Type_INTMEM:
                Sequence_IntMemory_data_destroy(s->data);
                break;
            case Sequence_Type_EXTMEM:
                Sequence_ExtMemory_data_destroy(s->data);
                break;
            case Sequence_Type_SUBSEQ:
                Sequence_Subseq_data_destroy(s->data);
                break;
            case Sequence_Type_REVCOMP:
                Sequence_revcomp_data_destroy(s->data);
                break;
            case Sequence_Type_FILTER:
                Sequence_Filter_data_destroy(s->data);
                break;
            case Sequence_Type_TRANSLATE:
                Sequence_Translation_data_destroy(s->data);
                break;
            default:
                g_error("Unknown Sequence type [%d]", s->type);
                break;
            }
        }
    Alphabet_destroy(s->alphabet);
#ifdef USE_PTHREADS
    pthread_mutex_destroy(&s->seq_lock);
#endif /* USE_PTHREADS */
    g_free(s);
    return;
    }

Sequence *Sequence_mask(Sequence *s){
    register gboolean ok = TRUE;
    register Sequence *curr_seq = s, *new_seq, *prev_seq;
    register gint i;
    register GPtrArray *seq_list = g_ptr_array_new();
    register Sequence_Subseq *seq_subseq;
    register Sequence_Filter *seq_filter;
    register Sequence_Translation *seq_translation;
    /* Find the base sequence */
    do {
        switch(curr_seq->type){
            case Sequence_Type_INTMEM:
            case Sequence_Type_EXTMEM:
                ok = FALSE;
                break;
            case Sequence_Type_SUBSEQ:
                seq_subseq = curr_seq->data;
                g_ptr_array_add(seq_list, curr_seq);
                curr_seq = seq_subseq->sequence;
                break;
            case Sequence_Type_REVCOMP:
                g_ptr_array_add(seq_list, curr_seq);
                curr_seq = curr_seq->data;
                break;
            case Sequence_Type_FILTER:
                seq_filter = curr_seq->data;
                g_ptr_array_add(seq_list, curr_seq);
                curr_seq = seq_filter->sequence;
                break;
            case Sequence_Type_TRANSLATE:
                seq_translation = curr_seq->data;
                g_ptr_array_add(seq_list, curr_seq);
                curr_seq = seq_translation->sequence;
                break;
            default:
                g_error("Unknown Sequence Type [%d]", curr_seq->type);
                break;
            }
    } while(ok);
    /* Apply masking filter to base sequence */
    new_seq = Sequence_filter(curr_seq, Alphabet_Filter_Type_MASKED);
    /* Apply other transformations to filtered sequence copy */
    for(i = seq_list->len-1; i >= 0; i--){
        curr_seq = seq_list->pdata[i];
        prev_seq = new_seq;
        switch(curr_seq->type){
            case Sequence_Type_SUBSEQ:
                seq_subseq = curr_seq->data;
                new_seq = Sequence_subseq(prev_seq,
                                          seq_subseq->start,
                                          seq_subseq->sequence->len);
                break;
            case Sequence_Type_REVCOMP:
                new_seq = Sequence_revcomp(prev_seq);
                break;
            case Sequence_Type_FILTER:
                seq_filter = curr_seq->data;
                new_seq = Sequence_filter(prev_seq, seq_filter->filter_type);
                break;
            case Sequence_Type_TRANSLATE:
                seq_translation = curr_seq->data;
                new_seq = Sequence_translate(prev_seq,
                                             seq_translation->translate,
                                             seq_translation->frame);
                break;
            case Sequence_Type_INTMEM:
            case Sequence_Type_EXTMEM:
                g_error("impossible");
                break;
            default:
                g_error("Unknown Sequence type");
                break;
            }
        Sequence_destroy(prev_seq);
        }
    g_ptr_array_free(seq_list, TRUE);
    return new_seq;
    }

gsize Sequence_memory_usage(Sequence *s){
    register SparseCache *cache;
    register gsize data_memory = 0;
    switch(s->type){
        case Sequence_Type_INTMEM:
            data_memory = sizeof(gchar)*s->len;
            break;
        case Sequence_Type_EXTMEM:
            cache = s->data;
            data_memory = SparseCache_memory_usage(cache);
            break;
        default:
            data_memory = 0;
            break;
        }
    return sizeof(Sequence)
         + sizeof(gchar)*strlen(s->id)
         + sizeof(gchar)*(s->def?strlen(s->def):0)
         + sizeof(Alphabet)
         + data_memory;
    }


void Sequence_lock(Sequence *s){
    register Sequence_Subseq *subseq;
    register Sequence_Filter *filter;
    register Sequence_Translation *translation;
    register Sequence *revcomp;
    g_assert(s);
#ifdef USE_PTHREADS
    pthread_mutex_lock(&s->seq_lock);
    if(s->data){
        switch(s->type){
            case Sequence_Type_INTMEM:
                break;
            case Sequence_Type_EXTMEM:
                break;
            case Sequence_Type_SUBSEQ:
                subseq = s->data;
                Sequence_lock(subseq->sequence);
                break;
            case Sequence_Type_REVCOMP:
                revcomp = s->data;
                Sequence_lock(revcomp);
                break;
            case Sequence_Type_FILTER:
                filter = s->data;
                Sequence_lock(filter->sequence);
                break;
            case Sequence_Type_TRANSLATE:
                translation = s->data;
                Sequence_lock(translation->sequence);
                break;
            default:
                g_error("Unknown Sequence type [%d]", s->type);
                break;
            }
        }
#endif /* USE_PTHREADS */
    return;
    }

void Sequence_unlock(Sequence *s){
    register Sequence_Subseq *subseq;
    register Sequence_Filter *filter;
    register Sequence_Translation *translation;
    register Sequence *revcomp;
    g_assert(s);
#ifdef USE_PTHREADS
    pthread_mutex_unlock(&s->seq_lock);
    if(s->data){
        switch(s->type){
            case Sequence_Type_INTMEM:
                break;
            case Sequence_Type_EXTMEM:
                break;
            case Sequence_Type_SUBSEQ:
                subseq = s->data;
                Sequence_unlock(subseq->sequence);
                break;
            case Sequence_Type_REVCOMP:
                revcomp = s->data;
                Sequence_unlock(revcomp);
                break;
            case Sequence_Type_FILTER:
                filter = s->data;
                Sequence_unlock(filter->sequence);
                break;
            case Sequence_Type_TRANSLATE:
                translation = s->data;
                Sequence_unlock(translation->sequence);
                break;
            default:
                g_error("Unknown Sequence type [%d]", s->type);
                break;
            }
        }
#endif /* USE_PTHREADS */
    return;
    }



