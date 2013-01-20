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

#ifndef INCLUDED_SEQUENCE_H
#define INCLUDED_SEQUENCE_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <stdio.h>
#include <glib.h>

#ifdef USE_PTHREADS
#include <pthread.h>
#endif /* USE_PTHREADS */

#include "alphabet.h"
#include "translate.h"
#include "sparsecache.h"

typedef enum {
  Sequence_Type_INTMEM,
  Sequence_Type_EXTMEM,
  Sequence_Type_SUBSEQ,
  Sequence_Type_REVCOMP,
  Sequence_Type_FILTER,
  Sequence_Type_TRANSLATE
} Sequence_Type;

typedef enum {
    Sequence_Strand_FORWARD,
    Sequence_Strand_REVCOMP,
    Sequence_Strand_UNKNOWN
} Sequence_Strand;

typedef struct {
              gchar *id;
    Sequence_Strand strand;
               gint cds_start;
               gint cds_length;
} Sequence_Annotation;

typedef struct {
    gchar *annotation_path;
    void *annotation_tree;
} Sequence_ArgumentSet;

Sequence_ArgumentSet *Sequence_ArgumentSet_create(Argument *arg);

/**/

Sequence_Strand Sequence_Strand_revcomp(Sequence_Strand strand);

#define Sequence_Strand_get_as_char(strand) \
    (((strand) == Sequence_Strand_UNKNOWN)?'.' \
    :((strand) == Sequence_Strand_FORWARD)?'+':'-')
/* Shows strand of sequence as [-.+] */

#define Sequence_Strand_get_as_string(strand) \
    (((strand) == Sequence_Strand_UNKNOWN)?"unknown" \
    :((strand) == Sequence_Strand_FORWARD)?"forward":"revcomp")

#define Sequence_get_strand_as_char(seq) \
    Sequence_Strand_get_as_char((seq)->strand)

typedef struct {
                  guint   ref_count;
                  gchar  *id;
                  gchar  *def;
                  guint   len;
        Sequence_Strand   strand;
               Alphabet  *alphabet;
    Sequence_Annotation  *annotation;
               gpointer   data;
          Sequence_Type   type;
                   gint (*get_symbol)(gpointer data, gint pos);
#ifdef USE_PTHREADS
        pthread_mutex_t   seq_lock;
#endif /* USE_PTHREADS */
} Sequence;

#define Sequence_get_symbol(sequence, pos) \
    ((sequence)->get_symbol((sequence)->data, pos))

#define Sequence_get_symbol_TEMP(sequence, pos) \
    ((sequence)->get_symbol((sequence)->data, pos))

Sequence *Sequence_create(gchar *id, gchar *def, gchar *seq, guint len,
                          Sequence_Strand strand, Alphabet *alphabet);
/* Any of the Sequence_create arguments, id, def and seq
 * may be NULL, and len may also be zero.
 *
 * If len is non-zero, the seq argument need not be NULL terminated.
 * (However it will be NULL terminated in the Sequence object).
 *
 * If alphabet is NULL, (UNKNOWN|UNMASKED) is assumed.
 */
Sequence *Sequence_create_extmem(gchar *id, gchar *def, guint len,
                          Sequence_Strand strand, Alphabet *alphabet,
                          SparseCache *cache);
    void  Sequence_preload_extmem(Sequence *s);
    void  Sequence_destroy(Sequence *s);
Sequence *Sequence_share(Sequence *s);
Sequence *Sequence_subseq(Sequence *s, guint start, guint length);
    void  Sequence_print_fasta(Sequence *s, FILE *fp,
                               gboolean show_info);
    gint  Sequence_print_fasta_block(Sequence *s, FILE *fp);

    void  Sequence_revcomp_in_place(gchar *seq, guint length);
Sequence *Sequence_revcomp(Sequence *s);
    void  Sequence_reverse_in_place(gchar *seq, guint length);
    gint  Sequence_checksum(Sequence *s);

    void  Sequence_strncpy(Sequence *s, gint start, gint length, gchar *dst);
    void  Sequence_strcpy(Sequence *s, gchar *dst);
   gchar *Sequence_get_substr(Sequence *s, gint start, gint length);
   gchar *Sequence_get_str(Sequence *s);
   gsize  Sequence_memory_usage(Sequence *s);

/**/

    void  Sequence_filter_in_place(gchar *seq, guint length,
                                   Alphabet *alphabet,
                                   Alphabet_Filter_Type filter_type);
Sequence *Sequence_filter(Sequence *s,
                          Alphabet_Filter_Type filter_type);

Sequence *Sequence_translate(Sequence *s, Translate *translate, gint frame);

Sequence *Sequence_mask(Sequence *s);
void Sequence_lock(Sequence *s);
void Sequence_unlock(Sequence *s);

/**/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* INCLUDED_SEQUENCE_H */


