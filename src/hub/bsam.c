/****************************************************************\
*                                                                *
*  BSAM: Big Sequence Alignment Manager                          *
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

#include <math.h>

#include "bsam.h"
#include "dejavu.h"
#include "comparison.h"

/**/

BSAM *BSAM_create(Comparison_Param *comparison_param,
                  gint saturate_threshold, gint verbosity){
    register BSAM *bsam = g_new(BSAM, 1);
    bsam->comparison_param = Comparison_Param_share(comparison_param);
    bsam->saturate_threshold = saturate_threshold;
    bsam->verbosity = verbosity;
    return bsam;
    }

void BSAM_destroy(BSAM *bsam){
    Comparison_Param_destroy(bsam->comparison_param);
    g_free(bsam);
    return;
    }

/**/

typedef struct {
       BSAM *bsam;
   Sequence *query;
   Sequence *target;
     HSPset *hspset;
     GArray *query_pos_list;
     GArray *target_pos_list;
     GArray *seed_list;
       gint  wordlen;
       gint  partition;
       gint  query_expect;
       gint  target_expect;
} BSAM_SeedSet;

static gint BSAM_get_expect(BSAM *bsam, Sequence *sequence,
                            gint wordlen){
    register gint expect;
    register gdouble expectation;
    register gint alphabet_size = 0;
    switch(sequence->alphabet->type){
        case Alphabet_Type_DNA:
            alphabet_size = 4;
            break;
        case Alphabet_Type_PROTEIN:
            alphabet_size = 20;
            break;
        default:
            g_error("Unknown alphabet type [%s]",
                    Alphabet_Type_get_name(sequence->alphabet->type));
            break;
        }
    expectation = 1.0/pow(alphabet_size, wordlen);
    expect = (gint)((expectation * (sequence->len - wordlen+1))
           + bsam->saturate_threshold);
    return expect;
    }

static BSAM_SeedSet *BSAM_SeedSet_create(BSAM *bsam,
        Comparison *comparison, HSPset *hspset){
    register BSAM_SeedSet *seedset = g_new(BSAM_SeedSet, 1);
    seedset->bsam = bsam;
    seedset->query = Sequence_share(comparison->query);
    seedset->target = Sequence_share(comparison->target);
    seedset->query_pos_list = g_array_new(FALSE, FALSE, sizeof(gint));
    seedset->target_pos_list = g_array_new(FALSE, FALSE, sizeof(gint));
    seedset->seed_list = g_array_new(FALSE, FALSE, sizeof(gint));
    if(hspset->param->match->query->is_translated)
       seedset->partition = 3 * ((comparison->query->len / 3)+1);
    else
        seedset->partition = comparison->query->len;
    seedset->hspset = HSPset_share(hspset);
    seedset->query_expect = BSAM_get_expect(bsam, comparison->query,
                                     seedset->hspset->param->wordlen);
    seedset->target_expect = BSAM_get_expect(bsam, comparison->target,
                                     seedset->hspset->param->wordlen);
    return seedset;
    }

static void BSAM_SeedSet_destroy(BSAM_SeedSet *seedset){
    Sequence_destroy(seedset->query);
    Sequence_destroy(seedset->target);
    HSPset_destroy(seedset->hspset);
    g_array_free(seedset->query_pos_list, TRUE);
    g_array_free(seedset->target_pos_list, TRUE);
    g_array_free(seedset->seed_list, TRUE);
    g_free(seedset);
    return;
    }

static void BSAM_SeedSet_empty_current(BSAM_SeedSet *seedset){
    register gint i, j;
    gint qpos, tpos;
    if(seedset->query_pos_list->len){
        if(seedset->target_pos_list->len){
            for(i = 0; i < seedset->query_pos_list->len; i++){
                qpos = g_array_index(seedset->query_pos_list,
                                     gint, i);
                for(j = 0; j < seedset->target_pos_list->len; j++){
                    tpos = g_array_index(seedset->target_pos_list,
                                         gint, j);
                    /* If seed is preceeded by a match, ignore.
                     * (Do not filter when using a saturate threshold
                     *  eg. poly-T at match start would prevent a seed).
                     */
                    if((!seedset->bsam->saturate_threshold)
                    && ((qpos > 0) && (tpos > 0)
                    && (Sequence_get_symbol(seedset->query, qpos-1)
                     == Sequence_get_symbol(seedset->target, tpos-1))))
                        continue;
                    g_array_append_val(seedset->seed_list, qpos);
                    g_array_append_val(seedset->seed_list, tpos);
                    }
                }
            g_array_set_size(seedset->target_pos_list, 0);
            }
        g_array_set_size(seedset->query_pos_list, 0);
        }
    /* (Can't have target positions without query positions) */
    return;
    }

static void BSAM_DejaVu_traverse(gint first_pos, gint curr_pos,
                                 gint length, gchar *seq, gint len,
                                 gpointer user_data){
    register BSAM_SeedSet *seedset = user_data;
    gint tpos;
    /* If first_pos is in target, ignore this word */
    if(first_pos > seedset->partition)
        return;
    if(first_pos == curr_pos) /* New word */
        BSAM_SeedSet_empty_current(seedset);
    if(curr_pos < seedset->partition){
        if((!seedset->bsam->saturate_threshold)
        || (seedset->query_pos_list->len <= seedset->query_expect))
           g_array_append_val(seedset->query_pos_list, curr_pos);
    } else {
        /* Collect 2nd list. (only if 1st list not empty) */
        /* Apply saturate threshold here
         * by emptying both lists.
         */
        if(seedset->bsam->saturate_threshold
        &&((seedset->query_pos_list->len
          * seedset->target_pos_list->len)
         > (seedset->query_expect * seedset->target_expect))){
            g_array_set_size(seedset->query_pos_list, 0);
            g_array_set_size(seedset->target_pos_list, 0);
            }
        /* FIXME: also need to apply saturatethreshold
         * to cap size of query_pos_list
         */
        if(seedset->query_pos_list->len){
            tpos = curr_pos - seedset->partition - 1;
            g_array_append_val(seedset->target_pos_list, tpos);
            }
        }
    return;
    }
/*
 *        <=>-<=>-<=>-<=>
 *   init:012-345-678-9AB
 *  trans:0369-147A-258B-
 *
 *  frame = (pos < len)?0:(pos < (len<<1))?1:2;
 *  frame_start = (frame * (aa_len + 1));
 *  codon = pos - frame_start
 *  orig_pos = (codon * 3) + frame
 *
 */

static guint BSAM_SeedSet_aapos2dnapos(guint pos, gint dna_len){
    register gint frame, frame_start, codon, aa_len;
    register guint dna_pos;
    aa_len = dna_len / 3;
    if(pos < ((aa_len+1) << 1))
        if(pos < (aa_len+1))
            frame = 0;
        else
            frame = 1;
    else
        frame = 2;
    frame_start = (frame * (aa_len + 1));
    codon = pos - frame_start;
    dna_pos = (codon * 3) + frame;
    g_assert(dna_pos >= 0);
    g_assert(dna_pos < dna_len);
    return dna_pos;
    }

static void BSAM_SeedSet_build_hspset(BSAM *bsam, BSAM_SeedSet *seedset,
                                      Comparison *comparison){
    register gint i;
    guint pos;
    /* Convert query coordinates */
    if(seedset->hspset->param->match->query->is_translated){
        for(i = 0; i < seedset->seed_list->len; i+=2){
            pos = g_array_index(seedset->seed_list, guint, i);
            pos = BSAM_SeedSet_aapos2dnapos(pos,
                                            comparison->query->len);
            g_array_index(seedset->seed_list, guint, i) = pos;
            }
        }
    /* Convert target coordinates */
    if(seedset->hspset->param->match->target->is_translated){
        for(i = 0; i < seedset->seed_list->len; i+=2){
            pos = g_array_index(seedset->seed_list, guint, i+1);
            pos = BSAM_SeedSet_aapos2dnapos(pos,
                                            comparison->target->len);
            g_array_index(seedset->seed_list, guint, i+1) = pos;
            }
        }
    g_assert(!(seedset->seed_list->len & 1));
    if(seedset->seed_list->len){
        HSPset_seed_all_hsps(seedset->hspset,
                         (guint*)seedset->seed_list->data,
                         (guint)(seedset->seed_list->len >> 1));
        }
    return;
    }

static void BSAM_BigSeq_add_sequence(GString *seq,
                                     Sequence *input,
                                     Translate *translate,
                                     gboolean is_translated,
                                     guchar *filter){
    register gint i, j, max_aa_len;
    register Sequence *aa_seq;
    if(is_translated){
        g_assert(translate);
        max_aa_len = input->len / 3;
        for(i = 0; i < 3; i++){
            aa_seq = Sequence_translate(input, translate, i+1);
            Sequence_strcpy(aa_seq, seq->str+seq->len);
            for(j = aa_seq->len; j < max_aa_len; j++)
                g_string_append_c(seq, '-');
            g_string_append_c(seq, '-');
            Sequence_destroy(aa_seq);
            }
    } else {
        Sequence_strcpy(input, seq->str+seq->len);
        seq->len += input->len;
        }
    return;
    }

static GString *BSAM_create_BigSeq(BSAM *bsam, Sequence *query,
                                               Sequence *target,
                                               BSAM_SeedSet *seedset){
    register Match *match = seedset->hspset->param->match;
    register GString *bigseq = g_string_sized_new(query->len+target->len+1);
    if(bsam->verbosity > 1)
        g_message("Preparing BigSeq comparison of [%s] and [%s]",
            query->id, target->id);
    BSAM_BigSeq_add_sequence(bigseq, query,
            match->mas->translate,
            match->query->is_translated,
            match->query->alphabet->masked);
    g_assert(bigseq->len == seedset->partition);
    g_string_append_c(bigseq, '-');
    BSAM_BigSeq_add_sequence(bigseq, target,
            match->mas->translate,
            match->target->is_translated,
            match->target->alphabet->masked);
    return bigseq;
    }

static void BSAM_build_HSPset(BSAM *bsam, Comparison *comparison,
                              HSPset *hspset){
    register DejaVu *dejavu;
    register GString *bigseq;
    register BSAM_SeedSet *seedset = BSAM_SeedSet_create(bsam,
                                             comparison, hspset);
    register guchar *filter = Alphabet_get_filter_by_type(
        seedset->hspset->param->match->comparison_alphabet,
        Alphabet_Filter_Type_NON_AMBIG); /* Will include softmasking */
    bigseq = BSAM_create_BigSeq(bsam, comparison->query,
                                      comparison->target, seedset);
    if(bsam->verbosity > 1)
        g_message("Finding BigSeq seeds");
    dejavu = DejaVu_create(bigseq->str, bigseq->len);
    DejaVu_traverse(dejavu, hspset->param->wordlen,
                            hspset->param->wordlen,
                            BSAM_DejaVu_traverse, seedset, filter,
                            bsam->verbosity);
    BSAM_SeedSet_empty_current(seedset);
    DejaVu_destroy(dejavu);
    g_string_free(bigseq, TRUE);
    if(bsam->verbosity > 1)
        g_message("Building BigSeq HSPs from [%d] seeds",
                 (seedset->seed_list->len >> 1));
    BSAM_SeedSet_build_hspset(bsam, seedset, comparison);
    BSAM_SeedSet_destroy(seedset);
    if(bsam->verbosity > 1)
        g_message("Building BigSeq alignment from [%d] hsps",
                  hspset->hsp_list->len);
    return;
    }

Comparison *BSAM_compare(BSAM *bsam, Sequence *query, Sequence *target){
    register Comparison *comparison
           = Comparison_create(bsam->comparison_param, query, target);
    if(bsam->comparison_param->dna_hsp_param)
        BSAM_build_HSPset(bsam, comparison, comparison->dna_hspset);
    if(bsam->comparison_param->protein_hsp_param)
        BSAM_build_HSPset(bsam, comparison, comparison->protein_hspset);
    if(bsam->comparison_param->codon_hsp_param)
        BSAM_build_HSPset(bsam, comparison, comparison->codon_hspset);
    if(!Comparison_has_hsps(comparison)){
        Comparison_destroy(comparison);
        return NULL;
        }
    return comparison;
    }

/**/

