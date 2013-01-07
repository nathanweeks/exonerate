/****************************************************************\
*                                                                *
*  Library for word-neighbourhood generation                     *
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

#include <string.h> /* For strlen() */
#include <ctype.h>  /* For toupper() */
#include "wordhood.h"

/* Uses pruned DFS on the implicit wordhood trie.
 */

static gint WordHood_Submat_score_func(WordHood_Alphabet *wha,
                                       gchar *seq_a, gchar *seq_b){
    return Submat_lookup(wha->submat, *seq_a, *seq_b);
    }

static gboolean WordHood_Submat_is_valid_func(WordHood_Alphabet *wha,
                                              gchar *seq){
    return (wha->input_index[(guchar)*seq] != -1);
    }

static gint WordHood_Submat_index_func(WordHood_Alphabet *wha,
                                       gchar *seq){
    return wha->output_index[(guchar)*seq];
    }

static void WordHood_Alphabet_fill_index(gint *index, gchar *alphabet,
                                         gboolean case_sensitive_input){
    register gint i;
    for(i = 0; i < ALPHABETSIZE; i++)
        index[i] = -1;
    for(i = 0; alphabet[i]; i++){
        if(case_sensitive_input){
            index[(guchar)alphabet[i]] = i;
        } else {
            index[tolower(alphabet[i])] = i;
            index[toupper(alphabet[i])] = i;
            }
        }
    return;
    }

static WordHood_Alphabet *WordHood_Alphabet_create(gint advance,
                                        gchar *input_alphabet,
                                        gchar *output_alphabet,
                                        gboolean case_sensitive_input){
    register WordHood_Alphabet *wha = g_new(WordHood_Alphabet, 1);
    g_assert(input_alphabet);
    wha->ref_count = 1;
    wha->advance = advance;
    wha->member_list = g_ptr_array_new();
    WordHood_Alphabet_fill_index(wha->input_index, input_alphabet,
                                 case_sensitive_input);
    WordHood_Alphabet_fill_index(wha->output_index, output_alphabet,
                                 FALSE);
    return wha;
    }

WordHood_Alphabet *WordHood_Alphabet_create_from_Submat(
                      gchar *input_alphabet,
                      gchar *output_alphabet, Submat *submat,
                      gboolean case_sensitive_input){
    register WordHood_Alphabet *wha = WordHood_Alphabet_create(1,
            input_alphabet, output_alphabet, case_sensitive_input);
    register gint i;
    g_assert(submat);
    g_assert(output_alphabet);
    wha->submat = Submat_share(submat);
    wha->codon_submat = NULL;
    for(i = 0; output_alphabet[i]; i++)
        g_ptr_array_add(wha->member_list,
                        g_strndup(output_alphabet+i, 1));
    wha->score_func = WordHood_Submat_score_func;
    wha->is_valid_func = WordHood_Submat_is_valid_func;
    wha->index_func = WordHood_Submat_index_func;
    return wha;
    }

static gint WordHood_CodonSubmat_score_func(WordHood_Alphabet *wha,
                                            gchar *seq_a, gchar *seq_b){
    /*
    g_message("Codon score [%c%c%c][%c%c%c] = [%d]",
           seq_a[0], seq_a[1], seq_a[2],
           seq_b[0], seq_b[1], seq_b[2],
           CodonSubmat_lookup(wha->codon_submat, (guchar*)seq_a,
                                                 (guchar*)seq_b));
    */
    return CodonSubmat_lookup(wha->codon_submat, (guchar*)seq_a,
                                                 (guchar*)seq_b);
    }

static gboolean WordHood_CodonSubmat_is_valid_func(
                WordHood_Alphabet *wha, gchar *seq){
    if((wha->input_index[(guchar)seq[0]] == -1)
    || (wha->input_index[(guchar)seq[1]] == -1)
    || (wha->input_index[(guchar)seq[2]] == -1))
        return FALSE;
    return TRUE;
    }

static gint WordHood_CodonSubmat_index_func(WordHood_Alphabet *wha,
                                            gchar *seq){
    g_assert(wha->output_index[(guchar)seq[0]] != -1);
    g_assert(wha->output_index[(guchar)seq[1]] != -1);
    g_assert(wha->output_index[(guchar)seq[2]] != -1);
    return (wha->output_index[(guchar)seq[0]] << 4)
         | (wha->output_index[(guchar)seq[1]] << 2)
         |  wha->output_index[(guchar)seq[2]];
    }

WordHood_Alphabet *WordHood_Alphabet_create_from_CodonSubmat(
                    CodonSubmat *cs, gboolean case_sensitive_input){
    register gint a, b, c;
    register gchar *input_alphabet = "ABCDGHKMNRSTVWY",
                   *output_alphabet = "ACGT";
    register WordHood_Alphabet *wha = WordHood_Alphabet_create(3,
            input_alphabet, output_alphabet, case_sensitive_input);
    wha->submat = NULL;
    wha->codon_submat = CodonSubmat_share(cs);
    for(a = 0; a < 4; a++)
        for(b = 0; b < 4; b++)
            for(c = 0; c < 4; c++)
                g_ptr_array_add(wha->member_list,
                        g_strdup_printf("%c%c%c", output_alphabet[a],
                                                  output_alphabet[b],
                                                  output_alphabet[c]));
    wha->score_func = WordHood_CodonSubmat_score_func;
    wha->is_valid_func = WordHood_CodonSubmat_is_valid_func;
    wha->index_func = WordHood_CodonSubmat_index_func;
    return wha;
    }

WordHood_Alphabet *WordHood_Alphabet_share(WordHood_Alphabet *wha){
    wha->ref_count++;
    return wha;
    }

void WordHood_Alphabet_destroy(WordHood_Alphabet *wha){
    register gint i;
    if(--wha->ref_count)
        return;
    if(wha->submat)
        Submat_destroy(wha->submat);
    if(wha->codon_submat)
        CodonSubmat_destroy(wha->codon_submat);
    for(i = 0; i < wha->member_list->len; i++){
        g_assert(wha->member_list->pdata[i]);
        g_free(wha->member_list->pdata[i]);
        }
    g_ptr_array_free(wha->member_list, TRUE);
    g_free(wha);
    return;
    }

/**/

WordHood *WordHood_create(WordHood_Alphabet *wha,
                          gint threshold, gboolean use_dropoff){
    register WordHood *wh = g_new(WordHood, 1);
    wh->wha = WordHood_Alphabet_share(wha);
    wh->threshold = threshold;
    wh->use_dropoff = use_dropoff;
    wh->orig_word = NULL;
    wh->alloc_len = 64;
    wh->curr_word = g_new0(gchar, wh->alloc_len);
    wh->depth_threshold = g_new0(gint, wh->alloc_len);
    wh->curr_len = 0;
    wh->word_pos = 0;
    wh->curr_score = 0;
    return wh;
    }

void WordHood_destroy(WordHood *wh){
    WordHood_Alphabet_destroy(wh->wha);
    g_free(wh->depth_threshold);
    g_free(wh->curr_word);
    g_free(wh);
    return;
    }

static void WordHood_set_word(WordHood *wh, gchar *word, gint len,
                              gint threshold){
    register gint i;
    g_assert(len > 0);
    g_assert(!(len % wh->wha->advance));
    wh->orig_word = word;
    wh->curr_len = len;
    if(wh->alloc_len < len){
        wh->alloc_len = len;
        wh->curr_word = g_renew(gchar, wh->curr_word, wh->alloc_len);
        wh->depth_threshold = g_renew(gint, wh->depth_threshold,
                                       wh->alloc_len);
        wh->depth_threshold = g_renew(gint, wh->depth_threshold,
                                       wh->alloc_len);
        }
    wh->depth_threshold[wh->curr_len - wh->wha->advance]
           = threshold;
    for(i = wh->curr_len-(wh->wha->advance << 1);
        i >= 0; i -= wh->wha->advance){
        wh->depth_threshold[i]
            = wh->depth_threshold[i+wh->wha->advance]
            - wh->wha->score_func(wh->wha,
                                  word+i+wh->wha->advance,
                                  word+i+wh->wha->advance);
        }
    return;
    }

static void WordHood_Alphabet_info(WordHood_Alphabet *wha){
    g_message("WordHood_Alphabet:\n"
     "         ref_count [%d]\n"
     "         advance [%d]\n"
     "         member_list->len [%d]\n"
     "--\n",
         wha->ref_count,
         wha->advance,
         wha->member_list->len);
    return;
    }

void WordHood_info(WordHood *wh){
    WordHood_Alphabet_info(wh->wha);
    g_message("WordHood:\n"
     "         threshold [%d]\n"
     "         use_dropoff [%s]\n"
     "--\n",
        wh->threshold,
        wh->use_dropoff?"yes":"no");
    return;
    }

/**/

static gint WordHood_score_pos(WordHood *wh){
    return wh->wha->score_func(wh->wha,
                               wh->orig_word+wh->word_pos,
                               wh->curr_word+wh->word_pos);
    }

static gchar *WordHood_copy_member(WordHood *wh, gint member){
    g_assert(member < wh->wha->member_list->len);
    return strncpy(wh->curr_word+wh->word_pos,
            wh->wha->member_list->pdata[member], wh->wha->advance);
    }

static gboolean WordHood_Traverser_next(WordHood *wh){
    /* Ascend while at end of alphabet */
    while(!strncmp(wh->curr_word+wh->word_pos,
                  wh->wha->member_list->pdata
                  [wh->wha->member_list->len-1],
                  wh->wha->advance)){
        wh->curr_score -= WordHood_score_pos(wh);
        if(!wh->word_pos)
            return TRUE;
        wh->word_pos -= wh->wha->advance;
        }
    /* traverse */
    wh->curr_score -= WordHood_score_pos(wh);
    WordHood_copy_member(wh,
                         wh->wha->index_func(wh->wha,
                         wh->curr_word+wh->word_pos)+1);
    wh->curr_score += WordHood_score_pos(wh);
    return FALSE;
    }

static void WordHood_traverse_word(WordHood *wh, WordHood_Traverse_Func whtf,
                                   gpointer user_data){
    WordHood_copy_member(wh, 0); /* Copy first member */
    wh->word_pos = 0;
    wh->curr_score = WordHood_score_pos(wh);
    do {
        if(wh->curr_score < wh->depth_threshold[wh->word_pos]){
            if(WordHood_Traverser_next(wh))
                break;
        } else if(wh->word_pos == (wh->curr_len - wh->wha->advance)){
            /* is at leaf */
            if(whtf(wh->curr_word, wh->curr_score, user_data))
                break;
            if(WordHood_Traverser_next(wh))
                break;
        } else { /* descend */
            wh->word_pos += wh->wha->advance;
            WordHood_copy_member(wh, 0);
            wh->curr_score += WordHood_score_pos(wh);
            }
    } while(TRUE);
    return;
    }

/**/

static gboolean WordHood_word_is_valid(WordHood *wh){
    register gint i;
    g_assert(wh->orig_word);
    for(i = 0; i < wh->curr_len; i += wh->wha->advance)
        if(!wh->wha->is_valid_func(wh->wha, wh->orig_word+i))
            return FALSE;
    return TRUE;
    }
/* This can reject softmasked words when the index is case-sensitive */

static gint WordHood_score_word(WordHood *wh){
    register gint i, score = 0;
    for(i = 0; i < wh->curr_len; i += wh->wha->advance)
        score += wh->wha->score_func(wh->wha, wh->orig_word+i,
                                              wh->orig_word+i);
    return score;
    }

void WordHood_traverse(WordHood *wh, WordHood_Traverse_Func whtf,
                       gchar *word, gint len, gpointer user_data){
    register gint actual_threshold;
    g_assert(wh);
    g_assert(whtf);
    g_assert(word);
    wh->orig_word = word;
    wh->curr_len = len;
    if(!WordHood_word_is_valid(wh))
        return;
    if(wh->use_dropoff){
        /* FIXME: supply word_score to obviate this */
        actual_threshold = WordHood_score_word(wh) - wh->threshold;
    } else {
        actual_threshold = wh->threshold;
        }
    WordHood_set_word(wh, word, len, actual_threshold);
    WordHood_traverse_word(wh, whtf, user_data);
    wh->orig_word = NULL;
    return;
    }
/* FIXME: dropoff wordhood:
 *        change to only exclude negatively scoring words
 *        from the wordhood, rather than any ones containing
 *        non-alphabet words.
 */

/* FIXME: optimisations
 *          o use submat ordering
 *            (this will make it run much faster for proteins)
 *          o detect redundant words
 *          o detect mappable words in context of the submat
 *            (will probably only work well with nucleotides)
 */

