/****************************************************************\
*                                                                *
*  Library for HSP sets (high-scoring segment pairs)             *
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

#include <string.h> /* For strlen()  */
#include <ctype.h>  /* For tolower() */
#include <stdlib.h> /* For qsort()   */

#include "exonerate_util.h"
#include "hspset.h"

HSPset_ArgumentSet *HSPset_ArgumentSet_create(Argument *arg){
    register ArgumentSet *as;
    static HSPset_ArgumentSet has = {0};
    if(arg){
        as = ArgumentSet_create("HSP creation options");
        ArgumentSet_add_option(as, '\0', "hspfilter", NULL,
                "Aggressive HSP filtering level", "0",
                Argument_parse_int, &has.filter_threshold);
        ArgumentSet_add_option(as, '\0', "useworddropoff", NULL,
                "Use word neighbourhood dropoff", "TRUE",
                Argument_parse_boolean, &has.use_wordhood_dropoff);
        ArgumentSet_add_option(as, '\0', "seedrepeat", NULL,
                "Seeds per diagonal required for HSP seeding", "1",
                Argument_parse_int, &has.seed_repeat);
        /**/
        ArgumentSet_add_option(as, 0, "dnawordlen", "bp",
            "Wordlength for DNA words", "12",
            Argument_parse_int, &has.dna_wordlen);
        ArgumentSet_add_option(as, 0, "proteinwordlen", "aa",
            "Wordlength for protein words", "6",
            Argument_parse_int, &has.protein_wordlen);
        ArgumentSet_add_option(as, 0, "codonwordlen", "bp",
            "Wordlength for codon words", "12",
            Argument_parse_int, &has.codon_wordlen);
        /**/
        ArgumentSet_add_option(as, 0, "dnahspdropoff", "score",
            "DNA HSP dropoff score", "30",
            Argument_parse_int, &has.dna_hsp_dropoff);
        ArgumentSet_add_option(as, 0, "proteinhspdropoff", "score",
            "Protein HSP dropoff score", "20",
            Argument_parse_int, &has.protein_hsp_dropoff);
        ArgumentSet_add_option(as, 0, "codonhspdropoff", "score",
            "Codon HSP dropoff score", "40",
            Argument_parse_int, &has.codon_hsp_dropoff);
        /**/
        ArgumentSet_add_option(as, 0, "dnahspthreshold", "score",
            "DNA HSP threshold score", "75",
            Argument_parse_int, &has.dna_hsp_threshold);
        ArgumentSet_add_option(as, 0, "proteinhspthreshold", "score",
            "Protein HSP threshold score", "30",
            Argument_parse_int, &has.protein_hsp_threshold);
        ArgumentSet_add_option(as, 0, "codonhspthreshold", "score",
            "Codon HSP threshold score", "50",
            Argument_parse_int, &has.codon_hsp_threshold);
        /**/
        ArgumentSet_add_option(as, 0, "dnawordlimit", "score",
            "Score limit for dna word neighbourhood", "0",
            Argument_parse_int, &has.dna_word_limit);
        ArgumentSet_add_option(as, 0, "proteinwordlimit", "score",
            "Score limit for protein word neighbourhood", "4",
            Argument_parse_int, &has.protein_word_limit);
        ArgumentSet_add_option(as, 0, "codonwordlimit", "score",
            "Score limit for codon word neighbourhood", "4",
            Argument_parse_int, &has.codon_word_limit);
        /**/
        ArgumentSet_add_option(as, '\0', "geneseed", "threshold",
             "Geneseed Threshold", "0",
             Argument_parse_int, &has.geneseed_threshold);
        ArgumentSet_add_option(as, '\0', "geneseedrepeat", "number",
             "Seeds per diagonal required for geneseed HSP seeding", "3",
             Argument_parse_int, &has.geneseed_repeat);
        /**/
        Argument_absorb_ArgumentSet(arg, as);
        }
    return &has;
    }

/**/

static gboolean HSP_check_positions(HSPset *hsp_set,
                                    gint query_pos, gint target_pos){
    g_assert(query_pos >= 0);
    g_assert(target_pos >= 0);
    g_assert(query_pos < hsp_set->query->len);
    g_assert(target_pos < hsp_set->target->len);
    return TRUE;
    }

static gboolean HSP_check(HSP *hsp){
    g_assert(hsp);
    g_assert(hsp->hsp_set);
    g_assert(hsp->length);
    g_assert(HSP_query_end(hsp) <= hsp->hsp_set->query->len);
    g_assert(HSP_target_end(hsp) <= hsp->hsp_set->target->len);
    return TRUE;
    }

void HSP_Param_set_wordlen(HSP_Param *hsp_param, gint wordlen){
    if(wordlen <= 0)
        g_error("Wordlength must be greater than zero");
    hsp_param->wordlen = wordlen;
    hsp_param->seedlen = hsp_param->wordlen
                       / hsp_param->match->query->advance;
    return;
    }

/**/

void HSP_Param_set_dna_hsp_threshold(HSP_Param *hsp_param,
                                     gint dna_hsp_threshold){
    if(hsp_param->match->type == Match_Type_DNA2DNA)
        hsp_param->threshold = dna_hsp_threshold;
    return;
    }

void HSP_Param_set_protein_hsp_threshold(HSP_Param *hsp_param,
                                         gint protein_hsp_threshold){
    if((hsp_param->match->type == Match_Type_PROTEIN2PROTEIN)
    || (hsp_param->match->type == Match_Type_PROTEIN2DNA)
    || (hsp_param->match->type == Match_Type_DNA2PROTEIN))
        hsp_param->threshold = protein_hsp_threshold;
    return;
    }

void HSP_Param_set_codon_hsp_threshold(HSP_Param *hsp_param,
                                       gint codon_hsp_threshold){
    if(hsp_param->match->type == Match_Type_CODON2CODON)
        hsp_param->threshold = codon_hsp_threshold;
    return;
    }

/**/

static void HSP_Param_refresh_wordhood(HSP_Param *hsp_param){
    register WordHood_Alphabet *wha = NULL;
    register Submat *submat;
    if(hsp_param->wordhood)
        WordHood_destroy(hsp_param->wordhood);
    if(hsp_param->has->use_wordhood_dropoff && (!hsp_param->wordlimit)){
        hsp_param->wordhood = NULL;
    } else {
        submat = (hsp_param->match->type == Match_Type_DNA2DNA)
               ? hsp_param->match->mas->dna_submat
               : hsp_param->match->mas->protein_submat;
        wha = WordHood_Alphabet_create_from_Submat(
                (gchar*)hsp_param->match->comparison_alphabet->member,
                (gchar*)hsp_param->match->comparison_alphabet->member,
                submat, FALSE);
        g_assert(wha);
        hsp_param->wordhood = WordHood_create(wha, hsp_param->wordlimit,
                                       hsp_param->has->use_wordhood_dropoff);
        WordHood_Alphabet_destroy(wha);
        }
    return;
    }
/* FIXME: optimisation: do not free/recreate wordhood when unnecessary */

void HSP_Param_set_dna_word_limit(HSP_Param *hsp_param,
                                  gint dna_word_limit){
    if(hsp_param->match->type == Match_Type_DNA2DNA){
        hsp_param->wordlimit = dna_word_limit;
        HSP_Param_refresh_wordhood(hsp_param);
        }
    return;
    }

void HSP_Param_set_protein_word_limit(HSP_Param *hsp_param,
                                      gint protein_word_limit){
    if((hsp_param->match->type == Match_Type_PROTEIN2PROTEIN)
    || (hsp_param->match->type == Match_Type_PROTEIN2DNA)
    || (hsp_param->match->type == Match_Type_DNA2PROTEIN)){
        hsp_param->wordlimit = protein_word_limit;
        HSP_Param_refresh_wordhood(hsp_param);
        }
    return;
    }

void HSP_Param_set_codon_word_limit(HSP_Param *hsp_param,
                                    gint codon_word_limit){
    if(hsp_param->match->type == Match_Type_CODON2CODON){
        hsp_param->threshold = codon_word_limit;
        HSP_Param_refresh_wordhood(hsp_param);
        }
    return;
    }

/**/

void HSP_Param_set_dna_hsp_dropoff(HSP_Param *hsp_param,
                                   gint dna_hsp_dropoff){
    if(hsp_param->match->type == Match_Type_DNA2DNA)
        hsp_param->dropoff = dna_hsp_dropoff;
    return;
    }

void HSP_Param_set_protein_hsp_dropoff(HSP_Param *hsp_param,
                                       gint protein_hsp_dropoff){
    if((hsp_param->match->type == Match_Type_PROTEIN2PROTEIN)
    || (hsp_param->match->type == Match_Type_PROTEIN2DNA)
    || (hsp_param->match->type == Match_Type_DNA2PROTEIN))
        hsp_param->dropoff  = protein_hsp_dropoff;
    return;
    }

void HSP_Param_set_codon_hsp_dropoff(HSP_Param *hsp_param,
                                     gint codon_hsp_dropoff){
    if(hsp_param->match->type == Match_Type_CODON2CODON)
        hsp_param->dropoff = codon_hsp_dropoff;
    return;
    }

/**/

void HSP_Param_set_hsp_threshold(HSP_Param *hsp_param,
                                 gint hsp_threshold){
    hsp_param->threshold = hsp_threshold;
    return;
    }

void HSP_Param_set_seed_repeat(HSP_Param *hsp_param,
                               gint seed_repeat){
    hsp_param->seed_repeat = seed_repeat;
    return;
    }

HSP_Param *HSP_Param_create(Match *match, gboolean use_horizon){
    register HSP_Param *hsp_param = g_new(HSP_Param, 1);
    hsp_param->thread_ref = ThreadRef_create();
    hsp_param->has = HSPset_ArgumentSet_create(NULL);
    hsp_param->match = match;
    hsp_param->seed_repeat = hsp_param->has->seed_repeat;
    switch(match->type){
        case Match_Type_DNA2DNA:
            hsp_param->dropoff = hsp_param->has->dna_hsp_dropoff;
            hsp_param->threshold = hsp_param->has->dna_hsp_threshold;
            HSP_Param_set_wordlen(hsp_param, hsp_param->has->dna_wordlen);
            hsp_param->wordlimit = hsp_param->has->dna_word_limit;
            break;
        case Match_Type_PROTEIN2PROTEIN: /*fallthrough*/
        case Match_Type_PROTEIN2DNA:     /*fallthrough*/
        case Match_Type_DNA2PROTEIN:
            hsp_param->dropoff = hsp_param->has->protein_hsp_dropoff;
            hsp_param->threshold
                = hsp_param->has->protein_hsp_threshold;
            HSP_Param_set_wordlen(hsp_param, hsp_param->has->protein_wordlen);
            hsp_param->wordlimit = hsp_param->has->protein_word_limit;
            break;
        case Match_Type_CODON2CODON:
            hsp_param->dropoff = hsp_param->has->codon_hsp_dropoff;
            hsp_param->threshold = hsp_param->has->codon_hsp_threshold;
            HSP_Param_set_wordlen(hsp_param, hsp_param->has->codon_wordlen);
            hsp_param->wordlimit = hsp_param->has->codon_word_limit;
            g_assert(!(hsp_param->wordlen % 3));
            break;
        default:
            g_error("Bad Match_Type [%d]", match->type);
            break;
        }
    hsp_param->use_horizon = use_horizon;
#ifdef USE_PTHREADS
    pthread_mutex_init(&hsp_param->hsp_recycle_lock, NULL);
#endif /* USE_PTHREADS */
    hsp_param->hsp_recycle = RecycleBin_create("HSP", sizeof(HSP), 64);
    hsp_param->wordhood = NULL;
    HSP_Param_refresh_wordhood(hsp_param);
    return hsp_param;
    }

void HSP_Param_destroy(HSP_Param *hsp_param){
    g_assert(hsp_param);
    if(ThreadRef_destroy(hsp_param->thread_ref))
        return;
    if(hsp_param->wordhood)
        WordHood_destroy(hsp_param->wordhood);
#ifdef USE_PTHREADS
    pthread_mutex_lock(&hsp_param->hsp_recycle_lock);
#endif /* USE_PTHREADS */
    RecycleBin_destroy(hsp_param->hsp_recycle);
#ifdef USE_PTHREADS
    pthread_mutex_unlock(&hsp_param->hsp_recycle_lock);
    pthread_mutex_destroy(&hsp_param->hsp_recycle_lock);
#endif /* USE_PTHREADS */
    g_free(hsp_param);
    return;
    }

HSP_Param *HSP_Param_share(HSP_Param *hsp_param){
    g_assert(hsp_param);
    ThreadRef_share(hsp_param->thread_ref);
    return hsp_param;
    }

HSP_Param *HSP_Param_swap(HSP_Param *hsp_param){
    g_assert(hsp_param);
    return HSP_Param_create(Match_swap(hsp_param->match),
                            hsp_param->use_horizon);
    }

HSPset *HSPset_create(Sequence *query, Sequence *target,
                      HSP_Param *hsp_param){
    register HSPset *hsp_set = g_new(HSPset, 1);
    g_assert(query);
    g_assert(target);
    hsp_set->ref_count = 1;
    hsp_set->query = Sequence_share(query);
    hsp_set->target = Sequence_share(target);
    hsp_set->param = HSP_Param_share(hsp_param);
    /**/
    if(hsp_param->use_horizon){
        hsp_set->horizon = (gint****)Matrix4d_create(
                                   1 + ((hsp_param->seed_repeat > 1)?2:0),
                                   query->len,
                                   hsp_param->match->query->advance,
                                   hsp_param->match->target->advance,
                                   sizeof(gint));
    } else {
        hsp_set->horizon = NULL;
        }
    hsp_set->hsp_list = g_ptr_array_new();
    hsp_set->is_finalised = FALSE;
    hsp_set->param->has = HSPset_ArgumentSet_create(NULL);
    if(hsp_set->param->has->filter_threshold){
        hsp_set->filter = g_new0(PQueue*, query->len);
        hsp_set->pqueue_set = PQueueSet_create();
    } else {
        hsp_set->filter = NULL;
        hsp_set->pqueue_set = NULL;
        }
    hsp_set->is_empty = TRUE;
    return hsp_set;
    }
/* horizon[score(,repeat_count,diag)]
 *        [query_len]
 *        [query_advance]
 *        [target_advance]
 */

/**/

HSPset *HSPset_share(HSPset *hsp_set){
    g_assert(hsp_set);
    hsp_set->ref_count++;
    return hsp_set;
    }

void HSPset_destroy(HSPset *hsp_set){
    register gint i;
    register HSP *hsp;
    g_assert(hsp_set);
    if(--hsp_set->ref_count)
        return;
    HSP_Param_destroy(hsp_set->param);
    if(hsp_set->filter)
        g_free(hsp_set->filter);
    if(hsp_set->pqueue_set)
        PQueueSet_destroy(hsp_set->pqueue_set);
    Sequence_destroy(hsp_set->query);
    Sequence_destroy(hsp_set->target);
    for(i = 0; i < hsp_set->hsp_list->len; i++){
        hsp = hsp_set->hsp_list->pdata[i];
        HSP_destroy(hsp);
        }
    g_ptr_array_free(hsp_set->hsp_list, TRUE);
    if(hsp_set->horizon)
        g_free(hsp_set->horizon);
    g_free(hsp_set);
    return;
    }

void HSPset_swap(HSPset *hsp_set, HSP_Param *hsp_param){
    register Sequence *query;
    register gint i, query_start;
    register HSP *hsp;
    g_assert(hsp_set->ref_count == 1);
    /* Swap query and target */
    query = hsp_set->query;
    hsp_set->query = hsp_set->target;
    hsp_set->target = query;
    /* Switch parameters */
    HSP_Param_destroy(hsp_set->param);
    hsp_set->param = HSP_Param_share(hsp_param);
    /* Swap HSPs coordinates */
    for(i = 0; i < hsp_set->hsp_list->len; i++){
        hsp = hsp_set->hsp_list->pdata[i];
        query_start = hsp->query_start;
        hsp->query_start = hsp->target_start;
        hsp->target_start = query_start;
        }
    return;
    }

void HSPset_revcomp(HSPset *hsp_set){
    register Sequence *rc_query = Sequence_revcomp(hsp_set->query),
                      *rc_target = Sequence_revcomp(hsp_set->target);
    register gint i;
    register HSP *hsp;
    g_assert(hsp_set);
    g_assert(hsp_set->is_finalised);
    g_assert(hsp_set->ref_count == 1);
    /**/
    Sequence_destroy(hsp_set->query);
    Sequence_destroy(hsp_set->target);
    hsp_set->query = rc_query;
    hsp_set->target = rc_target;
    for(i = 0; i < hsp_set->hsp_list->len; i++){
        hsp = hsp_set->hsp_list->pdata[i];
        hsp->query_start = hsp_set->query->len - HSP_query_end(hsp);
        hsp->target_start = hsp_set->target->len - HSP_target_end(hsp);
        }
    return;
    }


static gint HSP_find_cobs(HSP *hsp){
    register gint i, query_pos = hsp->query_start,
                     target_pos = hsp->target_start;
    register Match_Score score = 0;
    /* Find the HSP centre offset by score */
    for(i = 0; i < hsp->length; i++){
        g_assert(HSP_check_positions(hsp->hsp_set,
                                     query_pos, target_pos));
        score += HSP_get_score(hsp, query_pos, target_pos);
        if(score >= (hsp->score>>1))
            break;
        query_pos += HSP_query_advance(hsp);
        target_pos += HSP_target_advance(hsp);
        }
    return i;
    }

/**/
static void HSP_print_info(HSP *hsp){
    g_print("HSP info (%p)\n"
            "    query_start  = [%d]\n"
            "    target_start = [%d]\n"
            "    length       = [%d]\n"
            "    score        = [%d]\n"
            "    cobs         = [%d]\n",
            (gpointer)hsp,
            hsp->query_start,
            hsp->target_start,
            hsp->length,
            hsp->score,
            hsp->cobs);
    return;
    }

typedef struct {
        HSP *hsp;         /* not freed */
       gint  padding;
       gint  max_advance;
       gint  query_display_pad;
       gint  target_display_pad;
    GString *top;
    GString *mid;
    GString *low;
} HSP_Display;

static HSP_Display *HSP_Display_create(HSP *hsp, gint padding){
    register HSP_Display *hd = g_new(HSP_Display, 1);
    register gint approx_length;
    hd->hsp = hsp;
    hd->padding = padding;
    hd->max_advance = MAX(HSP_query_advance(hsp),
                          HSP_target_advance(hsp));
    hd->query_display_pad = (hd->max_advance
                             -HSP_query_advance(hsp))>>1;
    hd->target_display_pad = (hd->max_advance
                              -HSP_target_advance(hsp))>>1;
    approx_length = hd->max_advance * (hsp->length + 2);
    hd->top = g_string_sized_new(approx_length);
    hd->mid = g_string_sized_new(approx_length);
    hd->low = g_string_sized_new(approx_length);
    return hd;
    }

static void HSP_Display_destroy(HSP_Display *hd){
    g_assert(hd);
    g_string_free(hd->top, TRUE);
    g_string_free(hd->mid, TRUE);
    g_string_free(hd->low, TRUE);
    g_free(hd);
    return;
    }

static void HSP_Display_add(HSP_Display *hd,
                            gchar *top, gchar *mid, gchar *low){
    g_assert(hd);
    g_assert(top);
    g_assert(mid);
    g_assert(low);
    g_string_append(hd->top, top);
    g_string_append(hd->mid, mid);
    g_string_append(hd->low, low);
    g_assert(hd->top->len == hd->mid->len);
    g_assert(hd->mid->len == hd->low->len);
    return;
    }

static gchar HSP_Display_get_ruler_char(HSP_Display *hd, gint pos,
                                       gint advance){
    register gint stop;
    register gint pad_length = 3;
    stop = hd->padding * hd->max_advance;
    if(pos >= stop){
        if(pos < (stop+pad_length)){
            return '#';
            }
        stop = ((hd->padding+hd->hsp->length) * hd->max_advance)
             + pad_length;
        if(pos >= stop){
            if(pos < (stop+pad_length)){
                return '#';
                }
            pos -= pad_length;
            }
        pos -= pad_length;
        }
    if((pos/advance) & 1)
        return '=';
    return '-';
    }

static void HSP_Display_print_ruler(HSP_Display *hd, gint width,
                                    gint pos, gboolean is_query){
    register gint i, adv;
    if(is_query){
        if(HSP_target_advance(hd->hsp) == 1)
            return; /* opposite padding */
        adv = HSP_target_advance(hd->hsp);
    } else { /* Is target */
        if(HSP_query_advance(hd->hsp) == 1)
            return; /* opposite padding */
        adv = HSP_query_advance(hd->hsp);
        }
    g_print(" ruler:[");
    for(i = 0; i < width; i++)
        g_print("%c", HSP_Display_get_ruler_char(hd, pos+i, adv));
    g_print("]\n");
    return;
    }

static void HSP_Display_print(HSP_Display *hd){
    register gint pos, pause, width = 50;
    g_assert(hd);
    g_assert(hd->top->len == hd->mid->len);
    g_assert(hd->mid->len == hd->low->len);
    for(pos = 0, pause = hd->top->len-width; pos < pause; pos += width){
        HSP_Display_print_ruler(hd, width, pos, TRUE);
        g_print(" query:[%.*s]\n       [%.*s]\ntarget:[%.*s]\n",
                width, hd->top->str+pos,
                width, hd->mid->str+pos,
                width, hd->low->str+pos);
        HSP_Display_print_ruler(hd, width, pos, FALSE);
        g_print("\n");
        }
    HSP_Display_print_ruler(hd, hd->top->len-pos, pos, TRUE);
    g_print(" query:[%.*s]\n       [%.*s]\ntarget:[%.*s]\n",
            (gint)(hd->top->len-pos), hd->top->str+pos,
            (gint)(hd->mid->len-pos), hd->mid->str+pos,
            (gint)(hd->low->len-pos), hd->low->str+pos);
    HSP_Display_print_ruler(hd, hd->top->len-pos, pos, FALSE);
    g_print("\n");
    return;
    }

static void HSP_Display_insert(HSP_Display *hd, gint position){
    register gint query_pos, target_pos;
    register gboolean is_padding,
                      query_valid = TRUE, target_valid = TRUE;
    register Match *match = hd->hsp->hsp_set->param->match;
    gchar query_str[4] = {0},
          target_str[4] = {0},
          equiv_str[4] = {0};
    g_assert(hd);
    query_pos  = hd->hsp->query_start
               + (HSP_query_advance(hd->hsp) * position);
    target_pos = hd->hsp->target_start
               + (HSP_target_advance(hd->hsp) * position);
    /* If outside HSP, then is_padding */
    is_padding = ((position < 0) || (position >= hd->hsp->length));
    /* If outside seqs, then invalid */
    query_valid = ( (query_pos >= 0)
                  &&((query_pos+HSP_query_advance(hd->hsp))
                      <= hd->hsp->hsp_set->query->len));
    target_valid = ((target_pos >= 0)
          &&((target_pos+HSP_target_advance(hd->hsp))
              <= hd->hsp->hsp_set->target->len));
    /* Get equiv string */
    if(query_valid && target_valid){
        g_assert(HSP_check_positions(hd->hsp->hsp_set,
                                     query_pos, target_pos));
        HSP_get_display(hd->hsp, query_pos, target_pos, equiv_str);
    } else {
        strncpy(equiv_str, "###", hd->max_advance);
        equiv_str[hd->max_advance] = '\0';
        }
    /* Get query string */
    if(query_valid){
        Match_Strand_get_raw(match->query, hd->hsp->hsp_set->query,
                             query_pos, query_str);
        if((match->query->advance == 1) && (hd->max_advance == 3)){
            query_str[1] = query_str[0];
            query_str[0] = query_str[2] = ' ';
            query_str[3] = '\0';
            }
    } else {
        strncpy(query_str, "###", hd->max_advance);
        query_str[hd->max_advance] = '\0';
        }
    /* Get target string */
    if(target_valid){
        Match_Strand_get_raw(match->target, hd->hsp->hsp_set->target,
                             target_pos, target_str);
        if((match->target->advance == 1) && (hd->max_advance == 3)){
            target_str[1] = target_str[0];
            target_str[0] = target_str[2] = ' ';
            target_str[3] = '\0';
            }
    } else {
        strncpy(target_str, "###", hd->max_advance);
        target_str[hd->max_advance] = '\0';
        }
    /* Make lower case for padding */
    if(is_padding){
        strdown(query_str);
        strdown(target_str);
    } else {
        strup(query_str);
        strup(target_str);
        }
    HSP_Display_add(hd, query_str, equiv_str, target_str);
    return;
    }

static void HSP_print_alignment(HSP *hsp){
    register HSP_Display *hd = HSP_Display_create(hsp, 10);
    register gint i;
    for(i = 0; i < hd->padding; i++) /* Pre-padding */
        HSP_Display_insert(hd, i-hd->padding);
    /* Use pad_length == 3 */
    HSP_Display_add(hd, " < ", " < ", " < "); /* Start divider */
    for(i = 0; i < hsp->length; i++) /* The HSP itself */
        HSP_Display_insert(hd, i);
    HSP_Display_add(hd, " > ", " > ", " > "); /* End divider */
    for(i = 0; i < hd->padding; i++) /* Post-padding */
        HSP_Display_insert(hd, hsp->length+i);
    HSP_Display_print(hd);
    HSP_Display_destroy(hd);
    return;
    }
/*
 * HSP display style:
 *
 *  =-=-=-   =-=-=-=-=-=-=-=-=-   =-=-=-
 *  nnnnnn < ACGACGCCCACGATCGAT > nnn###
 *     ||| < |||:::|||   |||||| > |||
 *  ### x  <  A  R  N  D  C  Q  >  x ###
 *  ===---   ===---===---===---   ===---
 */

static void HSP_print_sugar(HSP *hsp){
    g_print("sugar: %s %d %d %c %s %d %d %c %d\n",
            hsp->hsp_set->query->id,
            hsp->query_start,
            hsp->length*HSP_query_advance(hsp),
            Sequence_get_strand_as_char(hsp->hsp_set->query),
            hsp->hsp_set->target->id,
            hsp->target_start,
            hsp->length*HSP_target_advance(hsp),
            Sequence_get_strand_as_char(hsp->hsp_set->target),
            hsp->score);
    return;
    }
/* Sugar output format:
 * sugar: <query_id> <query_start> <query_length> <query_strand>
 *        <target_id> <target_start> <target_start> <target_strand>
 *        <score>
 */

void HSP_print(HSP *hsp, gchar *name){
    g_print("draw_hsp(%d, %d, %d, %d, %d, %d, \"%s\")\n",
            hsp->query_start,
            hsp->target_start,
            hsp->length,
            hsp->cobs,
            HSP_query_advance(hsp),
            HSP_target_advance(hsp),
            name);
    HSP_print_info(hsp);
    HSP_print_alignment(hsp);
    HSP_print_sugar(hsp);
    return;
    }

void HSPset_print(HSPset *hsp_set){
    register gint i;
    register gchar *name;
    g_print("HSPset [%p] contains [%d] hsps\n",
            (gpointer)hsp_set, hsp_set->hsp_list->len);
    g_print("Comparison of [%s] and [%s]\n",
            hsp_set->query->id, hsp_set->target->id);
    for(i = 0; i < hsp_set->hsp_list->len; i++){
        name = g_strdup_printf("hsp [%d]", i+1);
        HSP_print(hsp_set->hsp_list->pdata[i], name);
        g_free(name);
        }
    return;
    }

/**/

static void HSP_init(HSP *nh){
    register gint i;
    register gint query_pos, target_pos;
    g_assert(HSP_check(nh));
    /* Initial hsp score */
    query_pos = nh->query_start;
    target_pos = nh->target_start;
    nh->score = 0;
    for(i = 0; i < nh->length; i++){
        g_assert(HSP_check_positions(nh->hsp_set,
                                     query_pos, target_pos));
        nh->score += HSP_get_score(nh, query_pos, target_pos);
        query_pos += HSP_query_advance(nh);
        target_pos += HSP_target_advance(nh);
        }
    if(nh->score < 0){
        HSP_print(nh, "Bad HSP seed");
        g_error("Initial HSP score [%d] less than zero", nh->score);
        }
    g_assert(HSP_check(nh));
    return;
    }

static void HSP_extend(HSP *nh, gboolean forbid_masked){
    register Match_Score score, maxscore;
    register gint query_pos, target_pos;
    register gint extend, maxext;
    g_assert(HSP_check(nh));
    /* extend left */
    maxscore = score = nh->score;
    query_pos = nh->query_start-HSP_query_advance(nh);
    target_pos = nh->target_start-HSP_target_advance(nh);
    for(extend = 1, maxext = 0;
       ((query_pos >= 0) && (target_pos >= 0));
       extend++){
        g_assert(HSP_check_positions(nh->hsp_set,
                                     query_pos, target_pos));
        if((forbid_masked)
        && (HSP_query_masked(nh, query_pos)
         || HSP_target_masked(nh, target_pos)))
            break;
        score += HSP_get_score(nh, query_pos, target_pos);
        if(maxscore <= score){
            maxscore = score;
            maxext   = extend;
        } else {
            if(score < 0) /* See note below */
                break;
            if((maxscore-score) >= nh->hsp_set->param->dropoff)
                break;
            }
        query_pos -= HSP_query_advance(nh);
        target_pos -= HSP_target_advance(nh);
        }
    query_pos = HSP_query_end(nh);
    target_pos = HSP_target_end(nh);
    nh->query_start -= (maxext * HSP_query_advance(nh));
    nh->target_start -= (maxext * HSP_target_advance(nh));
    nh->length += maxext;
    score = maxscore;
    /* extend right */
    for(extend = 1, maxext = 0;
       (  ((query_pos+HSP_query_advance(nh))
           <= nh->hsp_set->query->len)
       && ((target_pos+HSP_target_advance(nh))
           <= nh->hsp_set->target->len) );
       extend++){
        g_assert(HSP_check_positions(nh->hsp_set,
                                     query_pos, target_pos));
        if((forbid_masked)
        && (HSP_query_masked(nh, query_pos)
         || HSP_target_masked(nh, target_pos)))
            break;
        score += HSP_get_score(nh, query_pos, target_pos);
        if(maxscore <= score){
            maxscore = score;
            maxext = extend;
        } else {
            if(score < 0) /* See note below */
                break;
            if((maxscore-score) >= nh->hsp_set->param->dropoff)
                break;
            }
        query_pos += HSP_query_advance(nh);
        target_pos += HSP_target_advance(nh);
        }
    nh->score = maxscore;
    nh->length += maxext;
    g_assert(HSP_check(nh));
    return;
    }
/* The score cannot be allowed to drop below zero in the HSP,
 * as this can result in overlapping HSPs in some circumstances.
 */

static HSP *HSP_create(HSP *nh){
    register HSP *hsp;
#ifdef USE_PTHREADS
    pthread_mutex_lock(&nh->hsp_set->param->hsp_recycle_lock);
#endif /* USE_PTHREADS */
    hsp = RecycleBin_alloc(nh->hsp_set->param->hsp_recycle);
#ifdef USE_PTHREADS
    pthread_mutex_unlock(&nh->hsp_set->param->hsp_recycle_lock);
#endif /* USE_PTHREADS */
    hsp->hsp_set = nh->hsp_set;
    hsp->query_start = nh->query_start;
    hsp->target_start = nh->target_start;
    hsp->length = nh->length;
    hsp->score = nh->score;
    hsp->cobs = nh->cobs; /* Value can be set by HSPset_finalise(); */
    return hsp;
    }

void HSP_destroy(HSP *hsp){
    register HSPset *hsp_set = hsp->hsp_set;
#ifdef USE_PTHREADS
    pthread_mutex_lock(&hsp_set->param->hsp_recycle_lock);
#endif /* USE_PTHREADS */
    RecycleBin_recycle(hsp_set->param->hsp_recycle, hsp);
#ifdef USE_PTHREADS
    pthread_mutex_unlock(&hsp_set->param->hsp_recycle_lock);
#endif /* USE_PTHREADS */
    return;
    }

static void HSP_trim_ends(HSP *hsp){
    register gint i;
    register gint query_pos, target_pos;
    /* Trim left to first good match */
    g_assert(HSP_check(hsp));
    for(i = 0; i < hsp->length; i++){
        if(HSP_get_score(hsp, hsp->query_start, hsp->target_start) > 0)
            break;
        hsp->query_start += HSP_query_advance(hsp);
        hsp->target_start += HSP_target_advance(hsp);
        }
    hsp->length -= i;
    /**/
    g_assert(HSP_check(hsp));
    query_pos = HSP_query_end(hsp) - HSP_query_advance(hsp);
    target_pos = HSP_target_end(hsp) - HSP_target_advance(hsp);
    /* Trim right to last good match */
    while(hsp->length > 0){
        g_assert(HSP_check_positions(hsp->hsp_set,
                                     query_pos, target_pos));
        if(HSP_get_score(hsp, query_pos, target_pos) > 0)
            break;
        hsp->length--;
        query_pos -= HSP_query_advance(hsp);
        target_pos -= HSP_target_advance(hsp);
        }
    g_assert(HSP_check(hsp));
    return;
    }
/* This is to remove any unmatching ends from the HSP seed.
 */

static gboolean HSP_PQueue_comp_func(gpointer low, gpointer high,
                                     gpointer user_data){
    register HSP *hsp_low = (HSP*)low, *hsp_high = (HSP*)high;
    return hsp_low->score - hsp_high->score;
    }

static void HSP_store(HSP *nascent_hsp){
    register HSPset *hsp_set = nascent_hsp->hsp_set;
    register PQueue *pq;
    register HSP *hsp = NULL;
    g_assert(nascent_hsp);
    if(nascent_hsp->score < hsp_set->param->threshold)
        return;
    if(hsp_set->param->has->filter_threshold){ /* If have filter */
        /* Get cobs value */
        nascent_hsp->cobs = HSP_find_cobs(nascent_hsp);
        pq = hsp_set->filter[HSP_query_cobs(nascent_hsp)];
        if(pq){ /* Put in PQueue if better than worst */
            if(PQueue_total(pq)
               < hsp_set->param->has->filter_threshold){
                hsp = HSP_create(nascent_hsp);
                PQueue_push(pq, hsp);
            } else {
                g_assert(PQueue_total(pq));
                hsp = PQueue_top(pq);
                if(hsp->score < nascent_hsp->score){
                    hsp = PQueue_pop(pq);
                    HSP_destroy(hsp);
                    hsp = HSP_create(nascent_hsp);
                    PQueue_push(pq, hsp);
                    }
                }
        } else {
            pq = PQueue_create(hsp_set->pqueue_set,
                               HSP_PQueue_comp_func, NULL);
            hsp_set->filter[HSP_query_cobs(nascent_hsp)] = pq;
            hsp = HSP_create(nascent_hsp);
            PQueue_push(pq, hsp);
            hsp_set->is_empty = FALSE;
            }
    } else {
        hsp = HSP_create(nascent_hsp);
        g_ptr_array_add(hsp_set->hsp_list, hsp);
        hsp_set->is_empty = FALSE;
        }
    return;
    }
/* FIXME: optimisation: could store HSPs as a list up until
 *                      filter_threshold, then convert to a PQueue
 */

void HSPset_seed_hsp(HSPset *hsp_set,
                     guint query_start, guint target_start){
    register gint diag_pos
        = ((target_start * hsp_set->param->match->query->advance)
          -(query_start * hsp_set->param->match->target->advance));
    register gint query_frame = query_start
                              % hsp_set->param->match->query->advance,
                  target_frame = target_start
                               % hsp_set->param->match->target->advance;
    register gint section_pos = (diag_pos + hsp_set->query->len)
                              % hsp_set->query->len;
    HSP nascent_hsp;
    g_assert(!hsp_set->is_finalised);
    /**/
    g_assert(section_pos >= 0);
    g_assert(section_pos < hsp_set->query->len);
    /* Clear if diag has changed */
    if(hsp_set->param->seed_repeat > 1){
        if(hsp_set->horizon[2][section_pos][query_frame][target_frame]
               != (diag_pos + hsp_set->query->len)){
            hsp_set->horizon[0][section_pos][query_frame][target_frame] = 0;
            hsp_set->horizon[1][section_pos][query_frame][target_frame] = 0;
            hsp_set->horizon[2][section_pos][query_frame][target_frame]
                = diag_pos + hsp_set->query->len;
            }
        }
    /* Check whether we have seen this HSP already */
    if(target_start < hsp_set->horizon[0]
                                      [section_pos]
                                      [query_frame]
                                      [target_frame])
        return;
    if(hsp_set->param->seed_repeat > 1){
        if(++hsp_set->horizon[1][section_pos][query_frame][target_frame]
                < hsp_set->param->seed_repeat)
            return;
        hsp_set->horizon[1][section_pos][query_frame][target_frame] = 0;
        }
    /* Nascent HSP building: */
    nascent_hsp.hsp_set      = hsp_set;
    nascent_hsp.query_start  = query_start;
    nascent_hsp.target_start = target_start;
    nascent_hsp.length       = hsp_set->param->seedlen;
    nascent_hsp.cobs         = 0;
    g_assert(HSP_check(&nascent_hsp));
    HSP_trim_ends(&nascent_hsp);
    /* Score is irrelevant before HSP_init() */
    HSP_init(&nascent_hsp);
    /* Try to make above threshold HSP using masking */
    if(hsp_set->param->match->query->mask_func
    || hsp_set->param->match->target->mask_func){
        HSP_extend(&nascent_hsp, TRUE);
        if(nascent_hsp.score < hsp_set->param->threshold){
            hsp_set->horizon[0][section_pos][query_frame][target_frame]
                = HSP_target_end(&nascent_hsp);
            return;
            }
        }
    /* Extend the HSP again ignoring masking */
    HSP_extend(&nascent_hsp, FALSE);
    HSP_store(&nascent_hsp);
    hsp_set->horizon[0][section_pos][query_frame][target_frame]
            = HSP_target_end(&nascent_hsp);
    return;
    }

void HSPset_add_known_hsp(HSPset *hsp_set,
                          guint query_start, guint target_start,
                          guint length){
    HSP nascent_hsp;
    nascent_hsp.hsp_set      = hsp_set;
    nascent_hsp.query_start  = query_start;
    nascent_hsp.target_start = target_start;
    nascent_hsp.length       = length;
    nascent_hsp.cobs         = 0;
    /* Score is irrelevant before HSP_init() */
    HSP_init(&nascent_hsp);
    HSP_store(&nascent_hsp);
    return;
    }

static void HSPset_seed_hsp_sorted(HSPset *hsp_set,
                     guint query_start, guint target_start,
                     gint ***horizon){
    HSP nascent_hsp;
    register gint diag_pos
        = ((target_start * hsp_set->param->match->query->advance)
          -(query_start * hsp_set->param->match->target->advance));
    register gint query_frame = query_start
                              % hsp_set->param->match->query->advance,
                  target_frame = target_start
                               % hsp_set->param->match->target->advance;
    register gint section_pos = (diag_pos + hsp_set->query->len)
                              % hsp_set->query->len;
    g_assert(!hsp_set->is_finalised);
    g_assert(!hsp_set->horizon);
    g_assert(section_pos >= 0);
    g_assert(section_pos < hsp_set->query->len);
    /**/
    if(horizon[1][query_frame][target_frame] != section_pos){
        horizon[1][query_frame][target_frame] = section_pos;
        horizon[0][query_frame][target_frame] = 0;
        horizon[2][query_frame][target_frame] = 0;
        }
/* FIXME: seedrepeat overflow here */
    if(++horizon[2][query_frame][target_frame] < hsp_set->param->seed_repeat)
        return;
    horizon[2][query_frame][target_frame] = 0;
    /* Check whether we have seen this HSP already */
    if(target_start < horizon[0][query_frame][target_frame])
        return;
    /**/
    /* Nascent HSP building: */
    nascent_hsp.hsp_set      = hsp_set;
    nascent_hsp.query_start  = query_start;
    nascent_hsp.target_start = target_start;
    nascent_hsp.length       = hsp_set->param->seedlen;
    nascent_hsp.cobs         = 0;
    g_assert(HSP_check(&nascent_hsp));
    HSP_trim_ends(&nascent_hsp);
    /* Score is irrelevant before HSP_init() */
    HSP_init(&nascent_hsp);
    /* Try to make above threshold HSP using masking */
    if(hsp_set->param->match->query->mask_func
    || hsp_set->param->match->target->mask_func){
        HSP_extend(&nascent_hsp, TRUE);
        if(nascent_hsp.score < hsp_set->param->threshold){
            horizon[0][query_frame][target_frame]
                = HSP_target_end(&nascent_hsp);
            return;
            }
        }
    /* Extend the HSP again ignoring masking */
    HSP_extend(&nascent_hsp, FALSE);
    HSP_store(&nascent_hsp);
    /**/
    horizon[0][query_frame][target_frame] = HSP_target_end(&nascent_hsp);
    return;
    }
/* horizon[0] = diag
 * horizon[1] = target_end
 * horizon[2] = repeat_count
 */

/* Need to use the global to pass q,t advance to qsort compare func */
static HSPset *HSPset_seed_compare_hsp_set = NULL;

static int HSPset_seed_compare(const void *a, const void *b){
    register guint *seed_a = (guint*)a,
                   *seed_b = (guint*)b;
    register gint diag_a, diag_b;
    register HSPset *hsp_set = HSPset_seed_compare_hsp_set;
    g_assert(hsp_set);
    diag_a = ((seed_a[1] * hsp_set->param->match->query->advance)
           -  (seed_a[0] * hsp_set->param->match->target->advance)),
    diag_b = ((seed_b[1] * hsp_set->param->match->query->advance)
           -  (seed_b[0] * hsp_set->param->match->target->advance));
    if(diag_a == diag_b)
        return seed_a[0] - seed_b[0];
    return diag_a - diag_b;
    }

void HSPset_seed_all_hsps(HSPset *hsp_set,
                          guint *seed_list, guint seed_list_len){
    register gint i;
    register gint ***horizon;
    register gint qpos, tpos;
    if(seed_list_len > 1){
        HSPset_seed_compare_hsp_set = hsp_set;
        qsort(seed_list, seed_list_len, sizeof(guint) << 1,
              HSPset_seed_compare);
        HSPset_seed_compare_hsp_set = NULL;
        }
    if(seed_list_len){
        horizon = (gint***)Matrix3d_create(3,
                     hsp_set->param->match->query->advance,
                     hsp_set->param->match->target->advance,
                     sizeof(gint));
        for(i = 0; i < seed_list_len; i++){
            HSPset_seed_hsp_sorted(hsp_set,
                                   seed_list[(i << 1)],
                                   seed_list[(i << 1) + 1],
                                   horizon);
            qpos = seed_list[(i << 1)];
            tpos = seed_list[(i << 1) + 1];
            }
        g_free(horizon);
        }
    HSPset_finalise(hsp_set);
    return;
    }

/**/

HSPset *HSPset_finalise(HSPset *hsp_set){
    register gint i;
    register HSP *hsp;
    register PQueue *pq;
    g_assert(!hsp_set->is_finalised);
    hsp_set->is_finalised = TRUE;
    if(hsp_set->param->has->filter_threshold && (!hsp_set->is_empty)){
        /* Get HSPs from each PQueue */
        for(i = 0; i < hsp_set->query->len; i++){
            pq = hsp_set->filter[i];
            if(pq){
                while(PQueue_total(pq)){
                    hsp = PQueue_pop(pq);
                    g_ptr_array_add(hsp_set->hsp_list, hsp);
                    }
                }
            }
        }
    /* Set cobs for each HSP */
    if(!hsp_set->param->has->filter_threshold){
        for(i = 0; i < hsp_set->hsp_list->len; i++){
            hsp = hsp_set->hsp_list->pdata[i];
            hsp->cobs = HSP_find_cobs(hsp);
            }
        }
    hsp_set->is_finalised = TRUE;
    return hsp_set;
    }

/**/

static int HSPset_sort_by_diag_then_query_start(const void *a,
                                                const void *b){
    register HSP **hsp_a = (HSP**)a, **hsp_b = (HSP**)b;
    register gint diag_a = HSP_diagonal(*hsp_a),
                  diag_b = HSP_diagonal(*hsp_b);
    if(diag_a == diag_b)
        return (*hsp_a)->query_start - (*hsp_b)->query_start;
    return diag_a - diag_b;
    }

static Match_Score HSP_score_overlap(HSP *left, HSP *right){
    register Match_Score score = 0;
    register gint query_pos, target_pos;
    g_assert(left->hsp_set == right->hsp_set);
    g_assert(HSP_diagonal(left) == HSP_diagonal(right));
    query_pos = HSP_query_end(left) - HSP_query_advance(left);
    target_pos = HSP_target_end(left) - HSP_target_advance(left);
    while(query_pos >= right->query_start){
        score += HSP_get_score(left, query_pos, target_pos);
        query_pos -= HSP_query_advance(left);
        target_pos -= HSP_target_advance(left);
        }
    query_pos = right->query_start;
    target_pos = right->target_start;
    while(query_pos < (HSP_query_end(left)- HSP_query_advance(right))){
        score += HSP_get_score(right, query_pos, target_pos);
        query_pos += HSP_query_advance(right);
        target_pos += HSP_target_advance(right);
        }
    return score;
    }
/* Returns score for overlapping region of HSPs on same diagonal */

void HSPset_filter_ungapped(HSPset *hsp_set){
    register GPtrArray *new_hsp_list;
    register HSP *curr_hsp, *prev_hsp;
    register gboolean del_prev, del_curr;
    register gint i;
    register Match_Score score;
    /* Filter strongly overlapping HSPs on same diagonal
     * but different frames (happens with 3:3 HSPs only)
     */
    if((hsp_set->hsp_list->len > 1)
    && (hsp_set->param->match->query->advance == 3)
    && (hsp_set->param->match->target->advance == 3)){
        /* FIXME: should not sort when using all-at-once HSPset */
        qsort(hsp_set->hsp_list->pdata, hsp_set->hsp_list->len,
              sizeof(gpointer), HSPset_sort_by_diag_then_query_start);
        prev_hsp = hsp_set->hsp_list->pdata[0];
        del_prev = FALSE;
        del_curr = FALSE;
        new_hsp_list = g_ptr_array_new();
        for(i = 1; i < hsp_set->hsp_list->len; i++){
            curr_hsp = hsp_set->hsp_list->pdata[i];
            del_curr = FALSE;
            if((HSP_diagonal(prev_hsp) == HSP_diagonal(curr_hsp))
            && (HSP_query_end(prev_hsp) > curr_hsp->query_start)){
                score = HSP_score_overlap(prev_hsp, curr_hsp);
                if((score << 1) > (curr_hsp->score + prev_hsp->score)){
                    /* FIXME: use codon_usage scores here instead */
                    if(prev_hsp->score < curr_hsp->score){
                        del_prev = TRUE;
                    } else {
                        del_curr = TRUE;
                        }
                    }
                }
            if(del_prev)
                HSP_destroy(prev_hsp);
            else
                g_ptr_array_add(new_hsp_list, prev_hsp);
            prev_hsp = curr_hsp;
            del_prev = del_curr;
            }
        if(del_prev)
            HSP_destroy(prev_hsp);
        else
            g_ptr_array_add(new_hsp_list, prev_hsp);
        g_ptr_array_free(hsp_set->hsp_list, TRUE);
        hsp_set->hsp_list = new_hsp_list;
        }
    return;
    }

/**/

#define HSPset_SList_PAGE_BIT_WIDTH 10
#define HSPset_SList_PAGE_SIZE (1 << HSPset_SList_PAGE_BIT_WIDTH)

RecycleBin *HSPset_SList_RecycleBin_create(void){
    return RecycleBin_create("HSPset_Slist", sizeof(HSPset_SList_Node),
                             4096);
    }

HSPset_SList_Node *HSPset_SList_append(RecycleBin *recycle_bin,
                                       HSPset_SList_Node *next,
                                       gint query_pos, gint target_pos){
    register HSPset_SList_Node *node = RecycleBin_alloc(recycle_bin);
    node->next = next;
    node->query_pos = query_pos;
    node->target_pos = target_pos;
    return node;
    }

#if 0
typedef struct {
                 gint      page_alloc;
                 gint      page_total;
    HSPset_SList_Node    **diag_page_list;
                 gint  ****horizon;
                 gint     *page_used;
                 gint      page_used_total;
} HSPset_SeedData;

static HSPset_SeedData *HSPset_SeedData_create(HSP_Param *hsp_param,
                                               gint target_len){
    register HSPset_SeedData *seed_data = g_new(HSPset_SeedData, 1);
    seed_data->page_total = (target_len
                          >> HSPset_SList_PAGE_BIT_WIDTH) + 1;
    seed_data->page_alloc = seed_data->page_total;
    seed_data->diag_page_list = g_new0(HSPset_SList_Node*, seed_data->page_total);
    seed_data->page_used = g_new(gint, seed_data->page_total);
    seed_data->horizon = (gint****)Matrix4d_create(
                           2 + ((hsp_param->seed_repeat > 1)?1:0),
                           HSPset_SList_PAGE_SIZE,
                           hsp_param->match->query->advance,
                           hsp_param->match->target->advance,
                           sizeof(gint));
    seed_data->page_used_total = 0;
    return seed_data;
    }

static void HSPset_SeedData_destroy(HSPset_SeedData *seed_data){
    g_free(seed_data->diag_page_list);
    g_free(seed_data->page_used);
    g_free(seed_data->horizon);
    g_free(seed_data);
    return;
    }

static void HSPset_SeedData_set_target_len(HSPset_SeedData *seed_data,
                                           HSPset *hsp_set){
    register gint new_page_total = (hsp_set->target->len
                                    >> HSPset_SList_PAGE_BIT_WIDTH) + 1;
    register gint i, a, b, c, d;
    seed_data->page_total = new_page_total;
    if(seed_data->page_alloc < new_page_total){
        seed_data->page_alloc = seed_data->page_total;
        g_free(seed_data->diag_page_list);
        seed_data->diag_page_list = g_new(HSPset_SList_Node*,
                                          seed_data->page_total);
        g_free(seed_data->page_used);
        seed_data->page_used = g_new(gint, seed_data->page_total);
        }
    /* Clear diag_page_list */
    for(i = 0; i < seed_data->page_total; i++)
        seed_data->diag_page_list[i] = 0;
    /* Clear horizon */
    for(a = 2 + ((hsp_set->param->seed_repeat > 1)?1:0); a >= 0; a--)
        for(b = HSPset_SList_PAGE_SIZE; b >= 0; b--)
            for(c = hsp_set->param->match->query->advance; c >= 0; c--)
                for(d = hsp_set->param->match->target->advance; d >= 0; d--)
                    seed_data->horizon[a][b][c][d] = 0;
    seed_data->page_used_total = 0;
    return;
    }
#endif /* 0 */

void HSPset_seed_all_qy_sorted(HSPset *hsp_set, HSPset_SList_Node *seed_list){
    register gint page_total = (hsp_set->target->len
                                >> HSPset_SList_PAGE_BIT_WIDTH) + 1;
    register HSPset_SList_Node **diag_page_list
        = g_new0(HSPset_SList_Node*, page_total);
    register gint *page_used = g_new(gint, page_total);
    register gint ****horizon = (gint****)Matrix4d_create(
                           2 + ((hsp_set->param->seed_repeat > 1)?1:0),
                           HSPset_SList_PAGE_SIZE,
                           hsp_set->param->match->query->advance,
                           hsp_set->param->match->target->advance,
                           sizeof(gint));
    /*
    register HSPset_SeedData *seed_data = HSPset_SeedData_create(hsp_set->param,
                                                           hsp_set->target->len);
    */
    register gint i, page, diag_pos, query_frame, target_frame,
                   section_pos, page_pos;
    register HSPset_SList_Node *seed;
    register gint page_used_total = 0;
    HSP nascent_hsp;
    /* g_message("[%s] with [%d]", __FUNCTION__, hsp_set->target->len); */
    /* HSPset_SeedData_set_target_len(seed_data, hsp_set); */
    /* Bin on diagonal into pages */
    while(seed_list){
        seed = seed_list;
        seed_list = seed_list->next;
        /**/
        diag_pos = (seed->target_pos
                    * hsp_set->param->match->query->advance)
                 - (seed->query_pos
                    * hsp_set->param->match->target->advance);
        section_pos = ((diag_pos + hsp_set->target->len)
                    % hsp_set->target->len);
        page = (section_pos >> HSPset_SList_PAGE_BIT_WIDTH);
        g_assert(section_pos >= 0);
        g_assert(section_pos < hsp_set->target->len);
        g_assert(page >= 0);
        g_assert(page < page_total);
        /**/
        if(!diag_page_list[page])
            page_used[page_used_total++] = page;
        seed->next = diag_page_list[page];
        diag_page_list[page] = seed;
        }
    /* Seed each page using page horizon */
    for(i = 0; i < page_used_total; i++){
        page = page_used[i];
        for(seed = diag_page_list[page]; seed; seed = seed->next){
            g_assert((!seed->next)
                   || (seed->query_pos <= seed->next->query_pos));
            diag_pos = (seed->target_pos
                        * hsp_set->param->match->query->advance)
                     - (seed->query_pos
                        * hsp_set->param->match->target->advance);
            query_frame  = seed->query_pos
                         % hsp_set->param->match->query->advance;
            target_frame  = seed->target_pos
                          % hsp_set->param->match->target->advance;
            section_pos = ((diag_pos + hsp_set->target->len)
                        % hsp_set->target->len);
            page_pos = section_pos
                     - (page << HSPset_SList_PAGE_BIT_WIDTH);
            g_assert(page_pos >= 0);
            g_assert(page_pos < HSPset_SList_PAGE_SIZE);
            /* Clear if page has changed */
            if(horizon[1][page_pos][query_frame][target_frame] != page){
               horizon[0][page_pos][query_frame][target_frame] = 0;
               horizon[1][page_pos][query_frame][target_frame] = page;
               if(hsp_set->param->seed_repeat > 1)
                   horizon[2][page_pos][query_frame][target_frame] = 0;
               }
            if(seed->query_pos < horizon[0][page_pos][query_frame][target_frame])
                continue;
            if(hsp_set->param->seed_repeat > 1){
                if(++horizon[2][page_pos][query_frame][target_frame]
                        < hsp_set->param->seed_repeat){
                    continue;
                    }
                horizon[2][page_pos][query_frame][target_frame] = 0;
                }
            /* Nascent HSP building: */
            nascent_hsp.hsp_set      = hsp_set;
            nascent_hsp.query_start  = seed->query_pos;
            nascent_hsp.target_start = seed->target_pos;
            nascent_hsp.length       = hsp_set->param->seedlen;
            nascent_hsp.cobs         = 0;
            g_assert(HSP_check(&nascent_hsp));
            HSP_trim_ends(&nascent_hsp);
            /* Score is irrelevant before HSP_init() */
            HSP_init(&nascent_hsp);
            /* Try to make above threshold HSP using masking */
            if(hsp_set->param->match->query->mask_func
            || hsp_set->param->match->target->mask_func){
                HSP_extend(&nascent_hsp, TRUE);
                if(nascent_hsp.score < hsp_set->param->threshold){
                    horizon[0][page_pos][query_frame][target_frame]
                        = HSP_query_end(&nascent_hsp);
                    continue;
                    }
                }
            /* Extend the HSP again ignoring masking */
            HSP_extend(&nascent_hsp, FALSE);
            HSP_store(&nascent_hsp);
            /**/
            horizon[0][page_pos][query_frame][target_frame]
                    = HSP_query_end(&nascent_hsp);
            }
        }
    g_free(diag_page_list);
    g_free(page_used);
    g_free(horizon);
    HSPset_finalise(hsp_set);
    /* HSPset_SeedData_destroy(seed_data); */
    return;
    }
/* horizon[horizon][mailbox][seed_repeat]
 *        [page_size][qadv][tadv]
 */

/**/

