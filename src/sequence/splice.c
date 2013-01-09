/****************************************************************\
*                                                                *
*  Library for Splice site prediction.                           *
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
#include <math.h>   /* For log() */
#include <ctype.h>  /* For toupper() */
#include <string.h> /* For strlen() */
#include <strings.h> /* For strcasecmp() */

#include "splice.h"
#include "matrix.h"
#include "lineparse.h"

#ifndef SPLICE_IMPOSSIBLY_LOW_SCORE
#define SPLICE_IMPOSSIBLY_LOW_SCORE -987654321.0
#endif /* SPLICE_IMPOSSIBLY_LOW_SCORE */

Splice_ArgumentSet *Splice_ArgumentSet_create(Argument *arg){
    register ArgumentSet *as;
    static Splice_ArgumentSet sas = {NULL, NULL, FALSE};
    if(arg){
        as = ArgumentSet_create("Splice Site Prediction Options");
        ArgumentSet_add_option(as, 0, "splice3", NULL,
                "Supply frequency matrix for 3' splice sites", "primate",
                Argument_parse_string, &sas.splice3_data_path);
        ArgumentSet_add_option(as, 0, "splice5", NULL,
                "Supply frequency matrix for 5' splice sites", "primate",
                Argument_parse_string, &sas.splice5_data_path);
        ArgumentSet_add_option(as, 0, "forcegtag", NULL,
                "Force use of gt...ag splice sites", "FALSE",
                Argument_parse_boolean, &sas.force_gtag);
        Argument_absorb_ArgumentSet(arg, as);
        }
    return &sas;
    }

/* ------------------------------------------------------------ */
/* --- Start of Splice data ----------------------------------- */

/* Splice site data from:
 *
 * "Splice Junctions, Branch Point Sites, and Exons:
 *  Sequence Statistics, Identification, and Applications
 *  to Genome Project"
 *  Periannan Senapathy, Marvin B. Shapiro, and Nomi L. Harris
 *  Methods in Enzymology, Vol 183, pp 252-278.
 */

static void SplicePredictor_add_data_primate_5SS(SplicePredictor *sp){
    register gint i, j;
    gint splice_pssm_primate_5SS_data[9][4] = {
    /*            A    C    G   T                        */
    /*  -3 */ {  28,  40,  17,  14 },
    /*  -2 */ {  59,  14,  13,  14 },
    /*  -1 */ {   8,   5,  81,   6 },
                                      /* <-- Splice site */
    /*   1 */ {   0,   0, 100,   0 }, /* G               */
    /*   2 */ {   0,   0,   0, 100 }, /* T               */
    /*   3 */ {  54,   2,  42,   2 },
    /*   4 */ {  74,   8,  11,   8 },
    /*   5 */ {   5,   6,  85,   4 },
    /*   6 */ {  16,  18,  21,  45 }
    };
    sp->model_length = 9;
    sp->model_splice_after = 3;
    sp->model_data = (gfloat**)Matrix2d_create(sp->model_length,
                                               5, sizeof(gfloat));
    for(i = 0; i < sp->model_length; i++)
        for(j = 0; j < 4; j++)
            sp->model_data[i][j]
            = (gfloat)splice_pssm_primate_5SS_data[i][j];
    return;
    }

static void SplicePredictor_add_data_primate_3SS(SplicePredictor *sp){
    register gint i, j;
    gint splice_pssm_primate_3SS_data[15][4]  = {
    /*            A     C     G     T  */
    /* -14 */ {  10,  31,  14,  44 },
    /* -13 */ {   8,  36,  14,  43 },
    /* -12 */ {   6,  34,  12,  48 },
    /* -11 */ {   6,  34,   8,  52 },
    /* -10 */ {   9,  37,   9,  45 },
    /*  -9 */ {   9,  38,  10,  44 },
    /*  -8 */ {   8,  44,   9,  40 },
    /*  -7 */ {   9,  41,   8,  41 },
    /*  -6 */ {   6,  44,   6,  45 },
    /*  -5 */ {   6,  40,   6,  48 },
    /*  -4 */ {  23,  28,  26,  23 },
    /*  -3 */ {   2,  79,   1,  18 },
    /*  -2 */ { 100,   0,   0,   0 }, /* A               */
    /*  -1 */ {   0,   0, 100,   0 }, /* G               */
                                      /* <-- Splice site */
    /*   1 */ {  28,  14,  47,  11 }
    };
    sp->model_length = 15;
    sp->model_splice_after = 14;
    sp->model_data = (gfloat**)Matrix2d_create(sp->model_length,
                                               5, sizeof(gfloat));
    for(i = 0; i < sp->model_length; i++)
        for(j = 0; j < 4; j++)
            sp->model_data[i][j]
            = (gfloat)splice_pssm_primate_3SS_data[i][j];
    return;
    }

/* --- End of Splice data ------------------------------------- */
/* ------------------------------------------------------------ */

static void SplicePredictor_parse_data(SplicePredictor *sp, gchar *path){
    register FILE *fp = fopen(path, "r");
    register LineParse *lp;
    register gint i;
    gfloat num;
    register gchar *word;
    register GArray *number_list = g_array_new(FALSE, FALSE, sizeof(gfloat));
    if(!fp)
        g_error("Could not open splice frequency data [%s]", path);
    lp = LineParse_create(fp);
    /**/
    sp->model_length = 0;
    while(LineParse_word(lp) != EOF){
        if(lp->word->len){
            word = lp->word->pdata[0];
            if(word[0] == '#'){
                continue;
                }
            }
        switch(lp->word->len){
            case 0: /* ignore */
                break;
            case 1:
                word = lp->word->pdata[0];
                if(!strcasecmp(word, "splice"))
                    sp->model_splice_after = sp->model_length;
                else
                    g_error("bad line in splice data file");
                break;
            case 4:
                for(i = 0; i < 4; i++){
                    num = atof(lp->word->pdata[i]);
                    g_array_append_val(number_list, num);
                    }
                sp->model_length++;
                break;
            default:
                g_error("bad line in splice data file");
                break;
            }
        }
    sp->model_data = NULL;
    LineParse_destroy(lp);
    fclose(fp);
    sp->model_data = (gfloat**)Matrix2d_create(sp->model_length,
                                               5, sizeof(gfloat));
    for(i = 0; i < sp->model_length; i++){
        sp->model_data[i][0] = g_array_index(number_list, gfloat, (i*4));
        sp->model_data[i][1] = g_array_index(number_list, gfloat, (i*4)+1);
        sp->model_data[i][2] = g_array_index(number_list, gfloat, (i*4)+2);
        sp->model_data[i][3] = g_array_index(number_list, gfloat, (i*4)+3);
        }
    g_array_free(number_list, TRUE);
    return;
    }

#if 0
static void SplicePredictor_info(SplicePredictor *sp){
    register gint i;
    g_message("type [%d] len [%d] splice_after [%d]",
            sp->type, sp->model_length, sp->model_splice_after);
    for(i = 0; i < sp->model_length; i++)
        g_message("    [%d] [%2.2f][%2.2f][%2.2f][%2.2f]",
                i,
                sp->model_data[i][0],
                sp->model_data[i][1],
                sp->model_data[i][2],
                sp->model_data[i][3]);
    return;
    }
#endif /* 0 */

static void SplicePredictor_add_5SS_data(SplicePredictor *sp, gchar *path){
    if((!path) || (!strlen(path)) || (!strcasecmp(path, "primate")))
        SplicePredictor_add_data_primate_5SS(sp);
    else
        SplicePredictor_parse_data(sp, path);
    return;
    }

static void SplicePredictor_add_3SS_data(SplicePredictor *sp, gchar *path){
    if((!path) || (!strlen(path)) || (!strcasecmp(path, "primate")))
        SplicePredictor_add_data_primate_3SS(sp);
    else
        SplicePredictor_parse_data(sp, path);
    sp->model_splice_after -= 2;
    return;
    }

/**/

static SplicePredictor_GTAGonly *SplicePredictor_GTAGonly_create(
                                 SplicePredictor *sp){
    register SplicePredictor_GTAGonly *spgo
     = g_new(SplicePredictor_GTAGonly, 1);
    switch(sp->type){
        case SpliceType_ss5_forward:
            spgo->expect_one = 'G';
            spgo->expect_two = 'T';
            break;
        case SpliceType_ss3_forward:
            spgo->expect_one = 'A';
            spgo->expect_two = 'G';
            break;
        case SpliceType_ss5_reverse:
            spgo->expect_one = 'A';
            spgo->expect_two = 'C';
            break;
        case SpliceType_ss3_reverse:
            spgo->expect_one = 'C';
            spgo->expect_two = 'T';
            break;
        default:
            g_error("Unknown splice type [%d]", sp->type);
            break;
        }
    return spgo;
    }

SplicePredictor *SplicePredictor_create(SpliceType type){
    register SplicePredictor *sp = g_new(SplicePredictor, 1);
    register gint i, j, a, z;
    register gfloat swap;
    register Splice_ArgumentSet *sas = Splice_ArgumentSet_create(NULL);
    sp->ref_count = 1;
    sp->type = type;
    if((type == SpliceType_ss5_forward)
    || (type == SpliceType_ss5_reverse)){
        SplicePredictor_add_5SS_data(sp, sas->splice5_data_path);
    } else { /* use SpliceType_3_PRIME */
        SplicePredictor_add_3SS_data(sp, sas->splice3_data_path);
        }
    if((type == SpliceType_ss5_reverse)
    || (type == SpliceType_ss3_reverse)){
        for(a = 0, z = sp->model_length-1; a < z; a++, z--){
            for(j = 0; j < 4; j++){
                swap = sp->model_data[a][j];
                sp->model_data[a][j] = sp->model_data[z][j];
                sp->model_data[z][j] = swap;
                }
            }
        sp->model_splice_after = sp->model_length
                               - sp->model_splice_after
                               - 2;
        }
    for(i = 0; i < (1<<CHAR_BIT); i++)
        sp->index[i] = 4;
    if((type == SpliceType_ss5_forward)
    || (type == SpliceType_ss3_forward)){ /* forward */
        sp->index['A'] = sp->index['a'] = 0;
        sp->index['C'] = sp->index['c'] = 1;
        sp->index['G'] = sp->index['g'] = 2;
        sp->index['T'] = sp->index['t'] = 3;
    } else { /* reverse */
        sp->index['T'] = sp->index['t'] = 0;
        sp->index['G'] = sp->index['g'] = 1;
        sp->index['C'] = sp->index['c'] = 2;
        sp->index['A'] = sp->index['a'] = 3;
        }
    for(i = 0; i < sp->model_length; i++){
         for(j = 0; j < 4; j++){
             sp->model_data[i][j] = (((gfloat)(1+sp->model_data[i][j]))
                                  /((25.0+1.0)));
             sp->model_data[i][j] = log(sp->model_data[i][j]) * 1.5;
             }
         sp->model_data[i][4] = 0.0;
         }
    if(sas->force_gtag)
        sp->gtag_only = SplicePredictor_GTAGonly_create(sp);
    else
        sp->gtag_only = NULL;
    return sp;
    }

void SplicePredictor_destroy(SplicePredictor *sp){
    if(--sp->ref_count)
        return;
    if(sp->gtag_only)
        g_free(sp->gtag_only);
    g_free(sp->model_data);
    g_free(sp);
    return;
    }

SplicePredictor *SplicePredictor_share(SplicePredictor *sp){
    sp->ref_count++;
    return sp;
    }

static gboolean Splice_predict_is_on_GTAG(SplicePredictor *sp,
                          gchar *seq, guint seq_pos){
    register gchar base_one = toupper(seq[seq_pos]),
                   base_two = toupper(seq[seq_pos+1]);
    return ((base_one == sp->gtag_only->expect_one)
         && (base_two == sp->gtag_only->expect_two));
    }

static gfloat Splice_predict_position(SplicePredictor *sp,
                                      gchar *seq, gint seq_len, guint seq_pos){
    register gfloat pos_score, score = 0.0;
    register gint i,
                  seq_start = seq_pos-sp->model_splice_after,
                  model_start = 0,
                  calc_length = sp->model_length;
    if(seq_start < 0){
        model_start = -seq_start;
        seq_start = 0;
        calc_length -= model_start;
        }
    if((seq_start + calc_length) > seq_len){
        calc_length = seq_len-seq_start;
        }
    for(i = 0; i < calc_length; i++){
        pos_score = sp->model_data[model_start+i]
                                  [sp->index[(guchar)seq[seq_start+i]]];
        score += pos_score;
        }
    if(sp->gtag_only &&
       (!Splice_predict_is_on_GTAG(sp, seq, seq_pos)))
        return SPLICE_IMPOSSIBLY_LOW_SCORE;
    return score;
    }

gint SplicePredictor_predict(SplicePredictor *sp,
                             gchar *seq, gint seq_len,
                             guint start, guint length, gfloat *score){
    register guint i;
    register guint max_pos = -1;
    register gfloat curr_score, max_score = SPLICE_IMPOSSIBLY_LOW_SCORE;
    g_assert((start+length) <= seq_len);
    for(i = 0; i < length; i++){
        curr_score = Splice_predict_position(sp, seq, seq_len, start+i);
        if(max_score < curr_score){
            max_score = curr_score;
            max_pos = i;
            }
        }
    if(score)
        *score = max_score;
    return start+max_pos;
    }

void SplicePredictor_predict_array(SplicePredictor *sp,
                                   gchar *seq, guint seq_len,
                                   guint start, guint length,
                                   gfloat *pred){
    register guint i;
    g_assert((start+length) <= seq_len);
    g_assert(seq);
    g_assert(pred);
    for(i = 0; i < length; i++){
        pred[i] = Splice_predict_position(sp, seq, seq_len, start+i);
        }
    return;
    }

static gint SplicePredictor_round(gfloat num){
    return (gint)((num < 0)?(num - 0.5):(num + 0.5));
    }

void SplicePredictor_predict_array_int(SplicePredictor *sp,
                                       gchar *seq, guint seq_len,
                                       guint start, guint length,
                                       gint *pred){
    register guint i;
    register gfloat score;
    g_assert(seq);
    g_assert(length?(pred?TRUE:FALSE):TRUE);
    g_assert((start+length) <= seq_len);
    for(i = 0; i < length; i++){
        score = Splice_predict_position(sp, seq, seq_len, i+start);
        pred[i] = SplicePredictor_round(score);
        }
    return;
    }

gfloat SplicePredictor_get_max_score(SplicePredictor *sp){
    register gfloat score = 0.0, pos_score;
    register gint i, j;
    for(i = 0; i < sp->model_length; i++){
        pos_score = sp->model_data[i][0];
        for(j = 1; j < 4; j++)
            if(pos_score < sp->model_data[i][j])
                pos_score = sp->model_data[i][j];
        score += pos_score;
        }
    return score;
    }

SplicePredictorSet *SplicePredictorSet_create(void){
    register SplicePredictorSet *sps = g_new(SplicePredictorSet, 1);
    sps->ss5_forward = SplicePredictor_create(SpliceType_ss5_forward);
    sps->ss5_reverse = SplicePredictor_create(SpliceType_ss5_reverse);
    sps->ss3_forward = SplicePredictor_create(SpliceType_ss3_forward);
    sps->ss3_reverse = SplicePredictor_create(SpliceType_ss3_reverse);
    return sps;
    }

void SplicePredictorSet_destroy(SplicePredictorSet *sps){
    SplicePredictor_destroy(sps->ss5_forward);
    SplicePredictor_destroy(sps->ss5_reverse);
    SplicePredictor_destroy(sps->ss3_forward);
    SplicePredictor_destroy(sps->ss3_reverse);
    g_free(sps);
    return;
    }

/* TODO: o Maybe add support for branch signal prediction
 *       o Add support for different organisms
 */

/**/

static gpointer SplicePredictor_cache_get_func(gint pos, gpointer page_data,
                                               gpointer user_data){
    register gint *data = page_data;
    return GINT_TO_POINTER(data[pos]);
    }

static SparseCache_Page *SplicePrediction_cache_fill_func(gint start,
                                                          gpointer user_data){
    register SplicePrediction *spn = user_data;
    register SparseCache_Page *page = g_new(SparseCache_Page, 1);
    register gint seq_start, seq_length, pred_start, pred_length;
    register gchar *seq;
    page->get_func = SplicePredictor_cache_get_func;
    page->copy_func = NULL;
    pred_start = MAX(0, start);
    pred_length = MIN(spn->length, start+SparseCache_PAGE_SIZE)
                - pred_start;
    seq_start = MAX(0, pred_start-spn->sp->model_splice_after);
    seq_length = MIN(spn->length, pred_start+pred_length
                     +(spn->sp->model_length-spn->sp->model_splice_after))
               - seq_start;
    page->data = g_new(gint, pred_length);
    page->data_size = sizeof(gint)*pred_length;
    seq = spn->get_seq_func(seq_start, seq_length, spn->user_data);
    SplicePredictor_predict_array_int(spn->sp, seq, seq_length,
                                pred_start-seq_start, pred_length, page->data);
    g_free(seq);
    return page;
    }

SplicePrediction *SplicePrediction_create(SplicePredictor *sp, gint len,
                                  SplicePrediction_get_seq_func get_seq_func,
                                  gboolean use_single, gpointer user_data){
    register SplicePrediction *spn = g_new(SplicePrediction, 1);
    register gchar *seq;
    g_assert(get_seq_func);
    g_assert(len > 0);
    spn->sp = SplicePredictor_share(sp);
    spn->length = len;
    spn->get_seq_func = get_seq_func;
    spn->user_data = user_data;
    if(use_single){
        spn->single = g_new(gint, len);
        seq = get_seq_func(0, len, user_data);
        SplicePredictor_predict_array_int(sp, seq, len, 0, len, spn->single);
        g_free(seq);
        spn->cache = NULL;
    } else {
        spn->single = NULL;
        spn->cache = SparseCache_create(len,
                         SplicePrediction_cache_fill_func, NULL, NULL, spn);
        }
    return spn;
    }

void SplicePrediction_destroy(SplicePrediction *spn){
    SplicePredictor_destroy(spn->sp);
    if(spn->single)
        g_free(spn->single);
    if(spn->cache)
        SparseCache_destroy(spn->cache);
    g_free(spn);
    return;
    }

gint SplicePrediction_get(SplicePrediction *spn, gint pos){
    g_assert(spn);
    if(spn->single){ /* Have precomputed prediction */
        g_assert(pos >= 0);
        g_assert(pos < spn->length);
        return spn->single[pos];
        }
    g_assert(spn->cache);
    return GPOINTER_TO_INT(SparseCache_get(spn->cache, pos));
    }
/* FIXME: optimisation: remove conditional in codegen */

/**/

SplicePrediction_Set *SplicePrediction_Set_add(SplicePrediction_Set *spns,
                         SplicePredictorSet *sps, SpliceType type,
                         gint len, SplicePrediction_get_seq_func get_seq_func,
                         gboolean use_single, gpointer user_data){
    if(!spns)
        spns = g_new0(SplicePrediction_Set, 1);
    switch(type){
        case SpliceType_ss5_forward:
            if(!spns->ss5_forward)
                spns->ss5_forward = SplicePrediction_create(sps->ss5_forward,
                    len, get_seq_func, use_single, user_data);
            break;
        case SpliceType_ss3_forward:
            if(!spns->ss3_forward)
                spns->ss3_forward = SplicePrediction_create(sps->ss3_forward,
                    len, get_seq_func, use_single, user_data);
            break;
        case SpliceType_ss5_reverse:
            if(!spns->ss5_reverse)
                spns->ss5_reverse = SplicePrediction_create(sps->ss5_reverse,
                    len, get_seq_func, use_single, user_data);
            break;
        case SpliceType_ss3_reverse:
            if(!spns->ss3_reverse)
                spns->ss3_reverse = SplicePrediction_create(sps->ss3_reverse,
                    len, get_seq_func, use_single, user_data);
            break;
        default:
            g_error("Unknown SpliceType [%d]", type);
            break;
        }
    return spns;
    }

void SplicePrediction_Set_destroy(SplicePrediction_Set *sps){
    if(sps->ss5_forward)
        SplicePrediction_destroy(sps->ss5_forward);
    if(sps->ss5_reverse)
        SplicePrediction_destroy(sps->ss5_reverse);
    if(sps->ss3_forward)
        SplicePrediction_destroy(sps->ss3_forward);
    if(sps->ss3_reverse)
        SplicePrediction_destroy(sps->ss3_reverse);
    g_free(sps);
    return;
    }

/**/

