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

#ifndef INCLUDED_SPLICE_H
#define INCLUDED_SPLICE_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <glib.h>
#include <stdio.h>
#include <limits.h>

#include  "sparsecache.h"
#include  "argument.h"

#ifndef ALPHABETSIZE
#define ALPHABETSIZE (1<<(CHAR_BIT))
#endif /* ALPHABETSIZE */

typedef struct {
       gchar *splice3_data_path;
       gchar *splice5_data_path;
    gboolean  force_gtag;
} Splice_ArgumentSet;
/* When gtag_only is TRUE, only standard GT...AG splice sites
 * will be assigned a good prediction score.
 */

Splice_ArgumentSet *Splice_ArgumentSet_create(Argument *arg);

/**/

typedef enum {
    SpliceType_ss5_forward,
    SpliceType_ss3_forward,
    SpliceType_ss5_reverse,
    SpliceType_ss3_reverse
} SpliceType;

typedef struct {
    gchar expect_one;
    gchar expect_two;
} SplicePredictor_GTAGonly;

typedef struct {
                        gint   ref_count;
                  SpliceType   type;
                        gint   model_length;
                        gint   model_splice_after;
                      gfloat **model_data;
                      guchar   index[ALPHABETSIZE];
    SplicePredictor_GTAGonly  *gtag_only;
} SplicePredictor;

SplicePredictor *SplicePredictor_create(SpliceType type);
           void  SplicePredictor_destroy(SplicePredictor *sp);
SplicePredictor *SplicePredictor_share(SplicePredictor *sp);

gint SplicePredictor_predict(SplicePredictor *sp,
                             gchar *seq, gint seq_len,
                             guint start, guint length, gfloat *score);
/* Returns the position, and fills the score when provided */

void SplicePredictor_predict_array(SplicePredictor *sp,
                                   gchar *seq, guint seq_len,
                                   guint start, guint length,
                                   gfloat *pred);
/* Prediction score will be written to <pred>
 * which must contain space for <length> predictions.
 */

void SplicePredictor_predict_array_int(SplicePredictor *sp,
                                       gchar *seq, guint seq_len,
                                       guint start, guint length,
                                       gint *pred);
/* As SplicePredictor_predict_array(),
 * but results are cast as (int)(score*factor)
 */

gfloat SplicePredictor_get_max_score(SplicePredictor *sp);
/* Returns maximum predicton score for any sequence
 */

typedef struct {
    SplicePredictor *ss5_forward;
    SplicePredictor *ss5_reverse;
    SplicePredictor *ss3_forward;
    SplicePredictor *ss3_reverse;
} SplicePredictorSet;

SplicePredictorSet *SplicePredictorSet_create(void);
              void  SplicePredictorSet_destroy(SplicePredictorSet *sps);

/**/
typedef gchar *(*SplicePrediction_get_seq_func)
               (gint start, gint len, gpointer user_data);

typedef struct {
                 SplicePredictor *sp;
                            gint *single;
                            gint  length;
                     SparseCache *cache;
   SplicePrediction_get_seq_func  get_seq_func;
                        gpointer  user_data;
} SplicePrediction;

SplicePrediction *SplicePrediction_create(SplicePredictor *sp, gint len,
                                   SplicePrediction_get_seq_func get_seq_func,
                                   gboolean use_single, gpointer user_data);
void SplicePrediction_destroy(SplicePrediction *spn);
gint SplicePrediction_get(SplicePrediction *spn, gint pos);

typedef struct {
    SplicePrediction *ss5_forward;
    SplicePrediction *ss5_reverse;
    SplicePrediction *ss3_forward;
    SplicePrediction *ss3_reverse;
} SplicePrediction_Set;

SplicePrediction_Set *SplicePrediction_Set_add(SplicePrediction_Set *spns,
                                               SplicePredictorSet *sps,
                         SpliceType type,
                         gint len, SplicePrediction_get_seq_func get_seq_func,
                         gboolean use_single, gpointer user_data);
                void  SplicePrediction_Set_destroy(SplicePrediction_Set *sps);

/**/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* INCLUDED_SPLICE_H */

