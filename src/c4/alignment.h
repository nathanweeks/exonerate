/****************************************************************\
*                                                                *
*  C4 dynamic programming library - alignment code               *
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

#ifndef INCLUDED_ALIGNMENT_H
#define INCLUDED_ALIGNMENT_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <glib.h>

#include "c4.h"
#include "region.h"
#include "sequence.h"
#include "submat.h"
#include "translate.h"
#include "argument.h"

typedef struct {
        gint alignment_width;
    gboolean forward_strand_coords;
} Alignment_ArgumentSet;

Alignment_ArgumentSet *Alignment_ArgumentSet_create(Argument *arg);

typedef struct {
    C4_Transition *transition;
             gint  length;
} AlignmentOperation;

typedef struct {
       Region *region;
         gint  ref_count;
    GPtrArray *operation_list;
     C4_Score  score;
     C4_Model *model;
} Alignment;

Alignment *Alignment_create(C4_Model *model, Region *region,
                            C4_Score score);
Alignment *Alignment_share(Alignment *alignment);
     void  Alignment_destroy(Alignment *alignment);

void Alignment_add(Alignment *alignment, C4_Transition *transition,
                   gint length);

void Alignment_display(Alignment *alignment,
                       Sequence *query, Sequence *target,
                       Submat *dna_submat, Submat *protein_submat,
                       Translate *translate, FILE *fp);
/* submat and translate may be NULL
 * when not required by the alignment model
 */

void Alignment_display_sugar(Alignment *alignment,
                             Sequence *query, Sequence *target, FILE *fp);

void Alignment_display_cigar(Alignment *alignment,
                             Sequence *query, Sequence *target, FILE *fp);

void Alignment_display_vulgar(Alignment *alignment,
                              Sequence *query, Sequence *target, FILE *fp);

void Alignment_display_gff(Alignment *alignment,
                           Sequence *query, Sequence *target,
                           Translate *translate,
                           gboolean report_on_query,
                           gboolean report_on_genomic,
                           gint result_id, gpointer user_data, FILE *fp);

void Alignment_display_ryo(Alignment *alignment,
        Sequence *query, Sequence *target, gchar *format,
        Translate *translate, gint rank,
        gpointer user_data, gpointer self_data, FILE *fp);

/**/

gboolean Alignment_is_valid(Alignment *alignment, Region *region,
                            gpointer user_data);

void Alignment_import_derived(Alignment *alignment,
                              Alignment *to_add,
                              C4_DerivedModel *derived_model);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* INCLUDED_ALIGNMENT_H */

