/****************************************************************\
*                                                                *
*  C4 dynamic programming library - sub-alignment regions        *
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

#ifndef INCLUDED_SAR_H
#define INCLUDED_SAR_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <glib.h>

#include "c4.h"
#include "hpair.h"
#include "hspset.h"
#include "region.h"

typedef struct {
    gfloat hsp_quality;
} SAR_ArgumentSet;

SAR_ArgumentSet *SAR_ArgumentSet_create(Argument *arg);

/*  The five types of SAR used in C4:
 *
 *  +-----+
 *  |     | START
 *  |    \|
 *  +-----\
 *         x
 *          \-----+
 *          |\  \ | JOIN
 *          | \  \|
 *          +-----\
 *                 x
 *                  \-----+ ....... +-----+
 *                  |\    |         |     |
 *             SRC  | \   |  (span) |   \ |  DST
 *                  |     |         |    \|
 *                  +-----+ ....... +-----\
 *                                         x
 *                                          \-----+
 *                                          |\    |
 *                                          | \   | END
 *                                          |     |
 *                                          +-----+
 */

typedef struct {
        Region *region;
      C4_Score  component;
} SAR_Terminal;

SAR_Terminal *SAR_Terminal_create(HSP *hsp, HPair *hpair,
                Heuristic_Match *match, gboolean is_start);
void SAR_Terminal_destroy(SAR_Terminal *sar_terminal);

C4_Score SAR_Terminal_find_bound(SAR_Terminal *sar_terminal,
                                 Heuristic_Bound *heuristic_bound);

C4_Score SAR_Terminal_find_score(SAR_Terminal *sar_terminal,
                                 Optimal *optimal, HPair *hpair);

/**/

typedef struct {
            Region *region;
          C4_Score  src_component;
          C4_Score  dst_component;
    Heuristic_Pair *pair;
} SAR_Join;

SAR_Join *SAR_Join_create(HSP *src_hsp, HSP *dst_hsp, HPair *hpair,
                          Heuristic_Pair *pair);
    void  SAR_Join_destroy(SAR_Join *sar_join);
C4_Score  SAR_Join_find_bound(SAR_Join *sar_join);
C4_Score  SAR_Join_find_score(SAR_Join *sar_join, HPair *hpair);

/**/

typedef struct {
            Region *src_region;
            Region *dst_region;
          C4_Score  src_component;
          C4_Score  dst_component;
    Heuristic_Span *span;
} SAR_Span;

SAR_Span *SAR_Span_create(HSP *src_hsp, HSP *dst_hsp, HPair *hpair,
                          Heuristic_Span *span,
                          C4_Portal *src_portal, C4_Portal *dst_portal);
void SAR_Span_destroy(SAR_Span *sar_span);
C4_Score  SAR_Span_find_bound(SAR_Span *sar_span);
C4_Score  SAR_Span_find_score(SAR_Span *sar_span, HPair *hpair);

/**/

typedef struct {
          Alignment *alignment;
          Alignment *end_alignment;
             Region *end_region;
    Heuristic_Match *end_match;
              HPair *hpair;
                HSP *last_hsp;
             Region *last_region;
    Heuristic_Match *last_match;
} SAR_Alignment;

SAR_Alignment *SAR_Alignment_create(SAR_Terminal *sar_start,
                                    SAR_Terminal *sar_end,
                                    Heuristic_Match *start_match,
                                    Heuristic_Match *end_match,
                                    HPair *hpair, C4_Score score);
void SAR_Alignment_finalise(SAR_Alignment *sar_alignment);
void SAR_Alignment_destroy(SAR_Alignment *sar_alignment);

void SAR_Alignment_add_SAR_Join(SAR_Alignment *sar_alignment,
                                SAR_Join *sar_join);
void SAR_Alignment_add_SAR_Span(SAR_Alignment *sar_alignment,
                                SAR_Span *sar_span);
void SAR_Alignment_add_HSP(SAR_Alignment *sar_alignment, HSP *hsp,
                           Heuristic_Match *match);

/**/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* INCLUDED_SAR_H */

