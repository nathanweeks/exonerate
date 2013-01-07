/****************************************************************\
*                                                                *
*  C4 dynamic programming library - code for regions             *
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

#ifndef INCLUDED_REGION_H
#define INCLUDED_REGION_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <glib.h>

typedef struct {
     gint ref_count;    /* -1 when static region */
     gint query_start;
     gint target_start;
     gint query_length;
     gint target_length;
} Region;

gboolean  Region_is_valid(Region *region);
  Region *Region_create_blank(void);
  Region *Region_create(gint query_start, gint target_start,
                      gint query_length, gint target_length);
    void  Region_destroy(Region *region);
Region *Region_share(Region *region);
Region *Region_copy(Region *region);

  void  Region_set_static(Region *region);
  void  Region_init_static(Region *region,
                           gint query_start, gint target_start,
                           gint query_length, gint target_length);

gboolean Region_is_within(Region *outer, Region *inner);
gboolean Region_is_same(Region *region_a, Region *region_b);

void Region_print(Region *region, gchar *name);

#define Region_area(region) \
        ((region)->query_length * (region)->target_length)

#define Region_query_end(region) \
        ((region)->query_start+(region)->query_length)

#define Region_target_end(region) \
        ((region)->target_start+(region)->target_length)

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* INCLUDED_REGION_H */

