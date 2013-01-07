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


#include "region.h"

gboolean Region_is_valid(Region *region){
    g_assert(region);
    g_assert(region->query_start >= 0);
    g_assert(region->target_start >= 0);
    g_assert(region->query_length >= 0);
    g_assert(region->target_length >= 0);
    return TRUE;
    }

Region *Region_create_blank(void){
    register Region *region = g_new0(Region, 1);
    region->ref_count = 1;
    return region;
    }

Region *Region_create(gint query_start, gint target_start,
                      gint query_length, gint target_length){
    register Region *region = g_new(Region, 1);
    region->ref_count = 1;
    region->query_start = query_start;
    region->target_start = target_start;
    region->query_length = query_length;
    region->target_length = target_length;
    g_assert(Region_is_valid(region));
    return region;
    }

void Region_destroy(Region *region){
    g_assert(region);
    g_assert(region->ref_count != -1); /* Check not static */
    if(--region->ref_count)
        return;
    g_free(region);
    return;
    }

Region *Region_share(Region *region){
    g_assert(region);
    g_assert(region->ref_count != -1); /* Check not static */
    region->ref_count++;
    return region;
    }

Region *Region_copy(Region *region){
    g_assert(region);
    return Region_create(region->query_start, region->target_start,
                         region->query_length, region->target_length);
    }

void Region_set_static(Region *region){
    region->ref_count = -1;
    return;
    }

void Region_init_static(Region *region,
                        gint query_start, gint target_start,
                        gint query_length, gint target_length){
    region->ref_count = -1;
    region->query_start = query_start;
    region->target_start = target_start;
    region->query_length = query_length;
    region->target_length = target_length;
    g_assert(Region_is_valid(region));
    return;
    }

gboolean Region_is_within(Region *outer, Region *inner){
    g_assert(outer);
    g_assert(inner);
    if(outer->query_start > inner->query_start)
        return FALSE;
    if(outer->target_start > inner->target_start)
        return FALSE;
    if(Region_query_end(outer) < Region_query_end(inner))
        return FALSE;
    if(Region_target_end(outer) < Region_target_end(inner))
        return FALSE;
    return TRUE;
    }

gboolean Region_is_same(Region *region_a, Region *region_b){
    g_assert(region_a);
    g_assert(region_b);
    if(region_a->query_start != region_b->query_start)
        return FALSE;
    if(region_a->target_start != region_b->target_start)
        return FALSE;
    if(Region_query_end(region_a) != Region_query_end(region_b))
        return FALSE;
    if(Region_target_end(region_a) != Region_target_end(region_b))
        return FALSE;
    return TRUE;
    }

void Region_print(Region *region, gchar *name){
    g_print("Region [%s] q[%d(%d)%d] t[%d(%d)%d]\n",
            name,
            region->query_start,
            region->query_length,
            Region_query_end(region),
            region->target_start,
            region->target_length,
            Region_target_end(region));
    /* FIXME: temp */
    g_print("draw_region(%d, %d, %d, %d, \"%s\")\n",
            region->query_start, region->target_start,
            region->query_length, region->target_length,
            name);
    return;
    }

