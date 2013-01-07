/****************************************************************\
*                                                                *
*  Simple matrix creation routines                               *
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

#include "matrix.h"

static void Matrix_layout(gpointer *index, gint length,
                          gchar *data, gsize block){
    register gint i;
    register gchar *p;
    g_assert(index);
    g_assert(data);
    g_assert(length > 0);
    g_assert(block > 0);
    for(i = 0, p = data; i < length; i++){
        index[i] = p;
        p += block;
        }
    return;
    }

/**/

gsize Matrix2d_size(gint a, gint b, gsize cell){
    register gulong index_size = a * sizeof(gpointer),
                    row_size = b * cell,
                    data_size = a * row_size,
                    total_size = index_size + data_size;
    register gdouble check_index_size = a * sizeof(gpointer),
                     check_row_size = b * cell,
                     check_data_size = a * check_row_size,
                     check_total_size = check_index_size
                                      + check_data_size;
    g_assert(a > 0);
    g_assert(b > 0);
    g_assert(cell > 0);
    if((check_total_size - total_size) > 1) /* overflow */
        return 0;
    return total_size;
    }

void Matrix2d_init(gchar **m, gint a, gint b, gsize cell){
    Matrix_layout((gpointer*)m, a, (gchar*)&m[a], b*cell);
    return;
    }

gpointer *Matrix2d_create(gint a, gint b, gsize cell){
    register gchar **m;
    register gulong index_size = a * sizeof(gpointer),
                    row_size = b * cell,
                    data_size = a * row_size;
    g_assert(a > 0);
    g_assert(b > 0);
    g_assert(cell > 0);
    g_assert(Matrix2d_size(a, b, cell));
    m = g_malloc0(index_size + data_size);
    Matrix_layout((gpointer*)m, a, (gchar*)&m[a], row_size);
    return (gpointer*)m;
    }

/**/

gsize Matrix3d_size(gint a, gint b, gint c, gsize cell){
    register gulong primary_index_size = a * sizeof(gpointer),
                    secondary_index_size = b * sizeof(gpointer),
                    row_size  = c * cell,
                    block_size = secondary_index_size
                               + (b * row_size),
                    total_size;
    register gdouble check_primary_index_size = a * sizeof(gpointer),
                     check_secondary_index_size = b * sizeof(gpointer),
                     check_row_size  = c * cell,
                     check_block_size = check_secondary_index_size
                                      + (b * check_row_size),
                     check_total_size;
    g_assert(a > 0);
    g_assert(b > 0);
    g_assert(c > 0);
    g_assert(cell > 0);
    block_size += (block_size % sizeof(gpointer)); /* Pad */
    check_block_size += (block_size % sizeof(gpointer)); /* Pad */
    total_size = primary_index_size + (a * block_size);
    check_total_size = check_primary_index_size
                     + (a * check_block_size);
    if((check_total_size - total_size) > 1) /* overflow */
        return 0;
    return total_size;
    }

gpointer **Matrix3d_create(gint a, gint b, gint c, gsize cell){
    register gint i;
    register gchar ***m;
    register gulong primary_index_size = a * sizeof(gpointer),
                    secondary_index_size = b * sizeof(gpointer),
                    row_size  = c * cell,
                    block_size = secondary_index_size
                              + (b * row_size),
                    total_size;

    g_assert(a > 0);
    g_assert(b > 0);
    g_assert(c > 0);
    g_assert(cell > 0);
    g_assert(Matrix3d_size(a, b, c, cell));
    block_size += (block_size % sizeof(gpointer)); /* Pad */
    total_size = primary_index_size + (a * block_size);
    m = g_malloc0(total_size);
    Matrix_layout((gpointer*)m, a, (gchar*)&m[a], block_size);
    for(i = 0; i < a; i++)
        Matrix_layout((gpointer*)m[i], b,
                      (gchar*)m[i]+secondary_index_size, row_size);
    return (gpointer**)m;
    }

/**/

gsize Matrix4d_size(gint a, gint b, gint c, gint d, gsize cell){
    register gulong primary_index_size = a * sizeof(gpointer),
                    secondary_index_size = b * sizeof(gpointer),
                    tertiary_index_size = c * sizeof(gpointer),
                    row_size  = d * cell,
                    block_size = tertiary_index_size
                               + (c * row_size),
                    sheet_size, total_size;
    register gdouble check_primary_index_size = a * sizeof(gpointer),
                     check_secondary_index_size = b * sizeof(gpointer),
                     check_tertiary_index_size = c * sizeof(gpointer),
                     check_row_size  = d * cell,
                     check_block_size = check_tertiary_index_size
                                      + (c * check_row_size),
                     check_sheet_size, check_total_size;
    g_assert(a > 0);
    g_assert(b > 0);
    g_assert(c > 0);
    g_assert(d > 0);
    g_assert(cell > 0);
    block_size += (block_size % sizeof(gpointer)); /* Pad */
    check_block_size += (block_size % sizeof(gpointer)); /* Pad */
    sheet_size = secondary_index_size + (b * block_size),
    check_sheet_size = check_secondary_index_size
                     + (b * check_block_size),
    sheet_size += (sheet_size % sizeof(gpointer)); /* Pad */
    check_sheet_size += (sheet_size % sizeof(gpointer)); /* Pad */
    total_size = primary_index_size + (a * sheet_size);
    check_total_size = check_primary_index_size
                     + (a * check_sheet_size);
    if((check_total_size - total_size) > 1) /* overflow */
        return 0;
    return total_size;
    }

gpointer ***Matrix4d_create(gint a, gint b, gint c, gint d, gsize cell){
    register gint i, j;
    register gchar ***m;
    register gpointer *ptr;
    register gulong primary_index_size = a * sizeof(gpointer),
                    secondary_index_size = b * sizeof(gpointer),
                    tertiary_index_size = c * sizeof(gpointer),
                    row_size  = d * cell,
                    block_size = tertiary_index_size
                               + (c * row_size),
                    sheet_size, total_size;
    g_assert(a > 0);
    g_assert(b > 0);
    g_assert(c > 0);
    g_assert(d > 0);
    g_assert(cell > 0);
    g_assert(Matrix4d_size(a, b, c, d, cell));
    block_size += (block_size % sizeof(gpointer)); /* Pad */
    sheet_size = secondary_index_size + (b * block_size),
    sheet_size += (sheet_size % sizeof(gpointer)); /* Pad */
    total_size = primary_index_size + (a * sheet_size);
    m = g_malloc0(total_size);
    Matrix_layout((gpointer*)m, a, (gchar*)&m[a], sheet_size);
    for(i = 0; i < a; i++){
        ptr = (gpointer*)m[i];
        Matrix_layout(ptr, b,
                      (gchar*)ptr+secondary_index_size, block_size);
        for(j = 0; j < b; j++){
            Matrix_layout((gpointer*)ptr[j], c,
                          (gchar*)ptr[j]+tertiary_index_size, row_size);
            }
        }
    return (gpointer***)m;
    }

/* FIXME: tidy and reduce redundancy between create and size functions
 */

