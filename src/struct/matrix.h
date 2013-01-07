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

#ifndef INCLUDED_MATRIX_H
#define INCLUDED_MATRIX_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <glib.h>

gsize Matrix2d_size(gint a, gint b,                 gsize cell);
gsize Matrix3d_size(gint a, gint b, gint c,         gsize cell);
gsize Matrix4d_size(gint a, gint b, gint c, gint d, gsize cell);

/* Sizes may be unexpected due to word alignment padding */

gpointer   *Matrix2d_create(gint a, gint b, gsize cell);
    void    Matrix2d_init(gchar **m, gint a, gint b, gsize cell);
gpointer  **Matrix3d_create(gint a, gint b, gint c, gsize cell);
gpointer ***Matrix4d_create(gint a, gint b, gint c, gint d, gsize cell);

/* All the matrices are freed with a single call to g_free() */

/* FIXME: implement general N-dimentional versions */

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* INCLUDED_MATRIX_H */

