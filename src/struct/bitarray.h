/****************************************************************\
*                                                                *
*  A simple bitarray data structure                              *
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

#ifndef INCLUDED_BITARRAY_H
#define INCLUDED_BITARRAY_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <glib.h>
#include <stdio.h>
#include <limits.h> /* For CHAR_BIT */

typedef struct {
    guchar *data;
   guint64  length; /* Length in bits */
   guint64  alloc;
} BitArray;

#define BitArray_get_size(length)                    \
   (  ((length) / (sizeof(guchar) * CHAR_BIT))       \
    +(((length) % (sizeof(guchar) * CHAR_BIT))?1:0))

#define BitArray_memory_usage(ba)         \
            sizeof(BitArray)              \
            + (sizeof(gchar)*(ba)->alloc)

BitArray *BitArray_create(void);
    void  BitArray_destroy(BitArray *ba);
    void  BitArray_info(BitArray *ba);
    void  BitArray_empty(BitArray *ba);
    void  BitArray_write(BitArray *ba, FILE *fp);
BitArray *BitArray_read(FILE *fp, gsize size);

    void  BitArray_append_bit(BitArray *ba, gboolean bit);
    void  BitArray_append(BitArray *ba, guint64 data, guchar width);
gboolean  BitArray_get_bit(BitArray *ba, guint64 pos);
 guint64  BitArray_get(BitArray *ba, guint64 start, guchar width);

    void  BitArray_write_int(guint64 num, FILE *fp);
 guint64  BitArray_read_int(FILE *fp);
    void  BitArray_int_print(guint64 num);
   gchar *BitArray_int_sprintf(guint64 num);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* INCLUDED_BITARRAY_H */

