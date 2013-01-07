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

#include <string.h> /* For memset */
#include "bitarray.h"

/**/

#ifndef GetBit
#define GetBit(i,n) (((i)>>(n))&1)
#endif /* GetBit */

#ifndef BitOn
#define BitOn(i,n) ((i)|1<<(n))
#endif /* BitOn */

#ifndef BitOff
#define BitOff(i,n) ((i)&~(1<<(n)))
#endif /* BitOff */

#ifndef GetBits
#define GetBits(x,y,z) (((x)>>(y))&~(~0<<(z)))
#endif /* GetBits */

#ifndef SetBits
#define SetBits(w,x,y,z) (((x)&~(~(~0<<(z))<<(y)))|((w)<<(y)))
#endif /* SetBits */

/**/

BitArray *BitArray_create(void){
    register BitArray *ba = g_new(BitArray, 1);
    ba->alloc = 16;
    ba->length = 0;
    ba->data = g_new0(guchar, ba->alloc);
    return ba;
    }

void BitArray_destroy(BitArray *ba){
    g_free(ba->data);
    g_free(ba);
    return;
    }

void BitArray_info(BitArray *ba){
    register gint64 i;
    g_print("BitArray length [%d] alloc [%d] data [",
            (gint)ba->length, (gint)ba->alloc);
    for(i = 0; i < ba->length; i++)
        g_print("%d", BitArray_get_bit(ba, i));
    g_print("]\n");
    return;
    }

void BitArray_empty(BitArray *ba){
    ba->length = 0;
    return;
    }

void BitArray_write(BitArray *ba, FILE *fp){
    fwrite(ba->data, sizeof(guchar), BitArray_get_size(ba->length), fp);
    fflush(fp); /* Keep valgrind quiet ?? */
    return;
    }

BitArray *BitArray_read(FILE *fp, gsize size){
    register BitArray *ba = g_new(BitArray, 1);
    ba->alloc = size;
    ba->length = size * CHAR_BIT;
    ba->data = g_new(guchar, ba->alloc);
    fread(ba->data, sizeof(guchar), size, fp);
    return ba;
    }

void BitArray_append_bit(BitArray *ba, gboolean bit){
    register gint64 word = ba->length / CHAR_BIT,
                    pos = ba->length & (CHAR_BIT-1);
    if(ba->length == (ba->alloc*CHAR_BIT)){
        ba->alloc <<= 1;
        ba->data = g_realloc(ba->data, ba->alloc);
        memset(ba->data+(ba->alloc>>1), 0, (ba->alloc>>1));
        }
    ba->data[word] = bit?BitOn(ba->data[word], pos)
                        :BitOff(ba->data[word], pos);
    ba->length++;
    return;
    }

void BitArray_append(BitArray *ba, guint64 data, guchar width){
    register gint64 word = ba->length/ CHAR_BIT, input;
    register guchar todo = width,
                    bit = ba->length & (CHAR_BIT-1),
                    taken = CHAR_BIT-bit,
                    done = 0;
    ba->length += width;
    if(ba->length >= (ba->alloc*CHAR_BIT)){
        ba->alloc <<= 1;
        ba->data = g_realloc(ba->data, ba->alloc);
        memset(ba->data+(ba->alloc>>1), 0, (ba->alloc>>1));
        }
    do {
        if(taken >= todo){ /* final partial */
            input = GetBits(data, done, todo);
            ba->data[word] = SetBits(input, ba->data[word], bit, todo);
            break;
            }
        input = GetBits(data, done, taken);
        ba->data[word] = SetBits(input, ba->data[word], bit, taken);
        todo -= taken;
        done += taken;
        bit = 0;
        taken = CHAR_BIT;
        word++;
    } while(TRUE);
    return;
    }

/* Previous slow version */
#if 0
void BitArray_append(BitArray *ba, guint64 data, guchar width){
    register guchar i;
    g_assert(width <= (sizeof(guint64)*CHAR_BIT));
    for(i = 0; i < width; i++)
        BitArray_append_bit(ba, GetBit(data, i));
    return;
    }
#endif /* 0 */

gboolean BitArray_get_bit(BitArray *ba, guint64 pos){
    register gint64 word = pos / CHAR_BIT,
                    bit = pos & (CHAR_BIT-1);
    return GetBit(ba->data[word], bit);
    }

guint64 BitArray_get(BitArray *ba, guint64 start, guchar width){
    register gint64 word = start / CHAR_BIT,
                    data = 0;
    register guchar todo = width,
                    bit = start & (CHAR_BIT-1),
                    taken = CHAR_BIT-bit,
                    done = 0;
    do {
        if(taken >= todo){ /* final partial */
            data |= ((guint64)GetBits(ba->data[word], bit, todo) << done);
            break;
            }
        data |= ((guint64)GetBits(ba->data[word], bit, taken) << done);
        todo -= taken;
        done += taken;
        bit = 0;
        taken = CHAR_BIT;
        word++;
    } while(TRUE);
    return data;
    }

/* Previous slow version */
#if 0
guint64 BitArray_get(BitArray *ba, guint64 start, guchar width){
    register gint i;
    register guint64 data = 0;
    for(i = 0; i < width; i++){
        data |= ((guint64)BitArray_get_bit(ba, start+i) << i);
        }
    return data;
    }
#endif /* 0 */

void BitArray_write_int(guint64 num, FILE *fp){
    guint64 nbo_num = GUINT64_TO_BE(num); /* BigEndian == NBO */
    fwrite(&nbo_num, sizeof(guint64), 1, fp);
    return;
    }

guint64 BitArray_read_int(FILE *fp){
    guint64 num;
    fread(&num, sizeof(guint64), 1, fp);
    return GUINT64_FROM_BE(num);
    }

/**/

