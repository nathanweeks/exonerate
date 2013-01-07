/****************************************************************\
*                                                                *
*  Deja-vu library for fast linear space repeat finding          *
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

/* Exact repeat finding in large sequences.
 *
 * This is speeded up variant of Roger Sayle's linear space
 * repeat finding algorithm - this algorithm only requires
 * O(n+k) time, where n = seq length, k = number of repeats.
 */

#include <string.h> /* For memset() */
#include "dejavu.h"

DejaVu *DejaVu_create(gchar *seq, gint len){
    register DejaVu *dv = g_new(DejaVu, 1);
    g_assert(seq);
    g_assert(len >= 0);
    dv->seq = seq;
    dv->len = len;
    dv->next = g_new0(gint, len);
    dv->symbol_list_len = 0;
    return dv;
    }

void DejaVu_destroy(DejaVu *dv){
    g_free(dv->next);
    g_free(dv);
    return;
    }

void DejaVu_info(DejaVu *dv){
    register gint i;
    g_print(" Seq:");
    for(i = 0; i < dv->len; i++)
        g_print("  %c", dv->seq[i]);
    g_print("\n");
    g_print(" Pos:");
    for(i = 0; i < dv->len; i++)
        g_print("%3d", i);
    g_print("\n");
    g_print("Next:");
    for(i = 0; i < dv->len; i++)
        g_print("%3d", dv->next[i]);
    g_print("\n");
    return;
    }

static gint DejaVu_init(DejaVu *dv, guchar *filter){
    register gint i, ch, start = -1, end = -1;
    gint first[ALPHABET_SIZE];
    gint last[ALPHABET_SIZE];
    for(i = 0; i < ALPHABET_SIZE; i++)
        first[i] = last[i] = -1;
    for(i = 0; i < dv->len; i++){ /* Set next char */
        ch = dv->seq[i];
        if(filter[ch] == '-')
            continue;
        if(first[ch] == -1)
            first[ch] = i;
        if(last[ch] != -1)
            dv->next[last[ch]] = i;
        last[ch] = i;
        }
    for(i = 0; i < ALPHABET_SIZE; i++){
        if(first[i] != -1){ /* exists */
            dv->symbol_list[dv->symbol_list_len++] = i;
            if(first[i] != last[i]){ /* is repeated */
                if(end != -1)
                    dv->next[end] = -first[i];
                if(start == -1)
                    start = first[i];
                end = last[i];
                }
            }
        }
    if(end == -1) /* In case there are no repeats */
        return -1;
    dv->next[end] = -start;
    return start;
    }
/* dv->next holds the position of next repeat instance
 * except:
 *     when (next < 0)
 *         next holds -position of next repeat species
 *     when (next == -start)
 *         this is the last repeat of this length
 */

static gint DejaVu_promote_species(DejaVu *dv, gint curr,
                     gint length, gint min_wordlen,
                     DejaVu_TraverseFunc dvtf, guchar *filter,
                     gpointer user_data,
                     gint *first, gint *last){
    register gint ch, pos = curr;
    do {
        if(length >= min_wordlen)
            dvtf(curr, pos, length, dv->seq, dv->len, user_data); /* report */
        if((pos+length) < dv->len){
            ch = dv->seq[pos+length];
            if(filter[ch] != '-'){
                if(first[ch] == -1) /* record first */
                    first[ch] = pos;
                if(last[ch] != -1) /* link from prev */
                    dv->next[last[ch]] = pos;
                last[ch] = pos; /* record this */
                }
            }
        pos = dv->next[pos];
    } while(pos > 0);
    return -pos;
    }
/* FIXME: optimisation: could put dvrf() in separate loop */

static gint DejaVu_promote(DejaVu *dv, gint start,
                       gint length, gint min_wordlen,
                       DejaVu_TraverseFunc dvtf, guchar *filter,
                       gpointer user_data){
    register gint i, ch, pos = start, end = -1, new_start = -1;
    gint first[ALPHABET_SIZE];
    gint last[ALPHABET_SIZE];
    for(i = 0; i < dv->symbol_list_len; i++){
        ch = dv->symbol_list[i];
        first[ch] = last[ch] = -1;
        }
    do {
        pos = DejaVu_promote_species(dv, pos, length, min_wordlen,
                                 dvtf, filter, user_data, first, last);
        for(i = 0; i < dv->symbol_list_len; i++){
            ch = dv->symbol_list[i];
            if(first[ch] != -1){ /* exists */
                if(first[ch] != last[ch]){ /* is repeated */
                    if(end != -1)
                        dv->next[end] = -first[ch];
                    if(new_start == -1)
                        new_start = first[ch];
                    end = last[ch];
                    }
                first[ch] = -1;
                }
            last[ch] = -1;
            }
    } while(pos != start);
    if(end != -1)
        dv->next[end] = -new_start;
    return new_start;
    }

void DejaVu_traverse(DejaVu *dv,
                     gint min_wordlen, gint max_wordlen,
                     DejaVu_TraverseFunc dvtf, gpointer user_data,
                     guchar *filter, gint verbosity){
    register gint i, start, length = 1;
    register guchar *plain_filter = NULL;
    g_assert(min_wordlen >= 0);
    g_assert(min_wordlen <= max_wordlen);
    g_assert(dv);
    g_assert(dvtf);
    if(!filter){
        plain_filter = g_new(guchar, ALPHABET_SIZE);
        for(i = 0; i < ALPHABET_SIZE; i++)
            plain_filter[i] = i;
        filter = plain_filter;
        }
    start = DejaVu_init(dv, filter);
    if(verbosity > 0)
        g_print("Message: Processing [");
    while((start != -1) /* while there are more repeats */
        && (length <= max_wordlen)){
        if(verbosity > 0)
            g_print(".");
        start = DejaVu_promote(dv, start, length++, min_wordlen,
                               dvtf, filter, user_data);
        }
    if(verbosity > 0)
        g_print("]\n");
    if(plain_filter)
        g_free(plain_filter);
    return;
    }

/**/

