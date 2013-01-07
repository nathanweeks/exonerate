/****************************************************************\
*                                                                *
*  C4 dynamic programming library - SDP Lookahead Object         *
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

#include "lookahead.h"

/**/

#define Lookahead_index_pos(lookahead, advance) \
    (((lookahead)->pos + (advance)) & (Lookahead_Mask_WIDTH-1))

#define Lookahead_pos_is_set(lookahead, advance)         \
    ((lookahead)->mask                                   \
   & (1 << Lookahead_index_pos((lookahead), (advance))))

Lookahead *Lookahead_create(gint reset_pos, gint max_advance,
                            Lookahead_FreeFunc free_func,
                            gpointer user_data){
    register Lookahead *lookahead = g_new0(Lookahead, 1);
    g_assert(max_advance > 0);
    g_assert(max_advance <= Lookahead_Mask_WIDTH);
    lookahead->reset_pos = reset_pos;
    lookahead->pos = reset_pos;
    lookahead->max_advance = max_advance;
    lookahead->free_func = free_func;
    lookahead->user_data = user_data;
    return lookahead;
    }

static Lookahead_Mask Lookahead_orientate_mask(Lookahead *lookahead){
    register gint index_pos = Lookahead_index_pos(lookahead, 0);
    return (lookahead->mask >> index_pos)
         | (lookahead->mask << (Lookahead_Mask_WIDTH-index_pos));
    }

static gint Lookahead_next_pos(Lookahead *lookahead){
    g_assert(lookahead->mask); /* Not currently called when empty */
    return Lookahead_Mask_ffs(Lookahead_orientate_mask(lookahead)) - 1;
    }
/* Returns advance of next occupied or -1 when nothing is set
 */

void Lookahead_destroy(Lookahead *lookahead){
    if(lookahead->mask){ /* Not empty */
        /* Move to first occupied pos */
        Lookahead_move(lookahead, lookahead->pos
                                + Lookahead_next_pos(lookahead));
        while(lookahead->mask) /* Not empty */
            Lookahead_next(lookahead); /* Move to next position */
        }
    g_free(lookahead);
    return;
    }

void Lookahead_reset(Lookahead *lookahead){
     Lookahead_move(lookahead, lookahead->pos + Lookahead_Mask_WIDTH);
     g_assert(!lookahead->mask);
     lookahead->pos = lookahead->reset_pos;
     return;
     }

gpointer Lookahead_get(Lookahead *lookahead, gint advance){
    g_assert(advance >= 0);
    g_assert(advance <= lookahead->max_advance);
    return lookahead->index[Lookahead_index_pos(lookahead, advance)];
    }

void Lookahead_set(Lookahead *lookahead, gint advance, gpointer data){
    register gint index_pos = Lookahead_index_pos(lookahead, advance);
    g_assert(advance >= 0);
    g_assert(advance <= lookahead->max_advance);
    g_assert(!Lookahead_get(lookahead, advance));
    lookahead->index[index_pos] = data;
    lookahead->mask |= (1 << index_pos);
    return;
    }

static void Lookahead_clear(Lookahead *lookahead, gint advance){
    register gint index_pos = Lookahead_index_pos(lookahead, advance);
    g_assert(advance >= 0);
    g_assert(advance <= lookahead->max_advance);
    g_assert(lookahead->index[index_pos]);
    lookahead->free_func(lookahead->index[index_pos],
                         lookahead->user_data);
    lookahead->index[index_pos] = NULL;
    lookahead->mask &= (~(1 << index_pos));
    return;
    }

gint Lookahead_next(Lookahead *lookahead){
    g_assert(Lookahead_pos_is_set(lookahead, 0));
    Lookahead_clear(lookahead, 0);
    if(lookahead->mask) /* not empty */
        lookahead->pos += Lookahead_next_pos(lookahead);
    else
        lookahead->pos = lookahead->reset_pos; /* empty */
    return lookahead->pos;
    }

void Lookahead_move(Lookahead *lookahead, gint pos){
    register gint next_advance, max_advance = pos - lookahead->pos;
    g_assert(pos >= lookahead->pos);
    g_assert(max_advance >= 0);
    while(lookahead->mask){
        next_advance = Lookahead_next_pos(lookahead);
        g_assert(next_advance >= 0);
        g_assert(next_advance <= lookahead->max_advance);
        if(next_advance < max_advance){
            Lookahead_clear(lookahead, next_advance);
        } else {
            break;
            }
        }
    lookahead->pos = pos; /* update the position */
    return;
    }

/**/

