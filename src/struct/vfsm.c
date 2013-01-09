/****************************************************************\
*                                                                *
*  VFSM Library : complete trie navigation                       *
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

/* Based on Virtual Finite State Machine routines
 * from the ESTate package (by the same author).
 * http://www.hgmp.mrc.ac.uk/~gslater/ESTate/
 */

#include <string.h>  /* For strlen() */
#include <strings.h> /* For ffs() */
#include <math.h>    /* For pow() */

#include "vfsm.h"

#ifndef G_GNUC_EXTENSION
#define G_GNUC_EXTENSION
#endif /* G_GNUC_EXTENSION */

VFSM *VFSM_create(gchar *alphabet, guint depth){
    register VFSM *vfsm = g_new0(VFSM, 1);
    register gint i;
    register VFSM_Int tmp;
    g_assert(depth);
    g_assert(alphabet);
    vfsm->alphabet_size = strlen(alphabet);
    if(pow(vfsm->alphabet_size, depth) >= VFSM_MAX){
        g_free(vfsm);
        return NULL;
        }
    vfsm->alphabet  = g_strdup(alphabet);
    vfsm->depth     = depth;
    g_assert(vfsm->depth >= 3); /* FIXME: fix for depth < 3 */
    tmp = vfsm->prs = vfsm->alphabet_size;
    for(i = 3; i < depth; i++){
        tmp *= vfsm->alphabet_size;
        vfsm->prs += tmp;
        }
    vfsm->prs++;
    vfsm->lrs   = vfsm->prs + (tmp * vfsm->alphabet_size);
    vfsm->prw   = vfsm->lrs - vfsm->prs;
    vfsm->lrw   = vfsm->prw * vfsm->alphabet_size;
    vfsm->total = vfsm->lrw + vfsm->lrs;
    if(vfsm->alphabet_size & (vfsm->alphabet_size-1)){
        vfsm->is_poweroftwo = FALSE;
        vfsm->log_alphabet_size = 0;
    } else {
        vfsm->is_poweroftwo = TRUE;
        vfsm->log_alphabet_size = ffs(vfsm->alphabet_size)-1;
        }
    for(i = 0; i < vfsm->alphabet_size; i++)
        vfsm->index[(guchar)vfsm->alphabet[i]] = i+1;
    vfsm->branch_size = g_new(VFSM_Int, depth);
    tmp = 1;
    for(i = 0; i < depth; i++){
        vfsm->branch_size[depth-i-1] = tmp;
        tmp *= vfsm->alphabet_size;
        }
    return vfsm;
    }

void VFSM_destroy(VFSM *vfsm){
    g_free(vfsm->branch_size);
    g_free(vfsm->alphabet);
    g_free(vfsm);
    return;
    }

// allow compilation with "-Werror"
#pragma GCC diagnostic warning "-Wformat"
void VFSM_info(VFSM *vfsm){
    register gint i;
    G_GNUC_EXTENSION /* Allow %llu without pedantic warning */
    g_message("VFSM info:\n"
           " alphabet          = [%s]\n"
           " alphabet_size     = %d\n"
           " depth             = %d\n"
           " is_poweroftwo     = [%s]\n"
           " log_alphabet_size = %d\n"
           " prs           = %" VFSM_PRINT_FORMAT "\n"
           " prw           = %" VFSM_PRINT_FORMAT "\n"
           " lrs           = %" VFSM_PRINT_FORMAT "\n"
           " lrw           = %" VFSM_PRINT_FORMAT "\n"
           " total         = %" VFSM_PRINT_FORMAT "\n",
           vfsm->alphabet,
           vfsm->alphabet_size,
           vfsm->depth,
           vfsm->is_poweroftwo?"TRUE":"FALSE",
           vfsm->log_alphabet_size,
           vfsm->prs,
           vfsm->prw,
           vfsm->lrs,
           vfsm->lrw,
           vfsm->total
           );
    for(i = 0; i < vfsm->depth; i++){
        G_GNUC_EXTENSION /* Allow %llu without pedantic warning */
        g_print("  branch_size [%d] [%" VFSM_PRINT_FORMAT "]\n",
                i, vfsm->branch_size[i]);
        }
    g_print("--\n");
    return;
    }

static gboolean VFSM_word_is_valid_leaf(VFSM *vfsm, gchar *word){
    register gint i;
    for(i = 0; word[i]; i++)
        if(!vfsm->index[(guchar)word[i]])
            return FALSE;
    if(i != vfsm->depth)
        return FALSE;
    return TRUE;
    }

VFSM_Int VFSM_word2state(VFSM *vfsm, gchar *word){
    register gint i;
    register VFSM_Int pos;
    g_assert(VFSM_word_is_valid_leaf(vfsm, word));
    if(vfsm->is_poweroftwo)
        for(i = 0, pos = 0; i < vfsm->depth; i++){
            pos <<= vfsm->log_alphabet_size;
            pos |= vfsm->index[(guchar)word[i]]-1;
            }
    else
        for(i = 0, pos = 0; i < vfsm->depth; i++){
            pos *= vfsm->alphabet_size;
            pos += vfsm->index[(guchar)word[i]]-1;
            }
    return VFSM_leaf2state(vfsm, pos);
    }

gboolean VFSM_state2word(VFSM *vfsm, VFSM_Int state, gchar *word){
    register gint i, j;
    g_assert(word);
    if(!VFSM_state_is_leaf(vfsm, state))
        return FALSE;
    state = VFSM_state2leaf(vfsm, state);
    if(vfsm->is_poweroftwo)
        for(i = vfsm->depth-1, j = vfsm->alphabet_size-1; i >= 0; i--){
            word[i] = vfsm->alphabet[state & j];
            state >>= vfsm->log_alphabet_size;
            }
    else
        for(i = vfsm->depth-1; i >= 0; i--){
            word[i] = vfsm->alphabet[state % vfsm->alphabet_size];
            state /= vfsm->alphabet_size;
            }
    word[vfsm->depth] = '\0';
    return TRUE;
    }

VFSM_Int VFSM_change_state(VFSM *vfsm, VFSM_Int state, guchar ch){
    if(!vfsm->index[ch])
        return 0; /* Not a valid character: reset VFSM */
    if(VFSM_state_is_leaf(vfsm, state)) /* -> longest suffix */
        state = vfsm->prs + ((state-vfsm->lrs) % vfsm->prw);
    return (state * vfsm->alphabet_size)
         + vfsm->index[ch]; /* Descend */
    }

VFSM_Int VFSM_change_state_POW2(VFSM *vfsm, VFSM_Int state, guchar ch){
    g_assert(vfsm->is_poweroftwo);
    if(!vfsm->index[ch])
        return 0; /* NOT A VALID CHARACTER: RESET VFSM */
    if(VFSM_state_is_leaf(vfsm, state)) /* -> longest suffix */
        state = vfsm->prs + ((state-vfsm->lrs) & (vfsm->prw-1));
/* As alphabet_size is a power of 2, so will be vfsm->prw */
    return (state << vfsm->log_alphabet_size) /* Descend */
         + vfsm->index[ch];
    }

gchar VFSM_symbol_by_pos(VFSM *vfsm, VFSM_Int state, guint depth){
    register guint apos;
    g_assert(depth < vfsm->depth);
    apos = VFSM_state2leaf(vfsm, state) / vfsm->branch_size[depth]
         % vfsm->alphabet_size;
    g_assert(apos < vfsm->alphabet_size);
    return vfsm->alphabet[apos];
    }

gchar VFSM_symbol_by_pos_POW2(VFSM *vfsm, VFSM_Int state, guint depth){
    register guint apos;
    g_assert(depth < vfsm->depth);
    apos = (VFSM_state2leaf(vfsm, state)
         >> ((vfsm->depth-depth-1) << (vfsm->log_alphabet_size-1))
         & (vfsm->alphabet_size-1));
    g_assert(apos < vfsm->alphabet_size);
    return vfsm->alphabet[apos];
    }

VFSM_Int VFSM_jump(VFSM *vfsm, VFSM_Int state,
                   guint depth, guchar next){
    register guchar curr;
    g_assert(vfsm);
    curr = VFSM_symbol_by_pos(vfsm, state, depth);
    return state + (vfsm->branch_size[depth]
           * (vfsm->index[next]-vfsm->index[curr]));
    }

