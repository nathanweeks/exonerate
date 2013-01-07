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

#include "vfsm.h"

int main(void){
    register gint i, j;
    register VFSM *vfsm = VFSM_create("ACGT", 12);
    register guchar *seq = (guchar*)
                 "CGATCGATCTGATCGTAGNTAGCTCGATCGATGNAGCTAGC";
    register VFSM_Int state = 0;
    register gchar *testword;
    gchar word[13];
    VFSM_info(vfsm);
    for(i = j = 0; seq[i]; i++){
        if(!vfsm->index[seq[i]]){
            state = 0;
            continue;
            }
        state = VFSM_change_state_MPOW2(vfsm, state, seq[i]);
        if(VFSM_state2word(vfsm, state, word))
            g_message("pos=[%2d] word_number=[%d] word=[%s]",
                    i, ++j, word);
        }
    testword = "ACGTAACCGGTT";
    g_message("Using test word [%s]", testword);
    state = VFSM_word2state(vfsm, testword);
    g_message("Gives state [%d]", (gint)state);
    state = VFSM_jump(vfsm, state, 2, 'A');
    g_message("Jump to state [%d]", (gint)state);
    g_message("Convert status [%s]",
            VFSM_state2word(vfsm, state, word)?"true":"false");
    g_message("Jump 2:G->A [%s]", word);
    VFSM_destroy(vfsm);
    return 0;
    }

