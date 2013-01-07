/****************************************************************\
*                                                                *
*  Module for EST <-> genome alignments                          *
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

#include <ctype.h>  /* For tolower() */
#include <string.h> /* For strlen()  */

#include "est2genome.h"

void EST2Genome_Data_init(EST2Genome_Data *e2gd,
                          Sequence *query, Sequence *target){
    g_assert(e2gd);
    g_assert(query);
    g_assert(target);
    Affine_Data_init(&e2gd->ad, query, target, FALSE);
    if(!EST2Genome_Data_get_Intron_Data(e2gd))
        EST2Genome_Data_get_Intron_Data(e2gd) = Intron_Data_create();
    return;
    }

EST2Genome_Data *EST2Genome_Data_create(Sequence *query,
                                        Sequence *target){
    register EST2Genome_Data *e2gd = g_new0(EST2Genome_Data, 1);
    EST2Genome_Data_init(e2gd, query, target);
    return e2gd;
    }

void EST2Genome_Data_clear(EST2Genome_Data *e2gd){
    g_assert(e2gd);
    if(EST2Genome_Data_get_Intron_Data(e2gd)){
        Intron_Data_destroy(EST2Genome_Data_get_Intron_Data(e2gd));
        EST2Genome_Data_get_Intron_Data(e2gd) = NULL;
        }
    Affine_Data_clear(&e2gd->ad);
    return;
    }

void EST2Genome_Data_destroy(EST2Genome_Data *e2gd){
    g_assert(e2gd);
    EST2Genome_Data_clear(e2gd);
    g_free(e2gd);
    return;
    }

/**/

C4_Model *EST2Genome_create(void){
    register C4_Model *est2genome = Affine_create(
                                    Affine_Model_Type_LOCAL,
                                    Alphabet_Type_DNA,
                                    Alphabet_Type_DNA, FALSE);
    register C4_Model *intron_forward_model,
                      *intron_reverse_model;
    register C4_Transition *match_forward_transition,
                           *match_reverse_transition;
    register GPtrArray *match_transition_list;
    C4_Model_rename(est2genome, "est2genome");
    C4_Model_open(est2genome);
    C4_Model_make_stereo(est2genome, "forward", "reverse");
    match_transition_list = C4_Model_select_transitions(est2genome,
                                                        C4_Label_MATCH);
    g_assert(match_transition_list);
    g_assert(match_transition_list->len == 2);
    match_forward_transition = match_transition_list->pdata[0];
    match_reverse_transition = match_transition_list->pdata[1];
    g_ptr_array_free(match_transition_list, TRUE);
    intron_forward_model = Intron_create("forward", FALSE, TRUE, TRUE);
    intron_reverse_model = Intron_create("reverse", FALSE, TRUE, FALSE);
    C4_Model_insert(est2genome, intron_forward_model,
                    match_forward_transition->input,
                    match_forward_transition->input);
    C4_Model_insert(est2genome, intron_reverse_model,
                    match_reverse_transition->input,
                    match_reverse_transition->input);
    C4_Model_destroy(intron_forward_model);
    C4_Model_destroy(intron_reverse_model);
    C4_Model_append_codegen(est2genome,
     "#include \"est2genome.h\"\n",
     NULL, NULL);
    C4_Model_close(est2genome);
    return est2genome;
    }

