/****************************************************************\
*                                                                *
*  Coding <-> Genome comparison model                            *
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

#include "coding2genome.h"


void Coding2Genome_Data_init(Coding2Genome_Data *c2gd,
                             Sequence *query, Sequence *target){
    g_assert(query->alphabet->type == Alphabet_Type_DNA);
    g_assert(target->alphabet->type == Alphabet_Type_DNA);
    Coding2Coding_Data_init(&c2gd->c2cd, query, target);
    if(!Coding2Genome_Data_get_Intron_Data(c2gd))
        Coding2Genome_Data_get_Intron_Data(c2gd) = Intron_Data_create();
    return;
    }

Coding2Genome_Data *Coding2Genome_Data_create(
                         Sequence *query, Sequence *target){
    register Coding2Genome_Data *c2gd = g_new0(Coding2Genome_Data, 1);
    Coding2Genome_Data_init(c2gd, query, target);
    return c2gd;
    }

void Coding2Genome_Data_clear(Coding2Genome_Data *c2gd){
    Coding2Coding_Data_clear(&c2gd->c2cd);
    return;
    }

void Coding2Genome_Data_destroy(Coding2Genome_Data *c2gd){
    Coding2Genome_Data_clear(c2gd);
    if(Coding2Genome_Data_get_Intron_Data(c2gd)){
        Intron_Data_destroy(Coding2Genome_Data_get_Intron_Data(c2gd));
        Coding2Genome_Data_get_Intron_Data(c2gd) = NULL;
        }
    g_free(c2gd);
    return;
    }

/**/

C4_Model *Coding2Genome_create(void){
    register C4_Model *model = Coding2Coding_create();
    register C4_Model *phase_model;
    register C4_Transition *match_transition;
    register Match *match;
    g_assert(model);
    C4_Model_rename(model, "coding2genome");
    C4_Model_open(model);
    /**/
    match_transition = C4_Model_select_single_transition(model,
                                                C4_Label_MATCH);
    g_assert(match_transition);
    /* Target Introns */
    match = match_transition->label_data;
    phase_model = Phase_create("target intron", match, FALSE, TRUE);
    C4_Model_insert(model, phase_model, match_transition->input,
                                        match_transition->input);
    C4_Model_destroy(phase_model);
    /**/
    C4_Model_close(model);
    return model;
    }

