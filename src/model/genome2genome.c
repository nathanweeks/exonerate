/****************************************************************\
*                                                                *
*  Genome <-> Genome comparison model                            *
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

#include "genome2genome.h"

Genome2Genome_Data *Genome2Genome_Data_create(
                         Sequence *query, Sequence *target){
    register Genome2Genome_Data *g2gd = g_new0(Genome2Genome_Data, 1);
    g_assert(query->alphabet->type == Alphabet_Type_DNA);
    g_assert(target->alphabet->type == Alphabet_Type_DNA);
    CDNA2Genome_Data_init(&g2gd->cd2gd, query, target);
    return g2gd;
    }

void Genome2Genome_Data_destroy(Genome2Genome_Data *g2gd){
    CDNA2Genome_Data_clear(&g2gd->cd2gd);
    g_free(g2gd);
    return;
    }

/**/

C4_Model *Genome2Genome_create(void){
    register C4_Model *model = C4_Model_create("genome2genome");
    register C4_Model *cdna2genome = CDNA2Genome_create();
    register GPtrArray *transition_list;
    register gint i;
    register C4_Transition *transition;
    register C4_Model *query_intron_model, *joint_intron_model,
                      *query_phase_model, *joint_phase_model;
    register Match *codon_match;
    C4_Model_insert(model, cdna2genome, NULL, NULL);
    C4_Model_destroy(cdna2genome);
    /**/
    query_intron_model = Intron_create("query", TRUE, FALSE, TRUE);
    joint_intron_model = Intron_create("joint", TRUE, TRUE, TRUE);
    /**/
    codon_match = Match_find(Match_Type_CODON2CODON);
    query_phase_model = Phase_create("query", codon_match, TRUE, FALSE);
    joint_phase_model = Phase_create("joint", codon_match, TRUE, TRUE);
    /**/
    transition_list = C4_Model_select_transitions(model,
                                                  C4_Label_MATCH);
    g_assert(transition_list);
    g_assert(transition_list->len == 3);
    for(i = 0; i < transition_list->len; i++){
        transition = transition_list->pdata[i];
        if((transition->advance_query == 1)
        && (transition->advance_query == 1)){
            C4_Model_insert(model, query_intron_model,
                            transition->input, transition->output);
            C4_Model_insert(model, joint_intron_model,
                            transition->input, transition->output);
        } else {
            g_assert(transition->advance_query == 3);
            g_assert(transition->advance_target == 3);
            C4_Model_insert(model, query_phase_model,
                            transition->input, transition->output);
            C4_Model_insert(model, joint_phase_model,
                            transition->input, transition->output);
            }
        }
    C4_Model_destroy(query_intron_model);
    C4_Model_destroy(joint_intron_model);
    C4_Model_destroy(query_phase_model);
    C4_Model_destroy(joint_phase_model);
    g_ptr_array_free(transition_list, TRUE);
    C4_Model_close(model);
    return model;
    }

