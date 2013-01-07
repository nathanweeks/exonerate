/****************************************************************\
*                                                                *
*  cDNA <-> genomic alignment model                              *
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

#include "cdna2genome.h"

/**/

void CDNA2Genome_Data_init(CDNA2Genome_Data *cd2gd,
                           Sequence *query, Sequence *target){
    g_assert(query->alphabet->type == Alphabet_Type_DNA);
    g_assert(target->alphabet->type == Alphabet_Type_DNA);
    Coding2Genome_Data_init(&cd2gd->c2gd, query, target);
    EST2Genome_Data_init(&cd2gd->e2gd, query, target);
    return;
    }

CDNA2Genome_Data *CDNA2Genome_Data_create(Sequence *query,
                                          Sequence *target){
    register CDNA2Genome_Data *cd2gd = g_new0(CDNA2Genome_Data, 1);
    CDNA2Genome_Data_init(cd2gd, query, target);
    return cd2gd;
    }

void CDNA2Genome_Data_clear(CDNA2Genome_Data *cd2gd){
    Coding2Genome_Data_clear(&cd2gd->c2gd);
    EST2Genome_Data_clear(&cd2gd->e2gd);
    return;
    }

void CDNA2Genome_Data_destroy(CDNA2Genome_Data *cd2gd){
    CDNA2Genome_Data_clear(cd2gd);
    g_free(cd2gd);
    return;
    }

/**/

/* Create a special single-orientation e2g model */
/* FIXME: move this to utr module ? */
static C4_Model *CDNA2Genome_UTR_create(void){
    register C4_Model *model = Affine_create(
                Affine_Model_Type_LOCAL,
                Alphabet_Type_DNA,
                Alphabet_Type_DNA, FALSE);
    register C4_Model *intron_model = Intron_create("forward",
                                                    FALSE, TRUE, TRUE);
    register C4_Transition *match_transition
           = C4_Model_select_single_transition(model, C4_Label_MATCH);
    g_assert(match_transition);
    C4_Model_open(model);
    C4_Model_insert(model, intron_model, match_transition->input,
                                         match_transition->output);
    C4_Model_close(model);
    C4_Model_destroy(intron_model);
    return model;
    }

C4_Model *CDNA2Genome_create(void){
    register C4_Model *model = C4_Model_create("cdna2genome");
    register C4_Model *c2g_model = Coding2Genome_create();
    register C4_Model *utr_model = CDNA2Genome_UTR_create();
    register C4_Transition *transition;
    g_assert(model);
    /**/
    C4_Model_insert(model, c2g_model, NULL, NULL);
    transition = C4_Model_select_single_transition(model,
                                                   C4_Label_MATCH);
    g_assert(transition);
    g_assert(transition->input == transition->output);
    g_assert(transition->label_data);
    g_assert(((Match*)transition->label_data)->type
            == Match_Type_CODON2CODON);
    C4_Model_insert(model, utr_model, NULL, transition->input);
    C4_Model_insert(model, utr_model, transition->input, NULL);
    /* FIXME: add transition start->post_utr || post_utr->end */
    C4_Model_destroy(c2g_model);
    C4_Model_destroy(utr_model);
    C4_Model_close(model);
    /* C4_Model_dump_graphviz(model); */
    return model;
    }

/**/

