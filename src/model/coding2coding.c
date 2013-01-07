/****************************************************************\
*                                                                *
*  Coding DNA comparison model                                   *
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

#include "coding2coding.h"
#include "frameshift.h"

void Coding2Coding_Data_init(Coding2Coding_Data *c2cd,
                           Sequence *query, Sequence *target){
    g_assert(query->alphabet->type == Alphabet_Type_DNA);
    g_assert(target->alphabet->type == Alphabet_Type_DNA);
    Affine_Data_init(&c2cd->ad, query, target, TRUE);
    if(!Coding2Coding_Data_get_Frameshift_Data(c2cd))
        Coding2Coding_Data_get_Frameshift_Data(c2cd)
              = Frameshift_Data_create();
    return;
    }

Coding2Coding_Data *Coding2Coding_Data_create(
                      Sequence *query, Sequence *target){
    register Coding2Coding_Data *c2cd = g_new0(Coding2Coding_Data, 1);
    Coding2Coding_Data_init(c2cd, query, target);
    return c2cd;
    }

void Coding2Coding_Data_clear(Coding2Coding_Data *c2cd){
    Affine_Data_clear(&c2cd->ad);
    Frameshift_Data_destroy(
               Coding2Coding_Data_get_Frameshift_Data(c2cd));
    return;
    }

void Coding2Coding_Data_destroy(Coding2Coding_Data *c2cd){
    Coding2Coding_Data_clear(c2cd);
    g_free(c2cd);
    return;
    }

C4_Model *Coding2Coding_create(void){
    register gchar *name = "coding2coding";
    register C4_Model *model = NULL;
    register C4_Transition *match_transition;
    model = Affine_create(Affine_Model_Type_LOCAL,
                          Alphabet_Type_DNA, Alphabet_Type_DNA, TRUE);
    g_assert(model);
    C4_Model_rename(model, name);
    C4_Model_open(model);
    match_transition = C4_Model_select_single_transition(model,
                                                       C4_Label_MATCH);
    g_assert(match_transition);
    Frameshift_add(model, match_transition->input, "query", TRUE);
    Frameshift_add(model, match_transition->input, "target", FALSE);
    C4_Model_close(model);
    return model;
    }

