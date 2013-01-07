/****************************************************************\
*                                                                *
*  Protein <-> DNA comparison model                              *
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

#include <string.h> /* For strlen() */

#include "protein2dna.h"
#include "frameshift.h"

void Protein2DNA_Data_init(Protein2DNA_Data *p2dd,
                           Sequence *query, Sequence *target){
    g_assert(query->alphabet->type == Alphabet_Type_PROTEIN);
    g_assert(target->alphabet->type == Alphabet_Type_DNA);
    Affine_Data_init(&p2dd->ad, query, target, FALSE);
    if(!Protein2DNA_Data_get_Frameshift_Data(p2dd))
        Protein2DNA_Data_get_Frameshift_Data(p2dd)
                           = Frameshift_Data_create();
    return;
    }

Protein2DNA_Data *Protein2DNA_Data_create(
                      Sequence *query, Sequence *target){
    register Protein2DNA_Data *p2dd = g_new0(Protein2DNA_Data, 1);
    Protein2DNA_Data_init(p2dd, query, target);
    return p2dd;
    }

void Protein2DNA_Data_clear(Protein2DNA_Data *p2dd){
    Affine_Data_clear(&p2dd->ad);
    if(Protein2DNA_Data_get_Frameshift_Data(p2dd)){
        Frameshift_Data_destroy(
                   Protein2DNA_Data_get_Frameshift_Data(p2dd));
        Protein2DNA_Data_get_Frameshift_Data(p2dd) = NULL;
        }
    return;
    }

void Protein2DNA_Data_destroy(Protein2DNA_Data *p2dd){
    Protein2DNA_Data_clear(p2dd);
    g_free(p2dd);
    return;
    }

C4_Model *Protein2DNA_create(Affine_Model_Type type){
    register gchar *name = g_strdup_printf("protein2dna:%s",
                                           Affine_Model_Type_get_name(type));
    register C4_Model *model = NULL;
    register C4_Transition *match_transition;
    model = Affine_create(type, Alphabet_Type_PROTEIN,
                                Alphabet_Type_DNA, FALSE);
    g_assert(model);
    C4_Model_rename(model, name);
    g_free(name);
    C4_Model_open(model);
    match_transition = C4_Model_select_single_transition(model,
                                                 C4_Label_MATCH);
    g_assert(match_transition);
    Frameshift_add(model, match_transition->input, "p2d", FALSE);
    C4_Model_close(model);
    return model;
    }

