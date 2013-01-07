/****************************************************************\
*                                                                *
*  Interface for different types of alignment model              *
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

#ifndef INCLUDED_MODEL_TYPE_H
#define INCLUDED_MODEL_TYPE_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "c4.h"
#include "alphabet.h"
#include "sequence.h"

typedef enum {
    Model_Type_UNGAPPED,
    Model_Type_UNGAPPED_TRANS,
    Model_Type_AFFINE_GLOBAL,
    Model_Type_AFFINE_BESTFIT,
    Model_Type_AFFINE_LOCAL,
    Model_Type_AFFINE_OVERLAP,
    Model_Type_EST2GENOME,
    Model_Type_NER,
    Model_Type_PROTEIN2DNA,
    Model_Type_PROTEIN2DNA_BESTFIT,
    Model_Type_PROTEIN2GENOME,
    Model_Type_PROTEIN2GENOME_BESTFIT,
    Model_Type_CODING2CODING,
    Model_Type_CODING2GENOME,
    Model_Type_CDNA2GENOME,
    Model_Type_GENOME2GENOME,
    Model_Type_TOTAL /* just to store the total */
} Model_Type;

gchar *Model_Type_to_string(Model_Type type);
Model_Type Model_Type_from_string(gchar *str);

gboolean Model_Type_is_gapped(Model_Type type);
gboolean Model_Type_translate_both(Model_Type type);
gboolean Model_Type_has_dual_match(Model_Type type);
gboolean Model_Type_has_genomic_target(Model_Type type);

C4_Model *Model_Type_get_model(Model_Type type,
                               Alphabet_Type query_type,
                               Alphabet_Type target_type);

gpointer Model_Type_create_data(Model_Type type,
                                Sequence *query, Sequence *target);
void Model_Type_destroy_data(Model_Type type, gpointer model_data);

/**/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* INCLUDED_MODEL_TYPE_H */

