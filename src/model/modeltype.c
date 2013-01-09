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

#include <strings.h> /* For strcasecmp() */
#include "modeltype.h"

#include "ungapped.h"
#include "affine.h"
#include "est2genome.h"
#include "ner.h"
#include "protein2dna.h"
#include "protein2genome.h"
#include "coding2coding.h"
#include "coding2genome.h"
#include "cdna2genome.h"
#include "genome2genome.h"

gchar *Model_Type_to_string(Model_Type type){
    register gchar *name = NULL;
    switch(type){
        case Model_Type_UNGAPPED:
            name = "ungapped";
            break;
        case Model_Type_UNGAPPED_TRANS:
            name = "ungapped:trans";
            break;
        case Model_Type_AFFINE_GLOBAL:
            name = "affine:global";
            break;
        case Model_Type_AFFINE_BESTFIT:
            name = "affine:bestfit";
            break;
        case Model_Type_AFFINE_LOCAL:
            name = "affine:local";
            break;
        case Model_Type_AFFINE_OVERLAP:
            name = "affine:overlap";
            break;
        case Model_Type_EST2GENOME:
            name = "est2genome";
            break;
        case Model_Type_NER:
            name = "ner";
            break;
        case Model_Type_PROTEIN2DNA:
            name = "protein2dna";
            break;
        case Model_Type_PROTEIN2DNA_BESTFIT:
            name = "protein2dna:bestfit";
            break;
        case Model_Type_PROTEIN2GENOME:
            name = "protein2genome";
            break;
        case Model_Type_PROTEIN2GENOME_BESTFIT:
            name = "protein2genome:bestfit";
            break;
        case Model_Type_CODING2CODING:
            name = "coding2coding";
            break;
        case Model_Type_CODING2GENOME:
            name = "coding2genome";
            break;
        case Model_Type_CDNA2GENOME:
            name = "cdna2genome";
            break;
        case Model_Type_GENOME2GENOME:
            name = "genome2genome";
            break;
        default:
            g_error("Unknown Model Type [%d]", type);
            break;
        }
    return name;
    }

Model_Type Model_Type_from_string(gchar *str){
    gchar *name[Model_Type_TOTAL] = {
         "ungapped", "ungapped:trans",
         "affine:global", "affine:bestfit",
         "affine:local", "affine:overlap",
         "est2genome", "ner", "protein2dna", "protein2dna:bestfit",
         "protein2genome", "protein2genome:bestfit",
         "coding2coding", "coding2genome", "cdna2genome",
         "genome2genome"};
    gchar *short_name[Model_Type_TOTAL] = {
         "u", "u:t",
         "a:g", "a:b", "a:l", "a:o",
         "e2g", "ner",
         "p2d", "p2d:b", "p2g", "p2g:b",
         "c2c", "c2g", "cd2g", "g2g"};
    Model_Type type[Model_Type_TOTAL] = {
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
          Model_Type_GENOME2GENOME};
    register gint i;
    for(i = 0; i < Model_Type_TOTAL; i++)
        if(!strcasecmp(name[i], str))
            return type[i];
    for(i = 0; i < Model_Type_TOTAL; i++)
        if(!strcasecmp(short_name[i], str))
            return type[i];
    g_error("Unknown model type [%s]", str);
    return Model_Type_UNGAPPED; /* Not reached */
    }

gboolean Model_Type_is_gapped(Model_Type type){
    if((type == Model_Type_UNGAPPED)
    || (type == Model_Type_UNGAPPED_TRANS))
        return FALSE;
    return TRUE;
    }

gboolean Model_Type_translate_both(Model_Type type){
    if((type == Model_Type_UNGAPPED_TRANS)
    || (type == Model_Type_CODING2CODING)
    || (type == Model_Type_CODING2GENOME)
    || (type == Model_Type_CDNA2GENOME)
    || (type == Model_Type_GENOME2GENOME))
        return TRUE;
    return FALSE;
    }

gboolean Model_Type_has_dual_match(Model_Type type){
    if((type == Model_Type_CDNA2GENOME)
    || (type == Model_Type_GENOME2GENOME))
        return TRUE;
    return FALSE;
    }

gboolean Model_Type_has_genomic_target(Model_Type type){
    if((type == Model_Type_EST2GENOME)
    || (type == Model_Type_PROTEIN2GENOME)
    || (type == Model_Type_PROTEIN2GENOME_BESTFIT)
    || (type == Model_Type_CODING2GENOME)
    || (type == Model_Type_CDNA2GENOME)
    || (type == Model_Type_GENOME2GENOME))
        return TRUE;
    return FALSE;
    }

static void Model_Type_check_input(Model_Type type,
                                   Alphabet_Type query_type,
                                   Alphabet_Type target_type){
    switch(type){
        case Model_Type_UNGAPPED:
            break;
        case Model_Type_UNGAPPED_TRANS:
        case Model_Type_EST2GENOME:
        case Model_Type_CODING2CODING:
        case Model_Type_CODING2GENOME:
        case Model_Type_CDNA2GENOME:
        case Model_Type_GENOME2GENOME:
            if(query_type != Alphabet_Type_DNA)
                g_error("Expected DNA query (not %s) for model [%s]",
                        Alphabet_Type_get_name(query_type),
                        Model_Type_to_string(type));
            if(target_type != Alphabet_Type_DNA)
                g_error("Expected DNA target (not %s) for model [%s]",
                        Alphabet_Type_get_name(target_type),
                        Model_Type_to_string(type));
            break;
        case Model_Type_AFFINE_GLOBAL:
        case Model_Type_AFFINE_BESTFIT:
        case Model_Type_AFFINE_LOCAL:
        case Model_Type_AFFINE_OVERLAP:
        case Model_Type_NER:
            if(query_type != target_type)
                g_error("Expected similar sequence types for model"
                        " [%s] (not %s:%s)",
                       Model_Type_to_string(type),
                       Alphabet_Type_get_name(query_type),
                       Alphabet_Type_get_name(target_type));
            if(query_type == Alphabet_Type_UNKNOWN)
                g_error("Model [%s] cannot use unknown sequence type",
                        Model_Type_to_string(type));
            break;
        case Model_Type_PROTEIN2DNA:
        case Model_Type_PROTEIN2DNA_BESTFIT:
        case Model_Type_PROTEIN2GENOME:
        case Model_Type_PROTEIN2GENOME_BESTFIT:
            /* qy == AA, tg = NT */
            if(query_type != Alphabet_Type_PROTEIN)
                g_error(
                    "Expected protein query (not %s) for model [%s]",
                        Alphabet_Type_get_name(query_type),
                        Model_Type_to_string(type));
            if(target_type != Alphabet_Type_DNA)
                g_error("Expected DNA target (not %s) for model [%s]",
                        Alphabet_Type_get_name(target_type),
                        Model_Type_to_string(type));
            break;
        default:
            g_error("Unknown model type [%s]",
                    Model_Type_to_string(type));
            break;
        }
    return;
    }

C4_Model *Model_Type_get_model(Model_Type type,
                               Alphabet_Type query_type,
                               Alphabet_Type target_type){
    register C4_Model *model = NULL;
    register Match_Type match_type;
    Model_Type_check_input(type, query_type, target_type);
    switch(type){
        case Model_Type_UNGAPPED:
            match_type = Match_Type_find(query_type, target_type,
                                         FALSE);
            model = Ungapped_create(match_type);
            break;
        case Model_Type_UNGAPPED_TRANS:
            match_type = Match_Type_find(query_type, target_type,
                                         TRUE);
            model = Ungapped_create(match_type);
            break;
        case Model_Type_AFFINE_GLOBAL:
            model = Affine_create(Affine_Model_Type_GLOBAL,
                                  query_type, target_type, FALSE);
            break;
        case Model_Type_AFFINE_BESTFIT:
            model = Affine_create(Affine_Model_Type_BESTFIT,
                                  query_type, target_type, FALSE);
            break;
        case Model_Type_AFFINE_LOCAL:
            model = Affine_create(Affine_Model_Type_LOCAL,
                                  query_type, target_type, FALSE);
            break;
        case Model_Type_AFFINE_OVERLAP:
            model = Affine_create(Affine_Model_Type_OVERLAP,
                                  query_type, target_type, FALSE);
            break;
        case Model_Type_EST2GENOME:
            model = EST2Genome_create();
            break;
        case Model_Type_NER:
            model = NER_create(query_type, target_type);
            break;
        case Model_Type_PROTEIN2DNA:
            model = Protein2DNA_create(Affine_Model_Type_LOCAL);
            break;
        case Model_Type_PROTEIN2DNA_BESTFIT:
            model = Protein2DNA_create(Affine_Model_Type_BESTFIT);
            break;
        case Model_Type_PROTEIN2GENOME:
            model = Protein2Genome_create(Affine_Model_Type_LOCAL);
            break;
        case Model_Type_PROTEIN2GENOME_BESTFIT:
            model = Protein2Genome_create(Affine_Model_Type_BESTFIT);
            break;
        case Model_Type_CODING2CODING:
            model = Coding2Coding_create();
            break;
        case Model_Type_CODING2GENOME:
            model = Coding2Genome_create();
            break;
        case Model_Type_CDNA2GENOME:
            model = CDNA2Genome_create();
            break;
        case Model_Type_GENOME2GENOME:
            model = Genome2Genome_create();
            break;
        default:
            g_error("Unknown Model Type [%d]", type);
            break;
        }
    return model;
    }

gpointer Model_Type_create_data(Model_Type type,
                                Sequence *query, Sequence *target){
    register gpointer model_data = NULL;
    register Match_Type match_type;
    switch(type){
        case Model_Type_UNGAPPED:
            match_type = Match_Type_find(query->alphabet->type,
                                         target->alphabet->type,
                                         FALSE);
            model_data = Ungapped_Data_create(query, target,
                                              match_type);
            break;
        case Model_Type_UNGAPPED_TRANS:
            match_type = Match_Type_find(query->alphabet->type,
                                         target->alphabet->type,
                                         TRUE);
            model_data = Ungapped_Data_create(query, target,
                                              match_type);
            break;
        case Model_Type_AFFINE_GLOBAL:
            /*fallthrough*/
        case Model_Type_AFFINE_BESTFIT:
            /*fallthrough*/
        case Model_Type_AFFINE_LOCAL:
            /*fallthrough*/
        case Model_Type_AFFINE_OVERLAP:
            model_data = Affine_Data_create(query, target, FALSE);
            break;
        case Model_Type_EST2GENOME:
            model_data = EST2Genome_Data_create(query, target);
            break;
        case Model_Type_NER:
            model_data = NER_Data_create(query, target);
            break;
        case Model_Type_PROTEIN2DNA:
            /*fallthrough*/
        case Model_Type_PROTEIN2DNA_BESTFIT:
            model_data = Protein2DNA_Data_create(query, target);
            break;
        case Model_Type_PROTEIN2GENOME:
            /*fallthrough*/
        case Model_Type_PROTEIN2GENOME_BESTFIT:
            model_data = Protein2Genome_Data_create(query, target);
            break;
        case Model_Type_CODING2CODING:
            model_data = Coding2Coding_Data_create(query, target);
            break;
        case Model_Type_CODING2GENOME:
            model_data = Coding2Genome_Data_create(query, target);
            break;
        case Model_Type_CDNA2GENOME:
            model_data = CDNA2Genome_Data_create(query, target);
            break;
        case Model_Type_GENOME2GENOME:
            model_data = Genome2Genome_Data_create(query, target);
            break;
        default:
            g_error("Unknown Model Type [%d]", type);
        }
    return model_data;
    }

void Model_Type_destroy_data(Model_Type type, gpointer model_data){
    switch(type){
        case Model_Type_UNGAPPED:
            /*fallthrough*/
        case Model_Type_UNGAPPED_TRANS:
            Ungapped_Data_destroy(model_data);
            break;
        case Model_Type_AFFINE_GLOBAL:
            /*fallthrough*/
        case Model_Type_AFFINE_BESTFIT:
            /*fallthrough*/
        case Model_Type_AFFINE_LOCAL:
            /*fallthrough*/
        case Model_Type_AFFINE_OVERLAP:
            Affine_Data_destroy(model_data);
            break;
        case Model_Type_EST2GENOME:
            EST2Genome_Data_destroy(model_data);
            break;
        case Model_Type_NER:
            NER_Data_destroy(model_data);
            break;
        case Model_Type_PROTEIN2DNA:
            /*fallthrough*/
        case Model_Type_PROTEIN2DNA_BESTFIT:
            Protein2DNA_Data_destroy(model_data);
            break;
        case Model_Type_PROTEIN2GENOME:
            /*fallthrough*/
        case Model_Type_PROTEIN2GENOME_BESTFIT:
            Protein2Genome_Data_destroy(model_data);
            break;
        case Model_Type_CODING2CODING:
            Coding2Coding_Data_destroy(model_data);
            break;
        case Model_Type_CODING2GENOME:
            Coding2Genome_Data_destroy(model_data);
            break;
        case Model_Type_CDNA2GENOME:
            CDNA2Genome_Data_destroy(model_data);
            break;
        case Model_Type_GENOME2GENOME:
            Genome2Genome_Data_destroy(model_data);
            break;
        default:
            g_error("Unknown Model Type [%d]", type);
        }
    return;
    }

