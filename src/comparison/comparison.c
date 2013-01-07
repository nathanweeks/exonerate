/****************************************************************\
*                                                                *
*  Comparison : A module for pairwise sequence comparisons       *
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

#include "comparison.h"

/**/

static HSP_Param *Comparison_Param_get_hsp_param(HSP_Param *hsp_param,
                                                 gboolean swap_param){
    if(hsp_param){
        if(swap_param)
            return HSP_Param_swap(hsp_param);
        else
            return HSP_Param_share(hsp_param);
        }
    return NULL;
    }

static Comparison_Param *Comparison_Param_create_without_mirror(
                                          Alphabet_Type query_type,
                                          Alphabet_Type target_type,
                                          HSP_Param *dna_hsp_param,
                                          HSP_Param *protein_hsp_param,
                                          HSP_Param *codon_hsp_param,
                                          gboolean swap_param){
    register Comparison_Param *param = g_new(Comparison_Param, 1);
    param->thread_ref = NULL;
    param->query_type = query_type;
    param->target_type = target_type;
    param->dna_hsp_param = Comparison_Param_get_hsp_param(dna_hsp_param,
                                                          swap_param);
    param->protein_hsp_param = Comparison_Param_get_hsp_param(protein_hsp_param,
                                                              swap_param);
    param->codon_hsp_param = Comparison_Param_get_hsp_param(codon_hsp_param,
                                                            swap_param);
    param->mirror = NULL;
    return param;
    }

Comparison_Param *Comparison_Param_create(Alphabet_Type query_type,
                                          Alphabet_Type target_type,
                                          HSP_Param *dna_hsp_param,
                                          HSP_Param *protein_hsp_param,
                                          HSP_Param *codon_hsp_param){
    register Comparison_Param *param
        = Comparison_Param_create_without_mirror(
                            query_type, target_type,
                            dna_hsp_param, protein_hsp_param,
                            codon_hsp_param, FALSE);
    param->mirror = Comparison_Param_create_without_mirror(
                            target_type, query_type,
                            dna_hsp_param, protein_hsp_param,
                            codon_hsp_param, TRUE);
    param->mirror->mirror = param;
    param->mirror->thread_ref = param->thread_ref = ThreadRef_create();
    return param;
    }

static void Comparison_Param_destroy_without_mirror(
            Comparison_Param *param){
    g_assert(param);
    param->mirror = NULL;
    if(param->dna_hsp_param)
        HSP_Param_destroy(param->dna_hsp_param);
    if(param->protein_hsp_param)
        HSP_Param_destroy(param->protein_hsp_param);
    if(param->codon_hsp_param)
        HSP_Param_destroy(param->codon_hsp_param);
    g_free(param);
    return;
    }

void Comparison_Param_destroy(Comparison_Param *param){
    g_assert(param);
    if(ThreadRef_destroy(param->thread_ref))
        return;
    Comparison_Param_destroy_without_mirror(param->mirror);
    Comparison_Param_destroy_without_mirror(param);
    return;
    }

Comparison_Param *Comparison_Param_share(Comparison_Param *param){
    g_assert(param);
    ThreadRef_share(param->thread_ref);
    return param;
    }

Comparison_Param *Comparison_Param_swap(Comparison_Param *param){
    register Comparison_Param *mirror
           = Comparison_Param_share(param->mirror);
    g_assert(param);
    g_assert(param->mirror);
    Comparison_Param_destroy(param);
    return mirror;
    }

HSPset_ArgumentSet *Comparison_Param_get_HSPSet_Argument_Set(
                    Comparison_Param *param){
    if(param->dna_hsp_param)
        return param->dna_hsp_param->has;
    if(param->protein_hsp_param)
        return param->protein_hsp_param->has;
    if(param->codon_hsp_param)
        return param->codon_hsp_param->has;
    g_error("Comparison_Param does not contain any HSP_Param");
    return NULL;
    }

/**/

Comparison *Comparison_create(Comparison_Param *param,
                              Sequence *query, Sequence *target){
    register Comparison *comparison = g_new(Comparison, 1);
    g_assert(param);
    g_assert(query);
    g_assert(target);
    comparison->thread_ref = ThreadRef_create();
    comparison->param = Comparison_Param_share(param);
    comparison->query = Sequence_share(query);
    comparison->target = Sequence_share(target);
    ThreadRef_lock(comparison->param->thread_ref);
    comparison->dna_hspset = comparison->param->dna_hsp_param
                           ? HSPset_create(query, target,
                                           param->dna_hsp_param)
                           : NULL;
    comparison->protein_hspset = param->protein_hsp_param
                               ? HSPset_create(query, target,
                                               param->protein_hsp_param)
                               : NULL;
    comparison->codon_hspset = param->codon_hsp_param
                             ? HSPset_create(query, target,
                                             param->codon_hsp_param)
                             : NULL;
    ThreadRef_unlock(comparison->param->thread_ref);
    return comparison;
    }

Comparison *Comparison_share(Comparison *comparison){
    ThreadRef_share(comparison->thread_ref);
    return comparison;
    }

void Comparison_destroy(Comparison *comparison){
    g_assert(comparison);
    if(ThreadRef_destroy(comparison->thread_ref))
        return;
    Comparison_Param_destroy(comparison->param);
    Sequence_destroy(comparison->query);
    Sequence_destroy(comparison->target);
    if(comparison->dna_hspset)
        HSPset_destroy(comparison->dna_hspset);
    if(comparison->protein_hspset)
        HSPset_destroy(comparison->protein_hspset);
    if(comparison->codon_hspset)
        HSPset_destroy(comparison->codon_hspset);
    g_free(comparison);
    return;
    }

void Comparison_print(Comparison *comparison){
    g_message("Comparison qy [%s] (qylen:%d) tg [%s] (tglen:%d)",
            comparison->query->id, comparison->query->len,
            comparison->target->id, comparison->target->len);
    if(comparison->dna_hspset){
        g_message("DNA hspset");
        HSPset_print(comparison->dna_hspset);
        }
    if(comparison->protein_hspset){
        g_message("Protein hspset");
        HSPset_print(comparison->protein_hspset);
        }
    if(comparison->codon_hspset){
        g_message("Codon hspset");
        HSPset_print(comparison->codon_hspset);
        }
    return;
    }

gboolean Comparison_has_hsps(Comparison *comparison){
    if(comparison->dna_hspset)
        if(!HSPset_is_empty(comparison->dna_hspset))
            return TRUE;
    if(comparison->protein_hspset)
        if(!HSPset_is_empty(comparison->protein_hspset))
            return TRUE;
    if(comparison->codon_hspset)
        if(!HSPset_is_empty(comparison->codon_hspset))
            return TRUE;
    return FALSE;
    }

void Comparison_finalise(Comparison *comparison){
    if(comparison->dna_hspset)
        comparison->dna_hspset = HSPset_finalise(comparison->dna_hspset);
    if(comparison->protein_hspset)
        comparison->protein_hspset = HSPset_finalise(comparison->protein_hspset);
    if(comparison->codon_hspset)
        comparison->codon_hspset = HSPset_finalise(comparison->codon_hspset);
    return;
    }

void Comparison_swap(Comparison *comparison){
    register Sequence *query;
    g_assert(ThreadRef_get_count(comparison->thread_ref) == 1);
    /* Mirror parameters */
    comparison->param = Comparison_Param_share(comparison->param->mirror);
    Comparison_Param_destroy(comparison->param->mirror);
    /* Swap query and target */
    query = comparison->query;
    comparison->query = comparison->target;
    comparison->target = query;
    /* Swap each HSPset */
    if(comparison->dna_hspset)
        HSPset_swap(comparison->dna_hspset,
                    comparison->param->dna_hsp_param);
    if(comparison->protein_hspset)
        HSPset_swap(comparison->protein_hspset,
                    comparison->param->protein_hsp_param);
    if(comparison->codon_hspset)
        HSPset_swap(comparison->codon_hspset,
                    comparison->param->codon_hsp_param);
    return;
    }
/* Swaps query and target of comparison in place */

void Comparison_revcomp(Comparison *comparison){
    register Sequence *rc_query = Sequence_revcomp(comparison->query),
                      *rc_target = Sequence_revcomp(comparison->target);
    g_assert(ThreadRef_get_count(comparison->thread_ref) == 1);
    g_assert(comparison->dna_hspset);
    g_assert(!comparison->protein_hspset);
    g_assert(!comparison->codon_hspset);
    Sequence_destroy(comparison->query);
    Sequence_destroy(comparison->target);
    comparison->query = rc_query;
    comparison->target = rc_target;
    HSPset_revcomp(comparison->dna_hspset);
    return;
    }
/* Revcomps the comparison in place */

/**/

