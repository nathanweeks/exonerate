/****************************************************************\
*                                                                *
*  Seeder : A module for seeding pairwise alignments             *
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

#include "seeder.h"
#include "alphabet.h"
#include <ctype.h>    /* For toupper()  */
#include <string.h>   /* For strlen()   */
#include <strings.h>   /* For strcasecmp()   */
#include <math.h>     /* For pow()      */
#include <stddef.h>   /* For offsetof() */

#ifndef OFFSET_ITEM
#define OFFSET_ITEM(type, offset, instance) \
     (*(type*)((char*)(instance) + (offset)))
#endif /* OFFSET_ITEM */

#ifndef Swap
#define Swap(x,y,temp) ((temp)=(x),(x)=(y),(y)=(temp))
#endif /* Swap */

Seeder_ArgumentSet *Seeder_ArgumentSet_create(Argument *arg){
    register ArgumentSet *as;
    static Seeder_ArgumentSet sas;
    if(arg){
        as = ArgumentSet_create("Alignment Seeding Options");
        ArgumentSet_add_option(as, 'M', "fsmmemory", "Mb",
                "Memory limit for FSM scanning", "256",
                Argument_parse_int, &sas.fsm_memory_limit);
        ArgumentSet_add_option(as, '\0', "forcefsm", "fsm type",
            "Force FSM type ( normal | compact )", "none",
            Argument_parse_string, &sas.force_fsm);
        ArgumentSet_add_option(as, 0, "wordjump", NULL,
            "Jump between query words", "1",
            Argument_parse_int, &sas.word_jump);
        ArgumentSet_add_option(as, 0, "wordambiguity", NULL,
            "Number of ambiguous words to expand", "1",
            Argument_parse_int, &sas.word_ambiguity);
        Argument_absorb_ArgumentSet(arg, as);
        }
    return &sas;
    }

/**/

static Seeder_WordInfo *Seeder_WordInfo_create(Seeder *seeder){
    register Seeder_WordInfo *word_info
           = RecycleBin_alloc(seeder->recycle_wordinfo);
    word_info->seed_list = NULL;
    word_info->neighbour_list = NULL;
    if(seeder->saturate_threshold){
        word_info->match_count = 0;
        word_info->match_mailbox = -1;
        }
    return word_info;
    }

static void Seeder_WordInfo_empty(Seeder *seeder,
                                  Seeder_WordInfo *word_info){
    register Seeder_Seed *seed;
    register Seeder_Neighbour *neighbour;
    while(word_info->seed_list){
        seed = word_info->seed_list;
        word_info->seed_list = seed->next;
        RecycleBin_recycle(seeder->recycle_seed, seed);
        }
    while(word_info->neighbour_list){
        neighbour = word_info->neighbour_list;
        word_info->neighbour_list = neighbour->next;
        RecycleBin_recycle(seeder->recycle_neighbour, neighbour);
        }
    return;
    }

static void Seeder_WordInfo_add_Seed(Seeder *seeder,
                                     Seeder_WordInfo *word_info,
                                     Seeder_Context *context,
                                     Match_Score query_expect,
                                     gint qpos){
    register Seeder_Seed *seed;
    g_assert(word_info);
    if(seeder->saturate_threshold){
        if(!word_info->match_mailbox) /* Already blocked */
            return;
    /* Apply saturate_threshold to frequent query words */
        if(++word_info->match_count > query_expect){
            word_info->match_mailbox = 0; /* block */
            Seeder_WordInfo_empty(seeder, word_info);
            return;
            }
        }
    seed = RecycleBin_alloc(seeder->recycle_seed);
    seed->next = word_info->seed_list;
    seed->query_pos = qpos;
    seed->context = context;
    word_info->seed_list = seed;
    return;
    }

static void Seeder_WordInfo_add_Neighbour(Seeder *seeder,
            Seeder_WordInfo *word_info,
            Seeder_WordInfo *neighbour){
    register Seeder_Neighbour *n;
    g_assert(seeder);
    g_assert(word_info);
    g_assert(neighbour);
    if(seeder->saturate_threshold
    && (!word_info->match_mailbox)) /* Already blocked */
        return;
    n = RecycleBin_alloc(seeder->recycle_neighbour);
    n->next = word_info->neighbour_list;
    n->word_info = neighbour;
    word_info->neighbour_list = n;
    return;
    }

/**/

static Seeder_QueryInfo *Seeder_QueryInfo_create(Sequence *query){
    register Seeder_QueryInfo *query_info = g_new(Seeder_QueryInfo, 1);
    query_info->query = Sequence_share(query);
    query_info->curr_comparison = NULL;
    return query_info;
    }

static void Seeder_QueryInfo_destroy(Seeder_QueryInfo *query_info){
    Sequence_destroy(query_info->query);
    if(query_info->curr_comparison)
        Comparison_destroy(query_info->curr_comparison);
    g_free(query_info);
    return;
    }

/**/

static gpointer Seeder_FSM_merge_func(gpointer a, gpointer b,
                                      gpointer user_data){
    register Seeder_WordInfo *word_info_a = a,
                             *word_info_b = b;
    register Seeder *seeder = user_data;
    /* Assume word_info_b is always empty during merging */
    g_assert(!word_info_b->seed_list);
    g_assert(!word_info_b->neighbour_list);
    RecycleBin_recycle(seeder->recycle_wordinfo, word_info_b);
    return word_info_a;
    }

static gpointer Seeder_FSM_combine_func(gpointer a, gpointer b,
                                        gpointer user_data){
    g_error("Seeder implementation assumes words of same length");
    return NULL;
    }

static Seeder_FSM *Seeder_FSM_create(Seeder *seeder){
    register Seeder_FSM *seeder_fsm = g_new(Seeder_FSM, 1);
    register Match *match = seeder->any_hsp_param->match;
    seeder_fsm->fsm
            = FSM_create((gchar*)match->comparison_alphabet->member,
                                 Seeder_FSM_merge_func,
                                 Seeder_FSM_combine_func,
                                 seeder);
    return seeder_fsm;
    }

static void Seeder_FSM_destroy(Seeder_FSM *seeder_fsm){
    FSM_destroy(seeder_fsm->fsm);
    g_free(seeder_fsm);
    return;
    }

static gsize Seeder_FSM_memory_usage(Seeder_FSM *seeder_fsm){
    return sizeof(Seeder_FSM)
         + FSM_memory_usage(seeder_fsm->fsm);
    }

typedef struct {
             Seeder *seeder;
    Seeder_WordInfo *curr_word_info;
     Seeder_Context *context;
} Seeder_TraverseData;

static Seeder_WordInfo *Seeder_add_WordInfo(Seeder *seeder, gchar *word,
                                            Seeder_Context *context){
    register Seeder_WordInfo *word_info;
    register VFSM_Int state = 0, leaf;
    if(seeder->seeder_fsm){ /* Add to FSM */
        word_info = Seeder_WordInfo_create(seeder);
        word_info = FSM_add(seeder->seeder_fsm->fsm, word,
                            context->loader->hsp_param->wordlen, word_info);
    } else { /* Add to VFSM */
        g_assert(seeder->seeder_vfsm);
        state = VFSM_word2state(seeder->seeder_vfsm->vfsm, word);
        leaf = VFSM_state2leaf(seeder->seeder_vfsm->vfsm, state);
        word_info = seeder->seeder_vfsm->leaf[leaf];
        if(!word_info){
            word_info = Seeder_WordInfo_create(seeder);
            seeder->seeder_vfsm->leaf[leaf] = word_info;
            }
        }
    return word_info;
    }

static gboolean Seeder_word_is_valid(Match *match,
                                     Sequence *seq, gint pos,
                                     gint len){
    if(seq->annotation){
        if(match->type == Match_Type_DNA2DNA){
            if(((pos+len) > seq->annotation->cds_start)
            && (pos < (seq->annotation->cds_start
                      +seq->annotation->cds_length))){
                return FALSE;
                }
        } else {
            if(match->type == Match_Type_CODON2CODON){
                if((pos < seq->annotation->cds_start)
                || ((pos+len) >= (seq->annotation->cds_start
                                 +seq->annotation->cds_length))
                || ((pos % 3) != (seq->annotation->cds_start % 3))){
                    return FALSE;
                    }
                }
            }
        }
    return TRUE;
    }

/**/

static Seeder_VFSM *Seeder_VFSM_create(Seeder *seeder){
    register Seeder_VFSM *seeder_vfsm = g_new(Seeder_VFSM, TRUE);
    seeder_vfsm->vfsm = VFSM_create(
      (gchar*)seeder->any_hsp_param->match->comparison_alphabet->member,
              seeder->any_hsp_param->wordlen);
    if(!seeder_vfsm->vfsm)
        g_error("Could not create VFSM for alphabet [%s] depth [%d]",
          seeder->any_hsp_param->match->comparison_alphabet->member,
          seeder->any_hsp_param->wordlen);
    seeder_vfsm->leaf = g_new0(Seeder_WordInfo*, seeder_vfsm->vfsm->lrw);
    return seeder_vfsm;
    }

static void Seeder_VFSM_destroy(Seeder_VFSM *seeder_vfsm){
    VFSM_destroy(seeder_vfsm->vfsm);
    g_free(seeder_vfsm->leaf);
    g_free(seeder_vfsm);
    return;
    }

static gsize Seeder_VFSM_memory_usage(Seeder_VFSM *seeder_vfsm){
    return sizeof(Seeder_VFSM)
         + sizeof(VFSM)
         + (sizeof(Seeder_WordInfo*) * seeder_vfsm->vfsm->lrw);
    }

/**/

static gboolean Seeder_decide_fsm_type(gchar *force_fsm){
    if(!strcasecmp(force_fsm, "none"))
        return TRUE;
    /* FIXME: improve this by checking if VFSM is a better idea
     *        by looking at the query sizes ...
     */
    if(!strcasecmp(force_fsm, "normal"))
        return TRUE;
    if(!strcasecmp(force_fsm, "compact"))
        return FALSE;
    g_error("Unknown fsm type [%s]", force_fsm);
    return TRUE;
    }

static Seeder_Loader *Seeder_Loader_create(HSP_Param *hsp_param,
        Seeder_ArgumentSet *sas, gsize hspset_offset){
    register Seeder_Loader *loader = g_new(Seeder_Loader, 1);
    register Alphabet *comparison_alphabet
                      = hsp_param->match->comparison_alphabet;
    register gint alphabet_size
              = strlen((gchar*)comparison_alphabet->member);
    loader->hsp_param = HSP_Param_share(hsp_param);
    if(hsp_param->match->target->is_translated)
        loader->tpos_modifier = (hsp_param->wordlen * 3) - 3;
    else
        loader->tpos_modifier = hsp_param->wordlen - 1;
    loader->saturate_expectation = 1.0 / pow(alphabet_size,
                                             hsp_param->wordlen);
    loader->hspset_offset = hspset_offset;
    return loader;
    }

static void Seeder_Loader_destroy(Seeder_Loader *loader){
    HSP_Param_destroy(loader->hsp_param);
    g_free(loader);
    return;
    }

Seeder *Seeder_create(gint verbosity,
                      Comparison_Param *comparison_param,
                      Match_Score saturate_threshold,
                      Seeder_ReportFunc report_func,
                      gpointer user_data){
    register Seeder *seeder = g_new0(Seeder, 1);
    seeder->sas = Seeder_ArgumentSet_create(NULL);
    if(seeder->sas->word_ambiguity < 1)
        g_error("Word ambiguity cannot be less than one.");
    if(seeder->sas->word_ambiguity > 1){
        if(comparison_param->protein_hsp_param)
            g_error("Protein ambuiguity symbols not implemented");
        g_warning("setting --wordambiguity in exonerate will be SLOW.");
        g_warning("... use esd2esi --wordambiguity and exonerate-server instead");
        }
    seeder->is_prepared = FALSE;
    seeder->verbosity = verbosity;
    seeder->report_func = report_func;
    seeder->user_data = user_data;
    seeder->total_query_length = 0;
    seeder->comparison_count = 0;
    seeder->saturate_threshold = saturate_threshold;
    seeder->comparison_param = Comparison_Param_share(comparison_param);
    seeder->query_info_list = g_ptr_array_new();
    if(saturate_threshold)
        seeder->recycle_wordinfo = RecycleBin_create(
                "Seeder_WordInfo", sizeof(Seeder_WordInfo), 64);
    else
        seeder->recycle_wordinfo = RecycleBin_create(
                "Seeder_WordInfo", sizeof(Seeder_WordInfo_no_ST), 64);
    seeder->recycle_seed = RecycleBin_create(
            "Seeder_Seed", sizeof(Seeder_Seed), 64);
    seeder->recycle_neighbour = RecycleBin_create(
            "Seeder_Neighbour", sizeof(Seeder_Neighbour), 64);
    seeder->recycle_context = RecycleBin_create(
            "Seeder_Context", sizeof(Seeder_Context), 64);
    /**/
    seeder->any_hsp_param = NULL;
    if(comparison_param->dna_hsp_param){
        seeder->dna_loader
            = Seeder_Loader_create(comparison_param->dna_hsp_param,
                seeder->sas, offsetof(Comparison, dna_hspset));
        seeder->any_hsp_param
              = HSP_Param_share(comparison_param->dna_hsp_param);
        }
    if(comparison_param->protein_hsp_param){
        seeder->protein_loader
            = Seeder_Loader_create(comparison_param->protein_hsp_param,
                seeder->sas, offsetof(Comparison, protein_hspset));
        seeder->any_hsp_param
              = HSP_Param_share(comparison_param->protein_hsp_param);
        }
    if(comparison_param->codon_hsp_param){
        seeder->codon_loader
            = Seeder_Loader_create(comparison_param->codon_hsp_param,
                seeder->sas, offsetof(Comparison, codon_hspset));
        seeder->any_hsp_param
              = HSP_Param_share(comparison_param->codon_hsp_param);
        }
    g_assert(seeder->any_hsp_param);
    if(comparison_param->protein_hsp_param
    && (comparison_param->dna_hsp_param
     || comparison_param->codon_hsp_param))
        g_error("Mixed DNA and protein seeding not implemented");
    if(Seeder_decide_fsm_type(seeder->sas->force_fsm)){ /* Use FSM */
        seeder->seeder_fsm = Seeder_FSM_create(seeder);
        seeder->seeder_vfsm = NULL;
    } else { /* Use VFSM */
        if((seeder->dna_loader && seeder->codon_loader)
        && (seeder->dna_loader->hsp_param->wordlen
         != seeder->codon_loader->hsp_param->wordlen))
                g_error("Multi-wordlength VFSMs not implemented");
        seeder->seeder_fsm = NULL;
        seeder->seeder_vfsm = Seeder_VFSM_create(seeder);
        }
    seeder->active_queryinfo_list = g_ptr_array_new();
    return seeder;
    }

void Seeder_destroy(Seeder *seeder){
    register gint i;
    register Seeder_QueryInfo *query_info;
    Comparison_Param_destroy(seeder->comparison_param);
    for(i = 0; i < seeder->query_info_list->len; i++){
        query_info = seeder->query_info_list->pdata[i];
        Seeder_QueryInfo_destroy(query_info);
        }
    g_ptr_array_free(seeder->query_info_list, TRUE);
    RecycleBin_destroy(seeder->recycle_wordinfo);
    RecycleBin_destroy(seeder->recycle_seed);
    RecycleBin_destroy(seeder->recycle_neighbour);
    RecycleBin_destroy(seeder->recycle_context);
    if(seeder->seeder_fsm)
        Seeder_FSM_destroy(seeder->seeder_fsm);
    if(seeder->seeder_vfsm)
        Seeder_VFSM_destroy(seeder->seeder_vfsm);
    if(seeder->dna_loader)
       Seeder_Loader_destroy(seeder->dna_loader);
    if(seeder->protein_loader)
       Seeder_Loader_destroy(seeder->protein_loader);
    if(seeder->codon_loader)
       Seeder_Loader_destroy(seeder->codon_loader);
    HSP_Param_destroy(seeder->any_hsp_param);
    g_ptr_array_free(seeder->active_queryinfo_list, TRUE);
    g_free(seeder);
    return;
    }

gsize Seeder_memory_usage(Seeder *seeder){
    register gsize fsm_memory;
    g_assert(seeder);
    if(seeder->seeder_fsm)
        fsm_memory = Seeder_FSM_memory_usage(seeder->seeder_fsm);
    else
        fsm_memory = Seeder_VFSM_memory_usage(seeder->seeder_vfsm);
    return sizeof(Seeder)
         + fsm_memory
         + (seeder->query_info_list->len
            *(sizeof(gpointer)+sizeof(Seeder_QueryInfo)))
         + RecycleBin_memory_usage(seeder->recycle_wordinfo)
         + RecycleBin_memory_usage(seeder->recycle_seed)
         + RecycleBin_memory_usage(seeder->recycle_neighbour)
         + RecycleBin_memory_usage(seeder->recycle_context);
    }

void Seeder_memory_info(Seeder *seeder){
    g_message("Seeder Total = %dMb", (gint)(Seeder_memory_usage(seeder)>>20));
    g_assert(seeder);
    if(seeder->seeder_fsm)
        g_message(" -> FSM memory = %dMb",
                (gint)(Seeder_FSM_memory_usage(seeder->seeder_fsm)>>20));
    else
        g_message(" -> VFSM memory = %dMb",
                (gint)(Seeder_VFSM_memory_usage(seeder->seeder_vfsm)>>20));
    g_message(" -> QueryInfo memory = %dMb",
           (gint)(seeder->query_info_list->len
            *(sizeof(gpointer)+sizeof(Seeder_QueryInfo)))>>20);
    g_message(" -> Recycle_Wordinfo Memory = %dMb",
           (gint)(RecycleBin_memory_usage(seeder->recycle_wordinfo)>>20));
    g_message(" -> Recycle_Seed Memory = %dMb",
           (gint)(RecycleBin_memory_usage(seeder->recycle_seed)>>20));
    g_message(" -> Recycle_Neighbour Memory = %dMb",
           (gint)(RecycleBin_memory_usage(seeder->recycle_neighbour)>>20));
    g_message(" -> Recycle_Context Memory = %dMb",
           (gint)(RecycleBin_memory_usage(seeder->recycle_context)>>20));
    return;
    }

static gint Seeder_get_expect(Seeder *seeder, Seeder_Loader *loader,
                              gint len){
    return (gint)((loader->saturate_expectation
                 * (len - loader->hsp_param->wordlen+1))
                 + seeder->saturate_threshold);
    }

static gboolean Seeder_WordHood_traverse(gchar *word,
                 gint score, gpointer user_data){
    register Seeder_TraverseData *traverse_data = user_data;
    register Seeder *seeder = traverse_data->seeder;
    register Seeder_WordInfo *word_info;
    g_assert(seeder);
    word_info = Seeder_add_WordInfo(traverse_data->seeder, word,
                                    traverse_data->context);
    g_assert(word_info);
    if(word_info == traverse_data->curr_word_info)
        return FALSE; /* This is not a neighbour */
    /* Make curr seed word a neighbour of this word */
    Seeder_WordInfo_add_Neighbour(seeder, word_info,
                                  traverse_data->curr_word_info);
    return FALSE;
    }

static void Seeder_insert_query(Seeder *seeder, Seeder_Context *context,
                                Sequence *query, gint frame){
    register Match_Score query_expect = 0;
    register gint i, ch, pos, orig_pos, valid_count = 0, wj_ctr = 0;
    register Seeder_WordInfo *word_info;
    register Match *match = context->loader->hsp_param->match;
    register VFSM *vfsm = seeder->seeder_vfsm?seeder->seeder_vfsm->vfsm:NULL;
    register VFSM_Int state = 0, leaf;
    register Sequence *seq_masked = Sequence_mask(query);
    register gchar *seq = Sequence_get_str(seq_masked);
    Seeder_TraverseData traverse_data;
    if(seeder->saturate_threshold){
        seeder->total_query_length += query->len;
        query_expect = Seeder_get_expect(seeder, context->loader,
                                         seeder->total_query_length);
        }
    traverse_data.seeder = seeder;
    traverse_data.context = context;
    g_assert(context);
    g_assert(query);
    g_assert(query->len >= 0);
    for(i = 0; i < query->len; i++){
        if(seeder->seeder_fsm){
            if(!Alphabet_symbol_is_member(match->comparison_alphabet,
                                         (guchar)seq[i])){
                valid_count = 0;
                continue;
                }
            valid_count++;
            if(valid_count < context->loader->hsp_param->wordlen)
                continue;
        } else {
            ch = toupper(seq[i]); /* FIXME: should filter properly */
            if(!vfsm->index[ch]){
                state = 0;
                continue;
                }
            /* FIXME: use POW2 versions where applicable */
            state = VFSM_change_state_M(vfsm, state, ch);
            if(!VFSM_state_is_leaf(vfsm, state))
                continue;
            }
        if(wj_ctr--)
            continue;
        wj_ctr = seeder->sas->word_jump - 1;
        /* FIXME: apply sat_th here */
        pos = i - context->loader->hsp_param->wordlen + 1;
        if(!Seeder_word_is_valid(match, seq_masked, pos,
                   context->loader->hsp_param->wordlen))
            continue;
        if(seeder->seeder_fsm){
            word_info = Seeder_WordInfo_create(seeder);
            word_info = Seeder_add_WordInfo(seeder, seq+pos, context);
        } else {
            leaf = VFSM_state2leaf(vfsm, state);
            word_info = seeder->seeder_vfsm->leaf[leaf];
            if(!word_info){
                word_info = Seeder_WordInfo_create(seeder);
                seeder->seeder_vfsm->leaf[leaf] = word_info;
                }
            }
        g_assert(word_info);
        if(frame)
            orig_pos = (pos * 3) + frame - 1;
        else
            orig_pos = pos;
        Seeder_WordInfo_add_Seed(seeder, word_info, context,
                                 query_expect, orig_pos);
        if(word_info->seed_list /* not blocked */
        && (!word_info->seed_list->next)){ /* 1st seed */
            traverse_data.curr_word_info = word_info;
            if(context->loader->hsp_param->wordhood)
                WordHood_traverse(context->loader->hsp_param->wordhood,
                              Seeder_WordHood_traverse, seq+pos,
                              context->loader->hsp_param->wordlen, &traverse_data);
            }
        }
    Sequence_destroy(seq_masked);
    g_free(seq);
    return;
    }
/* FIXME: optimisation : efficient VFSM wordhood traversal ? */

static Seeder_Context *Seeder_Context_create(Seeder *seeder,
                       Seeder_QueryInfo *query_info,
                       Seeder_Loader *loader){
    register Seeder_Context *context
           = RecycleBin_alloc(seeder->recycle_context);
    context->loader = loader;
    context->query_info = query_info;
    return context;
    }

static void Seeder_load_query(Seeder *seeder,
                              Seeder_QueryInfo *query_info,
                              Seeder_Loader *loader){
    register Match *match = loader->hsp_param->match;
    register gint i;
    register Sequence *aa_seq;
    register Seeder_Context *context = Seeder_Context_create(seeder,
                                                  query_info, loader);
    g_assert(match);
    g_assert(query_info->query->alphabet->type
             == match->query->alphabet->type);
    if(match->query->is_translated){
        g_assert(match->mas->translate);
        for(i = 0; i < 3; i++){
            aa_seq = Sequence_translate(query_info->query,
                                        match->mas->translate, i+1);
            Seeder_insert_query(seeder, context, aa_seq, i+1);
            Sequence_destroy(aa_seq);
            }
    } else {
        Seeder_insert_query(seeder, context, query_info->query, 0);
        }
    return;
    }

gboolean Seeder_add_query(Seeder *seeder, Sequence *query){
    register Seeder_QueryInfo *query_info;
    g_assert(seeder);
    g_assert(query);
    g_assert(!seeder->is_prepared);
    if(seeder->verbosity > 2)
        g_message("Seeder Loading query [%s]", query->id);
    query_info = Seeder_QueryInfo_create(query);
    g_ptr_array_add(seeder->query_info_list, query_info);
    if(seeder->dna_loader)
        Seeder_load_query(seeder, query_info, seeder->dna_loader);
    if(seeder->protein_loader)
        Seeder_load_query(seeder, query_info, seeder->protein_loader);
    if(seeder->codon_loader)
        Seeder_load_query(seeder, query_info, seeder->codon_loader);
    /* Seeder_memory_info(seeder); */
    return Seeder_memory_usage(seeder)
        >  (seeder->sas->fsm_memory_limit << 20);
    }

typedef struct {
       Seeder *seeder;
     Sequence *target;
         gint  curr_frame;
         gint  target_expect;
} Seeder_TargetInfo;

static void Seeder_WordInfo_seed(Seeder_QueryInfo *query_info,
                                 Seeder_TargetInfo *target_info,
                                 gint query_pos, gint target_pos,
                                 Seeder_Loader *loader){
    register HSPset *hspset;
    g_assert(query_pos >= 0);
    g_assert(target_pos >= 0);
    target_info->seeder->comparison_count++;
    if(target_info->seeder->saturate_threshold)
        target_info->target_expect = Seeder_get_expect(
                target_info->seeder, loader, target_info->target->len);
    else
        target_info->target_expect = 0;
    if(!query_info->curr_comparison){
        query_info->curr_comparison = Comparison_create(
                target_info->seeder->comparison_param,
                query_info->query, target_info->target);
        g_ptr_array_add(target_info->seeder->active_queryinfo_list,
                        query_info);
        }
    hspset = OFFSET_ITEM(HSPset*, loader->hspset_offset,
                         query_info->curr_comparison);
    HSPset_seed_hsp(hspset, query_pos, target_pos);
    return;
    }

static void Seeder_FSM_traverse_func(guint seq_pos,
                                     gpointer node_data,
                                     gpointer user_data){
    register Seeder_TargetInfo *target_info = user_data;
    register Seeder_WordInfo *word_info = node_data;
    register gint tpos, target_pos;
    register Seeder *seeder = target_info->seeder;
    register Seeder_Seed *seed;
    register Seeder_Neighbour *neighbour;
    if(target_info->curr_frame)
        tpos = (seq_pos * 3) + target_info->curr_frame - 1;
    else
        tpos = seq_pos;
    /* Ignore if over saturate threshold */
    if(seeder->saturate_threshold){
        if(!word_info->match_mailbox) /* Already blocked */
            return;
        if(word_info->match_mailbox == seeder->comparison_count){
            if(++word_info->match_count > target_info->target_expect)
                return; /* Saturated */
        } else {
            word_info->match_mailbox = seeder->comparison_count;
            word_info->match_count = 1;
            }
        }
    g_assert(word_info->seed_list || word_info->neighbour_list);
    for(seed = word_info->seed_list; seed; seed = seed->next){
        target_pos = tpos - seed->context->loader->tpos_modifier;
        g_assert(target_pos >= 0);
        Seeder_WordInfo_seed(seed->context->query_info, target_info,
                             seed->query_pos, target_pos,
                             seed->context->loader);
        }
    for(neighbour = word_info->neighbour_list; neighbour;
        neighbour = neighbour->next){
        for(seed = neighbour->word_info->seed_list; seed;
            seed = seed->next){
            target_pos = tpos - seed->context->loader->tpos_modifier;
            g_assert(target_pos >= 0);
            Seeder_WordInfo_seed(seed->context->query_info, target_info,
                                 seed->query_pos, target_pos,
                                 seed->context->loader);
            }
        }
    return;
    }

static void Seeder_VFSM_traverse_single(Seeder *seeder, gchar *seq,
                                        Seeder_TargetInfo *target_info){
    register gint i, ch;
    register VFSM_Int state = 0, leaf;
    register VFSM *vfsm = seeder->seeder_vfsm->vfsm;
    register Seeder_WordInfo *word_info;
    g_assert(vfsm);
    for(i = 0; seq[i]; i++){
        ch = toupper(seq[i]); /* FIXME: filter properly */
        if(!vfsm->index[ch]){
            state = 0;
            continue;
            }
        /* FIXME: use POW2 versions where applicable */
        state = VFSM_change_state_M(vfsm, state, ch);
        if(VFSM_state_is_leaf(vfsm, state)){
            leaf = VFSM_state2leaf(vfsm, state);
            word_info = seeder->seeder_vfsm->leaf[leaf];
            if(word_info)
                Seeder_FSM_traverse_func(i, word_info, target_info);
            }
        }
    return;
    }

static void Seeder_VFSM_traverse_ambig(Seeder *seeder, gchar *seq,
                                       Seeder_TargetInfo *target_info){
    register gint i, j, k, ch;
    register VFSM_Int state, next_state, leaf;
    register VFSM *vfsm = seeder->seeder_vfsm->vfsm;
    register Seeder_WordInfo *word_info;
    register gchar *ambig;
    register VFSM_Int *temp_state_list,
                      *curr_state_list = g_new0(VFSM_Int,
                                                seeder->sas->word_ambiguity),
                      *next_state_list = g_new0(VFSM_Int,
                                                seeder->sas->word_ambiguity);
    register gint curr_state_list_len = 1, next_state_list_len = 0;
    g_assert(vfsm);
    for(i = 0; seq[i]; i++){
        ch = toupper(seq[i]); /* FIXME: filter properly */
        if((!(ambig = Alphabet_nt2ambig(ch)))
        || ((strlen(ambig) * curr_state_list_len)
        > seeder->sas->word_ambiguity)){
            next_state_list_len = 0;
            curr_state_list_len = 1;
            curr_state_list[0] = 0;
            continue;
            }
        for(j = 0; j < curr_state_list_len; j++){
            state = curr_state_list[j];
            for(k = 0; ambig[k]; k++){
                g_assert(vfsm->index[(guchar)ambig[k]]);
                next_state = VFSM_change_state_M(vfsm,
                                  state, (guchar)ambig[k]);
                /* FIXME: use POW2 versions where applicable */
                if(next_state_list_len
                && ((next_state == next_state_list[0])
                 || (next_state == next_state_list[next_state_list_len-1])))
                    break;
                next_state_list[next_state_list_len++] = next_state;
                if(VFSM_state_is_leaf(vfsm, next_state)){
                    leaf = VFSM_state2leaf(vfsm, next_state);
                    word_info = seeder->seeder_vfsm->leaf[leaf];
                    if(word_info)
                        Seeder_FSM_traverse_func(i, word_info, target_info);
                    }
                }
            }
        curr_state_list_len = next_state_list_len;
        next_state_list_len = 0;
        Swap(next_state_list, curr_state_list, temp_state_list);
        }
    g_free(curr_state_list);
    g_free(next_state_list);
    return;
    }

static void Seeder_VFSM_traverse(Seeder *seeder, gchar *seq,
                                 Seeder_TargetInfo *target_info){
    if(seeder->sas->word_ambiguity > 1)
        Seeder_VFSM_traverse_ambig(seeder, seq, target_info);
    else
        Seeder_VFSM_traverse_single(seeder, seq, target_info);
    return;
    }

static void Seeder_prepare(Seeder *seeder){
    g_assert(!seeder->is_prepared);
    if(seeder->seeder_fsm)
        FSM_compile(seeder->seeder_fsm->fsm);
    seeder->is_prepared = TRUE;
    return;
    }

static void Seeder_FSM_traverse_ambig(Seeder *seeder, gchar *seq,
                                      FSM_Traverse_Func ftf, gpointer user_data){
    register FSM *f = seeder->seeder_fsm->fsm;
    register FSM_Node *n, *next_state;
    register gint c;
    register gchar *ambig;
    register gint i, j, k;
    register FSM_Node **temp_state_list,
                      **curr_state_list = g_new0(FSM_Node*,
                                                seeder->sas->word_ambiguity),
                      **next_state_list = g_new0(FSM_Node*,
                                                seeder->sas->word_ambiguity);
    register gint curr_state_list_len = 1, next_state_list_len = 0;
    g_assert(f->is_compiled);
    curr_state_list[0] = f->root;
    for(i = 0; seq[i]; i++){
        if((!(ambig = Alphabet_nt2ambig(seq[i])))
        || ((strlen(ambig) * curr_state_list_len)
            > seeder->sas->word_ambiguity)){
            next_state_list_len = 0;
            curr_state_list_len = 1;
            curr_state_list[0] = f->root;
            continue;
            }
        for(j = 0; j < curr_state_list_len; j++){
            n = curr_state_list[j];
            for(k = 0; ambig[k]; k++){
                c = f->traversal_filter[(guchar)ambig[k]];
                if(n[c].data)
                    ftf(i, n[c].data, user_data);
                next_state = n[c].next;
                if(next_state_list_len
                && ((next_state == next_state_list[0])
                   || (next_state ==
                       next_state_list[next_state_list_len-1]))){
                    break;
                    }
                next_state_list[next_state_list_len++] = next_state;
                }
            }
        curr_state_list_len = next_state_list_len;
        next_state_list_len = 0;
        /* g_message("set curr [%d]", curr_state_list_len); */
        Swap(next_state_list, curr_state_list, temp_state_list);
        }
    g_free(curr_state_list);
    g_free(next_state_list);
    return;
    }
/* Slow (ambiguous) version of FSM_traverse */
/* FIXME: optimisation: remove strlen() etc */

static void Seeder_FSM_traverse(Seeder *seeder, gchar *seq, FSM_Traverse_Func ftf,
                                gpointer user_data){
    if(seeder->sas->word_ambiguity > 1)
        Seeder_FSM_traverse_ambig(seeder, seq, ftf, user_data);
    else
        FSM_traverse(seeder->seeder_fsm->fsm, seq, ftf, user_data);
    return;
    }

void Seeder_add_target(Seeder *seeder, Sequence *target){
    register gint i;
    register Sequence *aa_seq;
    register gchar *seq;
    register Seeder_QueryInfo *query_info;
    register Match *match = seeder->any_hsp_param->match;
    register Sequence *target_masked;
    Seeder_TargetInfo target_info;
    g_assert(seeder);
    g_assert(target);
    g_assert(target->alphabet->type == match->target->alphabet->type);
    target_info.seeder = seeder;
    target_info.target = Sequence_share(target);
    if(!seeder->is_prepared)
        Seeder_prepare(seeder);
    if(seeder->verbosity > 2)
        g_message("Seeder finding matches with target [%s]", target->id);
    if(match->target->is_translated){
        g_assert(match->mas->translate);
        for(i = 0; i < 3; i++){
            target_info.curr_frame = i+1;
            aa_seq = Sequence_translate(target, match->mas->translate, i+1);
            target_masked = Sequence_mask(aa_seq);
            Sequence_destroy(aa_seq);
            seq = Sequence_get_str(target_masked);
            Sequence_destroy(target_masked);
            if(seeder->seeder_fsm)
                Seeder_FSM_traverse(seeder, seq,
                                    Seeder_FSM_traverse_func, &target_info);
            else
                Seeder_VFSM_traverse(seeder, seq, &target_info);
            g_free(seq);
            }
    } else {
        target_info.curr_frame = 0;
        target_masked = Sequence_mask(target);
        seq = Sequence_get_str(target_masked);
        Sequence_destroy(target_masked);
        if(seeder->seeder_fsm){
            Seeder_FSM_traverse(seeder, seq,
                                Seeder_FSM_traverse_func, &target_info);
        } else {
            Seeder_VFSM_traverse(seeder, seq, &target_info);
            }
        g_free(seq);
        }
    Sequence_destroy(target_info.target);
    if(seeder->verbosity > 2)
        g_message("Seeder done finding matches with target [%s]", target->id);
    /* Report matches */
    for(i = 0; i < seeder->active_queryinfo_list->len; i++){
        query_info = seeder->active_queryinfo_list->pdata[i];
        g_assert(query_info->curr_comparison);
        if(Comparison_has_hsps(query_info->curr_comparison)){
            Comparison_finalise(query_info->curr_comparison);
            seeder->report_func(query_info->curr_comparison,
                                seeder->user_data);
            }
        Comparison_destroy(query_info->curr_comparison);
        query_info->curr_comparison = NULL;
        }
    g_ptr_array_set_size(seeder->active_queryinfo_list, 0);
    return;
    }

/**/

/* FIXME: replace FSM/VFSM conditional stuff with funcs */

