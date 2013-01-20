/****************************************************************\
*                                                                *
*  GAM: Gapped Alignment Manager                                 *
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

#include <search.h> /* For tdelete(), tfind(), tsearch(), twalk() */
#include <stdio.h>  /* For fopen() */
#include <stdlib.h> /* For qsort() */
#include <string.h> /* For strcmp() */
#include <strings.h> /* For strcasecmp() */
#include <unistd.h> /* For unlink(), getpid(), getppid() etc */

#include "gam.h"
#include "ungapped.h"
#include "opair.h"
#include "rangetree.h"

static GAM *_gam; /* file-scope variable for passing to twalk() */

static gchar *GAM_Argument_parse_Model_Type(gchar *arg_string,
                                            gpointer data){
    gchar *type_str;
    register gchar *ret_val = Argument_parse_string(arg_string,
                                                    &type_str);
    register Model_Type *model_type = (Model_Type*)data;
    if(ret_val)
        return ret_val;
    (*model_type) = Model_Type_from_string(type_str);
    return NULL;
    }

/**/

static gchar *GAM_Refinement_to_string(GAM_Refinement refinement){
    register gchar *name = NULL;
    switch(refinement){
        case GAM_Refinement_NONE:
            name = "none";
            break;
        case GAM_Refinement_FULL:
            name = "full";
            break;
        case GAM_Refinement_REGION:
            name = "region";
            break;
        default:
            g_error("Unknown GAM Refinement [%d]", refinement);
            break;
        }
    return name;
    }

static GAM_Refinement GAM_Refinement_from_string(gchar *str){
    gchar *name[GAM_Refinement_TOTAL] =
        {"none", "full", "region"};
    GAM_Refinement refinement[GAM_Refinement_TOTAL] =
       {GAM_Refinement_NONE,
        GAM_Refinement_FULL,
        GAM_Refinement_REGION};
    register gint i;
    for(i = 0; i < GAM_Refinement_TOTAL; i++)
        if(!strcasecmp(name[i], str))
            return refinement[i];
    g_error("Unknown refinement type [%s]", str);
    return GAM_Refinement_NONE; /* Not reached */
    }

static gchar *GAM_Argument_parse_refinement(gchar *arg_string,
                                                  gpointer data){
    register GAM_Refinement *refinement
          = (GAM_Refinement*)data;
    gchar *refine_str;
    register gchar *ret_val = Argument_parse_string(arg_string,
                                                    &refine_str);
    if(ret_val)
        return ret_val;
    (*refinement) = GAM_Refinement_from_string(refine_str);
    return NULL;
    }

/**/

GAM_ArgumentSet *GAM_ArgumentSet_create(Argument *arg){
    register ArgumentSet *as;
    static GAM_ArgumentSet gas;
    if(arg){
        as = ArgumentSet_create("Gapped Alignment Options");
        ArgumentSet_add_option(as, 'm', "model", "alignment model",
        "Specify alignment model type\n"
        "Supported types:\n"
        "    ungapped ungapped:trans\n"
        "    affine:global affine:bestfit affine:local affine:overlap\n"
        "    est2genome ner protein2dna protein2genome\n"
        "    protein2dna:bestfit protein2genome:bestfit\n"
        "    coding2coding coding2genome cdna2genome genome2genome\n",
         "ungapped",
         GAM_Argument_parse_Model_Type, &gas.type);
        ArgumentSet_add_option(as, 's', "score", "threshold",
         "Score threshold for gapped alignment", "100",
         Argument_parse_int, &gas.threshold);
        ArgumentSet_add_option(as, '\0', "percent", "threshold",
         "Percent self-score threshold", "0.0",
         Argument_parse_float, &gas.percent_threshold);
        ArgumentSet_add_option(as, 0, "showalignment", NULL,
         "Include (human readable) alignment in results", "TRUE",
         Argument_parse_boolean, &gas.show_alignment);
        ArgumentSet_add_option(as, 0, "showsugar", NULL,
         "Include 'sugar' format output in results", "FALSE",
         Argument_parse_boolean, &gas.show_sugar);
        ArgumentSet_add_option(as, 0, "showcigar", NULL,
         "Include 'cigar' format output in results", "FALSE",
         Argument_parse_boolean, &gas.show_cigar);
        ArgumentSet_add_option(as, 0, "showvulgar", NULL,
         "Include 'vulgar' format output in results", "TRUE",
         Argument_parse_boolean, &gas.show_vulgar);
        ArgumentSet_add_option(as, 0, "showquerygff", NULL,
         "Include GFF output on query in results", "FALSE",
         Argument_parse_boolean, &gas.show_query_gff);
        ArgumentSet_add_option(as, 0, "showtargetgff", NULL,
         "Include GFF output on target in results", "FALSE",
         Argument_parse_boolean, &gas.show_target_gff);
        ArgumentSet_add_option(as, 0, "ryo", "format",
         "Roll-your-own printf-esque output format", "NULL",
         Argument_parse_string, &gas.ryo);
        ArgumentSet_add_option(as, 'n', "bestn", "number",
         "Report best N results per query", "0",
         Argument_parse_int, &gas.best_n);
        ArgumentSet_add_option(as, 'S', "subopt", NULL,
         "Search for suboptimal alignments", "TRUE",
         Argument_parse_boolean, &gas.use_subopt);
        ArgumentSet_add_option(as, 'g', "gappedextension", NULL,
         "Use gapped extension (default is SDP)", "TRUE",
         Argument_parse_boolean, &gas.use_gapped_extension);
        /**/
        ArgumentSet_add_option(as, '\0', "refine", NULL,
        "Alignment refinement strategy [none|full|region]", "none",
        GAM_Argument_parse_refinement, &gas.refinement);
        ArgumentSet_add_option(as, '\0', "refineboundary", NULL,
        "Refinement region boundary", "32",
        Argument_parse_int, &gas.refinement_boundary);
        /**/
        Argument_absorb_ArgumentSet(arg, as);
        }
    return &gas;
    }

/**/

static gboolean GAM_StoredResult_compare(gpointer low,
                                         gpointer high,
                                         gpointer user_data){
    register GAM_StoredResult *gsr_low = (GAM_StoredResult*)low,
                              *gsr_high = (GAM_StoredResult*)high;
    return gsr_low->score < gsr_high->score;
    }

static void GAM_display_alignment(GAM *gam,
            Alignment *alignment, Sequence *query, Sequence *target,
            gint result_id, gint rank,
            gpointer user_data, gpointer self_data, FILE *fp);

static GAM_StoredResult *GAM_StoredResult_create(GAM *gam,
                         Sequence *query, Sequence *target,
                         Alignment *alignment,
                         gpointer user_data, gpointer self_data){
    register GAM_StoredResult *gsr = g_new(GAM_StoredResult, 1);
    gsr->score = alignment->score;
    gsr->pos = ftell(gam->bestn_tmp_file);
    GAM_display_alignment(gam, alignment, query, target,
                          0, -1, user_data, self_data, gam->bestn_tmp_file);
    gsr->len = ftell(gam->bestn_tmp_file)
             - gsr->pos;
    return gsr;
    }

static void GAM_StoredResult_destroy(GAM_StoredResult *gsr){
    g_free(gsr);
    return;
    }
/* FIXME: optimisation:
 *        could store unused positions in the tmp file for reuse
 */

static void GAM_StoredResult_display(GAM_StoredResult *gsr,
                                     GAM *gam, gint rank){
    register gint i, ch, tag_pos = 0;
    register gchar *rank_tag = "%_EXONERATE_BESTN_RANK_%";
    if(fseek(gam->bestn_tmp_file, gsr->pos, SEEK_SET))
        g_error("Could not seek in tmp file");
    for(i = 0; i < gsr->len; i++){
        ch = getc(gam->bestn_tmp_file);
        g_assert(ch != EOF);
        if(ch == rank_tag[tag_pos]){
            tag_pos++;
            if(!rank_tag[tag_pos]){
                printf("%d", rank);
                tag_pos = 0;
                }
        } else {
            if(tag_pos){
                printf("%.*s", tag_pos, rank_tag);
                tag_pos = 0;
                }
            putchar(ch);
            }
        }
    fflush(stdout);
    return;
    }

/**/

static int GAM_compare_id(const void *gqr1, const void *gqr2){
    return strcmp(((const GAM_QueryResult*)gqr1)->query_id, 
                  ((const GAM_QueryResult*)gqr2)->query_id);
    }

static GAM_QueryResult *GAM_QueryResult_create(GAM *gam,
                                               gchar *query_id){
    register GAM_QueryResult *gqr = g_new(GAM_QueryResult, 1);
    gqr->pq = PQueue_create(gam->pqueue_set,
                            GAM_StoredResult_compare, NULL);
    gqr->query_id = g_strdup(query_id);
    gqr->tie_count = 0;
    gqr->tie_score = C4_IMPOSSIBLY_LOW_SCORE;
    return gqr;
    }

static void GAM_QueryResult_pqueue_destroy_GAM_Result(gpointer data,
                                                 gpointer user_data){
    register GAM_StoredResult *gsr = data;
    GAM_StoredResult_destroy(gsr);
    return;
    }

static void GAM_QueryResult_destroy(GAM_QueryResult *gqr){
    PQueue_destroy(gqr->pq,
                   GAM_QueryResult_pqueue_destroy_GAM_Result, NULL);
    g_free(gqr->query_id);
    g_free(gqr);
    return;
    }

static void GAM_QueryResult_push(GAM_QueryResult *gqr,
                                 GAM_Result *gam_result,
                                 Alignment *alignment){
    register GAM_StoredResult *gsr = GAM_StoredResult_create(gam_result->gam,
                                          gam_result->query,
                                          gam_result->target,
                                          alignment,
                                          gam_result->user_data,
                                          gam_result->self_data);
    PQueue_push(gqr->pq, gsr);
    return;
    }

static void GAM_QueryResult_submit(GAM_QueryResult *gqr,
                                   GAM_Result *gam_result){
    register GAM_StoredResult *gsr;
    register Alignment *alignment;
    register gint i;
    register GPtrArray *tie_list = g_ptr_array_new();
    g_assert(!strcmp(gqr->query_id, gam_result->query->id));
    for(i = 0; i < gam_result->alignment_list->len; i++){
        alignment = gam_result->alignment_list->pdata[i];
        if(alignment->score == gqr->tie_score){
            GAM_QueryResult_push(gqr, gam_result, alignment);
            gqr->tie_count++;
        } else if(alignment->score < gqr->tie_score){
            if(PQueue_total(gqr->pq) < gam_result->gam->gas->best_n){
                GAM_QueryResult_push(gqr, gam_result, alignment);
                gqr->tie_count = 1;
                gqr->tie_score = alignment->score;
            } else {
                break; /* Other alignments are worse */
                }
        } else { /* (alignment->score > gqr->tie_score) */
            GAM_QueryResult_push(gqr, gam_result, alignment);
            if((PQueue_total(gqr->pq)-gqr->tie_count)
               >= gam_result->gam->gas->best_n){
                /* Remove old ties */
                while(gqr->tie_count){
                    gsr = PQueue_pop(gqr->pq);
                    GAM_StoredResult_destroy(gsr);
                    gqr->tie_count--;
                    }
                /* Count new ties */
                gsr = PQueue_top(gqr->pq);
                gqr->tie_score = gsr->score;
                do {
                    gsr = PQueue_top(gqr->pq);
                    if(gsr && (gsr->score == gqr->tie_score)){
                        gsr = PQueue_pop(gqr->pq);
                        g_ptr_array_add(tie_list, gsr);
                    } else {
                        break;
                        }
                } while(TRUE);
                gqr->tie_count = tie_list->len;
                /* Replace new ties */
                while(tie_list->len){
                    gsr = tie_list->pdata[tie_list->len-1];
                    PQueue_push(gqr->pq, gsr);
                    g_ptr_array_set_size(tie_list, tie_list->len-1);
                    }
            } else {
                if(PQueue_total(gqr->pq) == 1){ /* First alignment */
                    gqr->tie_count = 1;
                    gqr->tie_score = alignment->score;
                    }
                }
            }
        }
    g_ptr_array_free(tie_list, TRUE);
    return;
    }

/**/

static gboolean GAM_QueryResult_report_traverse_func(gpointer data,
                                                     gpointer user_data){
    register GAM_StoredResult *gsr = data;
    register GPtrArray *result_list = user_data;
    g_ptr_array_add(result_list, gsr);
    return FALSE;
    }

static int GAM_QueryResult_report_sort_func(const void *a,
                                            const void *b){
    register GAM_StoredResult **gsr_a = (GAM_StoredResult**)a,
                              **gsr_b = (GAM_StoredResult**)b;
    return (*gsr_a)->score - (*gsr_b)->score;
    }

static void GAM_QueryResult_report(GAM_QueryResult *gqr, GAM *gam){
    register GPtrArray *result_list = g_ptr_array_new();
    register gint i;
    register GAM_StoredResult *gsr;
    if(gam->verbosity > 2)
        g_message("Reporting [%d/%d] results for query [%s]",
              PQueue_total(gqr->pq), gam->gas->best_n,
              gqr->query_id);
    PQueue_traverse(gqr->pq, GAM_QueryResult_report_traverse_func,
                    result_list);
    qsort(result_list->pdata, result_list->len,
          sizeof(gpointer), GAM_QueryResult_report_sort_func);
    for(i = result_list->len-1; i >= 0; i--){
        gsr = result_list->pdata[i];
        GAM_StoredResult_display(gsr, gam, result_list->len-i);
        }
    g_ptr_array_free(result_list, TRUE);
    return;
    }

/**/

#define Combined_Alphabet_Type(qt,tt) ((qt) << 8 | (tt))

static GPtrArray *GAM_build_match_list(C4_Model *model){
    register GPtrArray *match_list = g_ptr_array_new();
    Match *match_array[Match_Type_TOTAL] = {0};
    register GPtrArray *match_transition_list
        = C4_Model_select_transitions(model, C4_Label_MATCH);
    register gint i;
    register C4_Transition *transition;
    register Match *match;
    for(i = 0; i < match_transition_list->len; i++){
        transition = match_transition_list->pdata[i];
        g_assert(transition->label == C4_Label_MATCH);
        g_assert(transition->label_data);
        match = transition->label_data;
        if(!match_array[match->type]){
            match_array[match->type] = match;
            g_ptr_array_add(match_list, match);
            }
        }
    g_assert(match_list->len);
    g_ptr_array_free(match_transition_list, TRUE);
    return match_list;
    }

GAM *GAM_create(Alphabet_Type query_type, Alphabet_Type target_type,
                Submat *dna_submat, Submat *protein_submat,
                Translate *translate, gboolean use_exhaustive,
                gint verbosity){
    register GAM *gam = g_new0(GAM, 1);
    register gint i;
    register C4_Span *span;
    gam->thread_ref = ThreadRef_create();
    gam->dna_submat = Submat_share(dna_submat);
    gam->protein_submat = Submat_share(protein_submat);
    gam->translate = Translate_share(translate);
    gam->gas = GAM_ArgumentSet_create(NULL);
    if(use_exhaustive && gam->gas->use_subopt)
        g_warning("Exhaustively generating suboptimal alignments"
                  " will be VERY SLOW: use -S no");
    gam->query_type = query_type;
    gam->target_type = target_type;
    if(gam->gas->best_n){
        gam->bestn_tmp_file = tmpfile();
        gam->pqueue_set = PQueueSet_create();
        }
    gam->translate_both = Model_Type_translate_both(gam->gas->type);
    gam->dual_match = Model_Type_has_dual_match(gam->gas->type);
    gam->model = Model_Type_get_model(gam->gas->type,
                                      query_type, target_type);
    if((!use_exhaustive) && (!C4_Model_is_local(gam->model)))
        g_error("Cannot perform heuristic alignments using non-local models: use -E");
    gam->match_list = GAM_build_match_list(gam->model);
    if(use_exhaustive){
        gam->optimal = Optimal_create(gam->model, NULL,
                       Optimal_Type_SCORE
                      |Optimal_Type_PATH
                      |Optimal_Type_REDUCED_SPACE, TRUE);
        if(gam->gas->refinement != GAM_Refinement_NONE)
            g_error("Exhaustive alignments cannot be refined");
    } else {
        if(gam->gas->refinement != GAM_Refinement_NONE){
            gam->optimal = Optimal_create(gam->model, NULL,
                           Optimal_Type_SCORE
                          |Optimal_Type_PATH
                          |Optimal_Type_REDUCED_SPACE, TRUE);
            }
        if(Model_Type_is_gapped(gam->gas->type)){
            if(gam->gas->use_gapped_extension){
                gam->sdp = SDP_create(gam->model);
            } else {
                gam->heuristic = Heuristic_create(gam->model);
                }
            }
        }
    gam->verbosity = verbosity;
    /* Find max_{query,target}_span */
    gam->max_query_span = gam->max_target_span = 0;
    for(i = 0; i < gam->model->span_list->len; i++){
        span = gam->model->span_list->pdata[i];
        if(gam->max_query_span < span->max_query)
            gam->max_query_span = span->max_query;
        if(gam->max_target_span < span->max_target)
            gam->max_target_span = span->max_target;
        }
#ifdef USE_PTHREADS
    pthread_mutex_init(&gam->gam_lock, NULL);
#endif /* USE_PTHREADS */
    return gam;
    }

GAM *GAM_share(GAM *gam){
    g_assert(gam);
    ThreadRef_share(gam->thread_ref);
    return gam;
    }

/**/

static GAM_QueryInfo *GAM_QueryInfo_create(Sequence *query, GAM *gam){
    register GAM_QueryInfo *gqi = g_new(GAM_QueryInfo, 1);
    register gint i, j;
    register Match *match;
    register Match_Score th;
    gqi->query_id = g_strdup(query->id);
    gqi->threshold = 0;
    /* Calculate best threshold for each query match */
    for(i = 0; i < gam->match_list->len; i++){
        match = gam->match_list->pdata[i];
        th = 0;
        for(j = 0; j < query->len; j += match->query->advance)
            th += match->query->self_func(match->query, query, j);
        if(gqi->threshold < th)
            gqi->threshold = th;
        }
    gqi->threshold *= gam->gas->percent_threshold;
    gqi->threshold /= 100;
    if(gqi->threshold < gam->gas->threshold)
        gqi->threshold = gam->gas->threshold;
    return gqi;
    }

static void GAM_QueryInfo_destroy(GAM_QueryInfo *gqi){
    g_free(gqi->query_id);
    g_free(gqi);
    return;
    }

/**/

static void GAM_bestn_tree_destroy(void *bestn_tree){
    while (bestn_tree) {
        GAM_QueryResult *gqr = *(GAM_QueryResult **)bestn_tree;
        tdelete((void *)gqr, &bestn_tree, GAM_compare_id);
        GAM_QueryResult_destroy(gqr);
        }
    }

static void GAM_percent_threshold_tree_destroy(void *percent_threshold_tree){
    while (percent_threshold_tree) {
        GAM_QueryInfo *gqi = *(GAM_QueryInfo **)percent_threshold_tree;
        tdelete((void *)gqi, &percent_threshold_tree, GAM_compare_id);
        GAM_QueryInfo_destroy(gqi);
        }
    }

void GAM_destroy(GAM *gam){
    g_assert(gam);
    if(ThreadRef_destroy(gam->thread_ref))
        return;
#ifdef USE_PTHREADS
    pthread_mutex_destroy(&gam->gam_lock);
#endif /* USE_PTHREADS */
    g_assert(gam->model);
    g_ptr_array_free(gam->match_list, TRUE);
    C4_Model_destroy(gam->model);
    if(gam->optimal)
        Optimal_destroy(gam->optimal);
    if(gam->heuristic)
        Heuristic_destroy(gam->heuristic);
    if(gam->sdp)
        SDP_destroy(gam->sdp);
    if(gam->translate)
        Translate_destroy(gam->translate);
    if(gam->bestn_tree)
        GAM_bestn_tree_destroy(gam->bestn_tree);
    if(gam->percent_threshold_tree)
        GAM_percent_threshold_tree_destroy(gam->percent_threshold_tree);
    if(gam->bestn_tmp_file)
        fclose(gam->bestn_tmp_file);
    if(gam->pqueue_set)
        PQueueSet_destroy(gam->pqueue_set);
    g_free(gam);
    return;
    }

void GAM_bestn_tree_report_traverse(const void *gqr,
                                    VISIT order,
                                    int level) {
    if (order == leaf || order == postorder)
        GAM_QueryResult_report(*(GAM_QueryResult**)gqr, _gam);
    }

void GAM_report(GAM *gam){
    if(gam->bestn_tree) {
        _gam = gam;
        twalk(gam->bestn_tree, GAM_bestn_tree_report_traverse);
    }
    return;
    }

/**/

static C4_Portal *GAM_Pair_find_portal(C4_Model *model,
                                       HSPset *hspset){
    register gint i;
    register C4_Portal *portal = NULL;
    register C4_Transition *transition;
    for(i = 0; i < model->portal_list->len; i++){
        g_assert(model->portal_list);
        portal = model->portal_list->pdata[i];
        g_assert(portal->transition_list->len);
        transition = portal->transition_list->pdata[0];
        g_assert(transition);
        if((transition->advance_query
            == hspset->param->match->query->advance)
        && (transition->advance_target
            == hspset->param->match->target->advance)){
            return portal;
            }
        }
    if(!portal)
        g_error("No compatible portal found for hspset");
    return portal;
    }

static GAM_Result *GAM_Result_create(GAM *gam,
                                     Sequence *query,
                                     Sequence *target){
    register GAM_Result *gam_result = g_new(GAM_Result, 1);
    g_assert(gam);
    g_assert(query);
    g_assert(target);
    g_assert(query->alphabet->type == gam->query_type);
    g_assert(target->alphabet->type == gam->target_type);
    gam_result->ref_count = 1;
    gam_result->gam = GAM_share(gam);
    gam_result->alignment_list = NULL;
    gam_result->query = Sequence_share(query);
    gam_result->target = Sequence_share(target);
    gam_result->user_data = Model_Type_create_data(gam->gas->type,
                                                   query, target);
    gam_result->self_data = Model_Type_create_data(gam->gas->type,
                                                   query, query);
    gam_result->subopt = SubOpt_create(query->len, target->len);
    return gam_result;
    }

static Alignment *GAM_Result_refine_alignment(GAM_Result *gam_result,
                                              Alignment *alignment){
    register Alignment *refined_alignment = NULL;
    register Region *region;
    register gint query_region_start, target_region_start;
    g_assert(gam_result->gam->optimal);
    if(gam_result->gam->verbosity > 1)
        g_message("Refining alignment ... (%d)", alignment->score);
    switch(gam_result->gam->gas->refinement){
        case GAM_Refinement_FULL:
            region = Region_create(0, 0, gam_result->query->len,
                                         gam_result->target->len);
            refined_alignment = Optimal_find_path(
                    gam_result->gam->optimal, region,
                    gam_result->user_data, 0, gam_result->subopt);
            g_assert(refined_alignment);
            Region_destroy(region);
            break;
        case GAM_Refinement_REGION:
            query_region_start
                = MAX(0, alignment->region->query_start
                       - gam_result->gam->gas->refinement_boundary);
            target_region_start
                = MAX(0, alignment->region->target_start
                       - gam_result->gam->gas->refinement_boundary);
            region = Region_create(query_region_start,
                                   target_region_start,
                MIN(gam_result->query->len,
                    Region_query_end(alignment->region)
                  + gam_result->gam->gas->refinement_boundary)
                - query_region_start,
                MIN(gam_result->target->len,
                    Region_target_end(alignment->region)
                  + gam_result->gam->gas->refinement_boundary)
                - target_region_start);
            refined_alignment = Optimal_find_path(
                    gam_result->gam->optimal, region,
                    gam_result->user_data, 0, gam_result->subopt);
            g_assert(refined_alignment);
            Region_destroy(region);
            break;
        default:
            g_error("Bad type for refinement [%s]",
                GAM_Refinement_to_string(gam_result->gam->gas->refinement));
            break;
        }
    g_assert(refined_alignment);
    if(gam_result->gam->verbosity > 1)
        g_message("Refined alignment score [%d]", refined_alignment->score);
    return refined_alignment;
    }

static void GAM_Result_add_alignment(GAM_Result *gam_result,
                                     Alignment *alignment,
                                     C4_Score threshold){
    register Alignment *refined_alignment;
    if(!gam_result->alignment_list)
        gam_result->alignment_list = g_ptr_array_new();
    if(gam_result->gam->gas->refinement != GAM_Refinement_NONE){
        refined_alignment = GAM_Result_refine_alignment(gam_result, alignment);
        /* Use refined alignment only if it has a higher score */
        if (refined_alignment->score >= alignment->score) {
            Alignment_destroy(alignment);
            alignment = refined_alignment;
        } else
            Alignment_destroy(refined_alignment);
    }
    g_ptr_array_add(gam_result->alignment_list, alignment);
    SubOpt_add_alignment(gam_result->subopt, alignment);
    return;
    }

static C4_Score GAM_get_query_threshold(GAM *gam, Sequence *query){
    register GAM_QueryResult *gqr;
    register GAM_StoredResult *gsr;
    register GAM_QueryInfo *gqi;
    GAM_QueryResult gq_lookup = {.query_id = query->id};
    void *tree_node;
    if(gam->gas->best_n){
        tree_node = tfind((void*)&gq_lookup, &gam->bestn_tree, GAM_compare_id);
        gqr = tree_node ? *(GAM_QueryResult **)tree_node : NULL;
        if(gqr && (PQueue_total(gqr->pq) >= gam->gas->best_n)){
            gsr = PQueue_top(gqr->pq);
            if(gam->verbosity > 2)
                g_message("Using threshold [%d] for query [%s]",
                      gsr->score, query->id);
            return gsr->score;
            }
        }
    if(gam->gas->percent_threshold){
        tree_node = tfind((void*)&gq_lookup, &gam->percent_threshold_tree, 
                          GAM_compare_id);
        gqi = tree_node ? *(GAM_QueryInfo **)tree_node : NULL;
        if(!gqi){
            gqi = GAM_QueryInfo_create(query, gam);
            tsearch((void*)gqi, &gam->percent_threshold_tree, GAM_compare_id);
            }
        return gqi->threshold;
        }
    return gam->gas->threshold;
    }

static int GAM_Result_ungapped_create_sort_func(const void *a,
                                                const void *b){
    register Alignment **alignment_a = (Alignment**)a,
                       **alignment_b = (Alignment**)b;
    return (*alignment_b)->score - (*alignment_a)->score;
    }

static void GAM_Result_ungapped_add_HSPset(GAM_Result *gam_result,
                                           C4_Model *model,
                                           HSPset *hspset){
    register gint i;
    register C4_Score threshold;
    register Alignment *alignment;
    register HSP *hsp;
    HSPset_filter_ungapped(hspset);
    for(i = 0; i < hspset->hsp_list->len; i++){
        hsp = hspset->hsp_list->pdata[i];
        threshold = GAM_get_query_threshold(gam_result->gam, hspset->query);
        if(hsp->score >= threshold){
            alignment = Ungapped_Alignment_create(model,
                                                  gam_result->user_data,
                                                  hsp);
            g_assert(alignment);
            GAM_Result_add_alignment(gam_result, alignment, threshold);
            }
        }
    return;
    }

GAM_Result *GAM_Result_ungapped_create(GAM *gam,
                                       Comparison *comparison){
    register GAM_Result *gam_result;
    g_assert(comparison);
    if(!Comparison_has_hsps(comparison))
        return NULL;
    gam_result = GAM_Result_create(gam, comparison->query,
                                        comparison->target);
    /**/
    if(comparison->dna_hspset)
        GAM_Result_ungapped_add_HSPset(gam_result, gam->model,
                                       comparison->dna_hspset);
    if(comparison->protein_hspset)
        GAM_Result_ungapped_add_HSPset(gam_result, gam->model,
                                       comparison->protein_hspset);
    if(comparison->codon_hspset)
        GAM_Result_ungapped_add_HSPset(gam_result, gam->model,
                                       comparison->codon_hspset);
    /**/
    if(!gam_result->alignment_list){
        GAM_Result_destroy(gam_result);
        return NULL;
        }
    /* Sort results, best scores first */
    qsort(gam_result->alignment_list->pdata,
          gam_result->alignment_list->len,
          sizeof(gpointer), GAM_Result_ungapped_create_sort_func);
    return gam_result;
    }

/**/

static void GAM_Result_BSDP_add_HSPset(HSPset *hspset, GAM *gam,
                                       HPair *hpair){
    register C4_Portal *portal;
    portal = GAM_Pair_find_portal(gam->model, hspset);
    HPair_add_hspset(hpair, portal, hspset);
    if(gam->verbosity > 2)
        g_message("Added [%d] hsps to HPair",
                  hspset->hsp_list->len);
    return;
    }

static gboolean GAM_Result_is_full(GAM_Result *gam_result){
    register Alignment *penult, *last;
    if((gam_result->gam->gas->best_n)
    && (gam_result->alignment_list->len >= gam_result->gam->gas->best_n)
    && (gam_result->alignment_list->len > 1)){
        penult = gam_result->alignment_list->pdata
                [gam_result->alignment_list->len-2];
        last   = gam_result->alignment_list->pdata
                [gam_result->alignment_list->len-1];
        if(penult->score != last->score)
            return TRUE;
        }
    return FALSE;
    }
/* GAM_Result is only full when best_n threshold has been met
 * and all tie-breakers have been admitted into the GAM_Result
 */

static GAM_Result *GAM_Result_BSDP_create(GAM *gam,
                                          Comparison *comparison){
    register HPair *hpair;
    register Alignment *alignment;
    register C4_Score threshold
             = GAM_get_query_threshold(gam, comparison->query);
    register GAM_Result *gam_result = GAM_Result_create(gam,
                                        comparison->query,
                                        comparison->target);
    g_assert(gam->heuristic);
    if(gam->verbosity > 2)
        g_message("Preparing HPair for [%s][%s]",
                  comparison->query->id, comparison->target->id);
    hpair = HPair_create(gam->heuristic, gam_result->subopt,
                         comparison->query->len,
                         comparison->target->len,
                         gam->verbosity, gam_result->user_data);
    /**/
    g_assert(comparison);
    g_assert(Comparison_has_hsps(comparison));
    if(comparison->dna_hspset)
        GAM_Result_BSDP_add_HSPset(comparison->dna_hspset, gam, hpair);
    if(comparison->protein_hspset)
        GAM_Result_BSDP_add_HSPset(comparison->protein_hspset, gam,
                                   hpair);
    if(comparison->codon_hspset)
        GAM_Result_BSDP_add_HSPset(comparison->codon_hspset, gam,
                                   hpair);
    HPair_finalise(hpair, threshold);
    if(gam->verbosity > 2)
        g_message("Finalised HPair");
    do {
        threshold = GAM_get_query_threshold(gam, comparison->query);
        alignment = HPair_next_path(hpair, threshold);
        if(!alignment)
            break;
        g_assert(alignment->score >= threshold);
        if(gam->verbosity > 2)
            g_message("Found BSDP alignment score [%d]",
                      alignment->score);
        GAM_Result_add_alignment(gam_result, alignment, threshold);
        if(GAM_Result_is_full(gam_result))
            break;
    } while(gam->gas->use_subopt);
    HPair_destroy(hpair);
    if(!gam_result->alignment_list){ /* No alignments */
        GAM_Result_destroy(gam_result);
        return NULL;
        }
    return gam_result;
    }
/* FIXME: optimisation: change to update threshold between
 *        suboptimal alignments.
 */

static GAM_Result *GAM_Result_SDP_create(GAM *gam,
                                         Comparison *comparison){
    register SDP_Pair *sdp_pair;
    register Alignment *alignment;
    register C4_Score threshold;
    register GAM_Result *gam_result = GAM_Result_create(gam,
                                        comparison->query,
                                        comparison->target);
    if(gam->verbosity > 2)
        g_message("Preparing SDP_Pair for [%s][%s]",
                comparison->query->id, comparison->target->id);
    g_assert(gam);
    g_assert(gam->sdp);
    g_assert(comparison);
    g_assert(Comparison_has_hsps(comparison));
    sdp_pair = SDP_Pair_create(gam->sdp, gam_result->subopt,
                               comparison, gam_result->user_data);
    do {
        threshold = GAM_get_query_threshold(gam, comparison->query);
        alignment = SDP_Pair_next_path(sdp_pair, threshold);
        if(!alignment)
            break;
        g_assert(alignment->score >= threshold);
        if(gam->verbosity > 2)
            g_message("Found SDP alignment score [%d]",
                      alignment->score);
        GAM_Result_add_alignment(gam_result, alignment, threshold);
        if(GAM_Result_is_full(gam_result))
            break;
    } while(gam->gas->use_subopt);
    SDP_Pair_destroy(sdp_pair);
    if(!gam_result->alignment_list){ /* No alignments */
        GAM_Result_destroy(gam_result);
        return NULL;
        }
    return gam_result;
    }

/**/

typedef struct {
     gboolean *fwd_keep;
     gboolean *rev_keep;
    GPtrArray *hsp_list;
    RangeTree *rangetree;
          HSP *max_cobs_hsp;
          GAM *gam;
} GAM_Geneseed_Data;

static void GAM_Result_geneseed_add(HSPset *hspset, GAM_Geneseed_Data *gsd){
    register gint i;
    register HSP *hsp;
    if(!hspset)
        return;
    for(i = 0; i < hspset->hsp_list->len; i++){
        hsp = hspset->hsp_list->pdata[i];
        if(!RangeTree_check_pos(gsd->rangetree,
                                HSP_query_cobs(hsp), HSP_target_cobs(hsp)))
            RangeTree_add(gsd->rangetree,
                          HSP_query_cobs(hsp), HSP_target_cobs(hsp),
                          GINT_TO_POINTER(gsd->hsp_list->len));
        g_ptr_array_add(gsd->hsp_list, hsp);
        if((!gsd->max_cobs_hsp) || (gsd->max_cobs_hsp->cobs < hsp->cobs))
                gsd->max_cobs_hsp = hsp;
        }
    return;
    }

static void GAM_Result_geneseed_visit_hsp_fwd(GAM_Geneseed_Data *gsd,
                                              gint hsp_id);
static void GAM_Result_geneseed_visit_hsp_rev(GAM_Geneseed_Data *gsd,
                                              gint hsp_id);

static gboolean GAM_Result_geneseed_report_fwd_func(gint x, gint y,
                          gpointer info, gpointer user_data){
    register GAM_Geneseed_Data *gsd = user_data;
    register gint hsp_id = GPOINTER_TO_INT(info);
    GAM_Result_geneseed_visit_hsp_fwd(gsd, hsp_id);
    return FALSE;
    }

static gboolean GAM_Result_geneseed_report_rev_func(gint x, gint y,
                          gpointer info, gpointer user_data){
    register GAM_Geneseed_Data *gsd = user_data;
    register gint hsp_id = GPOINTER_TO_INT(info);
    GAM_Result_geneseed_visit_hsp_rev(gsd, hsp_id);
    return FALSE;
    }

static void GAM_Result_geneseed_visit_hsp_fwd(GAM_Geneseed_Data *gsd,
                                              gint hsp_id){
    register HSP *hsp;
    register gint query_range, target_range;
    g_assert(hsp_id < gsd->hsp_list->len);
    if(!gsd->fwd_keep[hsp_id]){
        gsd->fwd_keep[hsp_id] = TRUE;
        hsp = gsd->hsp_list->pdata[hsp_id];
        query_range = gsd->gam->max_query_span
                    + (((HSP_query_end(hsp) - HSP_query_cobs(hsp))
                       + (HSP_query_cobs(gsd->max_cobs_hsp)
                         - gsd->max_cobs_hsp->query_start)) * 2);
        target_range = gsd->gam->max_target_span
                     + (((HSP_target_end(hsp) - HSP_target_cobs(hsp))
                       + (HSP_target_cobs(gsd->max_cobs_hsp)
                         - gsd->max_cobs_hsp->target_start)) * 2);
        /* Search forwards */
        /*
        g_message("find fwd (%d,%d,%d) (%d,%d,%d) [%d,%d]",
                    (HSP_query_end(hsp) - HSP_query_cobs(hsp)),
                    gsd->gam->max_query_span,
                    (HSP_query_cobs(gsd->max_cobs_hsp)
                      - gsd->max_cobs_hsp->query_start),
                    (HSP_target_end(hsp) - HSP_target_cobs(hsp)),
                    gsd->gam->max_target_span,
                    (HSP_target_cobs(gsd->max_cobs_hsp)
                      - gsd->max_cobs_hsp->target_start),
                    query_range, target_range);
        g_message("RangeTree_find fwd");
        */
        RangeTree_find(gsd->rangetree,
                       HSP_query_cobs(hsp), query_range,
                       HSP_target_cobs(hsp), target_range,
                       GAM_Result_geneseed_report_fwd_func, gsd);
        }
    return;
    }

static void GAM_Result_geneseed_visit_hsp_rev(GAM_Geneseed_Data *gsd,
                                              gint hsp_id){
    register HSP *hsp;
    register gint query_range, target_range;
    g_assert(hsp_id < gsd->hsp_list->len);
    if(!gsd->rev_keep[hsp_id]){
        gsd->rev_keep[hsp_id] = TRUE;
        hsp = gsd->hsp_list->pdata[hsp_id];
        query_range = gsd->gam->max_query_span
                    + (((HSP_query_end(hsp) - HSP_query_cobs(hsp))
                       + (HSP_query_cobs(gsd->max_cobs_hsp)
                         - gsd->max_cobs_hsp->query_start)) * 2);
        target_range = gsd->gam->max_target_span
                     + (((HSP_target_end(hsp) - HSP_target_cobs(hsp))
                       + (HSP_target_cobs(gsd->max_cobs_hsp)
                         - gsd->max_cobs_hsp->target_start)) * 2);
        /*
        g_message("find rev (%d,%d,%d) (%d,%d,%d) [%d,%d]",
                    (HSP_query_end(hsp) - HSP_query_cobs(hsp)),
                    gsd->gam->max_query_span,
                    (HSP_query_cobs(gsd->max_cobs_hsp)
                      - gsd->max_cobs_hsp->query_start),
                    (HSP_target_end(hsp) - HSP_target_cobs(hsp)),
                    gsd->gam->max_target_span,
                    (HSP_target_cobs(gsd->max_cobs_hsp)
                      - gsd->max_cobs_hsp->target_start),
                    query_range, target_range);
        g_message("RangeTree_find rev");
        */
        /* Search backwards */
        RangeTree_find(gsd->rangetree,
                       HSP_query_cobs(hsp)-query_range, query_range,
                       HSP_target_cobs(hsp)-target_range, target_range,
                       GAM_Result_geneseed_report_rev_func, gsd);
        }
    return;
    }

static gint GAM_Result_geneseed_select(HSPset *hspset,
                                       gboolean *fwd_keep,
                                       gboolean *rev_keep,
                                       gint hsp_id){
    register gint i;
    register GPtrArray *keep_list = g_ptr_array_new();
    register HSP *hsp;
    if(!hspset)
        return hsp_id;
    for(i = 0; i < hspset->hsp_list->len; i++){
        hsp = hspset->hsp_list->pdata[i];
        if(fwd_keep[hsp_id] || rev_keep[hsp_id]){
            g_ptr_array_add(keep_list, hsp);
        } else {
            HSP_destroy(hsp);
            }
        hsp_id++;
        }
    /*
    g_message("geneseed selected [%d] from [%d]", keep_list->len,
                                                  hspset->hsp_list->len);
    */
    g_ptr_array_free(hspset->hsp_list, TRUE);
    hspset->hsp_list = keep_list;
    return hsp_id;
    }

static void GAM_Result_geneseed_filter(GAM *gam, Comparison *comparison){
    register gint i, hsp_id = 0;
    register HSP *hsp;
    GAM_Geneseed_Data gsd;
    gsd.hsp_list = g_ptr_array_new();
    gsd.rangetree = RangeTree_create();
    gsd.max_cobs_hsp = NULL;
    gsd.gam = gam;
    /* Build the rangetree */
    GAM_Result_geneseed_add(comparison->dna_hspset, &gsd);
    GAM_Result_geneseed_add(comparison->protein_hspset, &gsd);
    GAM_Result_geneseed_add(comparison->codon_hspset, &gsd);
    g_assert(gsd.hsp_list->len);
    g_assert(gsd.max_cobs_hsp);
    gsd.fwd_keep = g_new0(gboolean, gsd.hsp_list->len);
    gsd.rev_keep = g_new0(gboolean, gsd.hsp_list->len);
    /* Find HSPs reachable from geneseed HSPs */
    for(i = 0; i < gsd.hsp_list->len; i++){
        hsp = gsd.hsp_list->pdata[i];
        if(hsp->score >= hsp->hsp_set->param->has->geneseed_threshold){
            GAM_Result_geneseed_visit_hsp_fwd(&gsd, i);
            GAM_Result_geneseed_visit_hsp_rev(&gsd, i);
            }
        }
    /* Keep HSPs marked as keep */
    hsp_id = GAM_Result_geneseed_select(comparison->dna_hspset,
                 gsd.fwd_keep, gsd.rev_keep, hsp_id);
    hsp_id = GAM_Result_geneseed_select(comparison->protein_hspset,
                 gsd.fwd_keep, gsd.rev_keep, hsp_id);
    hsp_id = GAM_Result_geneseed_select(comparison->codon_hspset,
                 gsd.fwd_keep, gsd.rev_keep, hsp_id);
    /**/
    if(comparison->dna_hspset
    && (!comparison->dna_hspset->hsp_list->len)){
        HSPset_destroy(comparison->dna_hspset);
        comparison->dna_hspset = NULL;
        }
    if(comparison->protein_hspset
    && (!comparison->protein_hspset->hsp_list->len)){
        HSPset_destroy(comparison->protein_hspset);
        comparison->protein_hspset = NULL;
        }
    if(comparison->codon_hspset
    && (!comparison->codon_hspset->hsp_list->len)){
        HSPset_destroy(comparison->codon_hspset);
        comparison->codon_hspset = NULL;
        }
    /**/
    /**/
    /*
    g_message("start with [%d] FILTER down to [%d]",
            gsd.hsp_list->len,
        (comparison->dna_hspset?comparison->dna_hspset->hsp_list->len:0)
       +(comparison->protein_hspset?comparison->protein_hspset->hsp_list->len:0)
       +(comparison->codon_hspset?comparison->codon_hspset->hsp_list->len:0));
       */
    RangeTree_destroy(gsd.rangetree, NULL, NULL);
    g_free(gsd.fwd_keep);
    g_free(gsd.rev_keep);
    g_ptr_array_free(gsd.hsp_list, TRUE);
    return;
    }

GAM_Result *GAM_Result_heuristic_create(GAM *gam,
                                        Comparison *comparison){
    register GAM_Result *gam_result;
    register HSPset_ArgumentSet *has
        = Comparison_Param_get_HSPSet_Argument_Set(comparison->param);
    if(has->geneseed_threshold){
        /* Raise score threshold to geneseed
         * to prevent low-scoring subopt alignments.
         */
        GAM_lock(gam);
        if(gam->gas->threshold < has->geneseed_threshold)
            gam->gas->threshold = has->geneseed_threshold;
        GAM_unlock(gam);
        GAM_Result_geneseed_filter(gam, comparison);
        }
    if(!Comparison_has_hsps(comparison))
        return NULL;
    /*
    g_message("heuristic create with [%d,%d,%d]",
        comparison->dna_hspset?comparison->dna_hspset->hsp_list->len:0,
        comparison->protein_hspset?comparison->protein_hspset->hsp_list->len:0,
        comparison->codon_hspset?comparison->codon_hspset->hsp_list->len:0);
    Comparison_print(comparison);
    */
    if(gam->gas->use_gapped_extension)
        gam_result = GAM_Result_SDP_create(gam, comparison);
    else
        gam_result = GAM_Result_BSDP_create(gam, comparison);
    return gam_result;
    }

/**/

GAM_Result *GAM_Result_exhaustive_create(GAM *gam,
                                         Sequence *query,
                                         Sequence *target){
    register Alignment *alignment;
    register C4_Score threshold;
    register GAM_Result *gam_result;
    register OPair *opair;
    g_assert(gam->optimal);
    GAM_lock(gam);
    Sequence_lock(query);
    Sequence_lock(target);
    gam_result = GAM_Result_create(gam, query, target);
    Sequence_unlock(query);
    Sequence_unlock(target);
    opair = OPair_create(gam->optimal, gam_result->subopt,
                         query->len, target->len, gam_result->user_data);
    GAM_unlock(gam);
    if(gam->verbosity > 1)
        g_message("Exhaustive alignment of [%s] [%s]",
                  query->id, target->id);
    /**/
    do {
        GAM_lock(gam);
        threshold = GAM_get_query_threshold(gam, query);
        GAM_unlock(gam);
        alignment = OPair_next_path(opair, threshold);
        if(!alignment)
            break;
        g_assert(alignment->score >= threshold);
        GAM_Result_add_alignment(gam_result, alignment, threshold);
        if(gam->verbosity > 2)
            g_message("Found alignment number [%d] score [%d]",
                gam_result->alignment_list->len, alignment->score);
    } while(gam->gas->use_subopt);
    OPair_destroy(opair);
    if(!gam_result->alignment_list){ /* No alignments */
        GAM_Result_destroy(gam_result);
        return NULL;
        }
    return gam_result;
    }

GAM_Result *GAM_Result_share(GAM_Result *gam_result){
    gam_result->ref_count++;
    return gam_result;
    }

void GAM_Result_destroy(GAM_Result *gam_result){
    register gint i;
    if(--gam_result->ref_count)
        return;
    if(gam_result->alignment_list){
        for(i = 0; i < gam_result->alignment_list->len; i++)
            Alignment_destroy(gam_result->alignment_list->pdata[i]);
        g_ptr_array_free(gam_result->alignment_list, TRUE);
        }
    Model_Type_destroy_data(gam_result->gam->gas->type,
                            gam_result->user_data);
    Model_Type_destroy_data(gam_result->gam->gas->type,
                            gam_result->self_data);
    Sequence_destroy(gam_result->query);
    Sequence_destroy(gam_result->target);
    GAM_lock(gam_result->gam);
    GAM_destroy(gam_result->gam);
    SubOpt_destroy(gam_result->subopt);
    GAM_unlock(gam_result->gam);
    g_free(gam_result);
    return;
    }

static void GAM_display_alignment(GAM *gam, Alignment *alignment,
            Sequence *query, Sequence *target,
            gint result_id, gint rank,
            gpointer user_data, gpointer self_data, FILE *fp){
    if(gam->gas->show_alignment)
        Alignment_display(alignment, query, target,
                          gam->dna_submat, gam->protein_submat,
                          gam->translate, fp);
    if(gam->gas->show_sugar)
        Alignment_display_sugar(alignment, query, target, fp);
    if(gam->gas->show_cigar)
        Alignment_display_cigar(alignment, query, target, fp);
    if(gam->gas->show_vulgar)
        Alignment_display_vulgar(alignment, query, target, fp);
    if(gam->gas->show_query_gff)
        Alignment_display_gff(alignment, query, target, gam->translate,
                              TRUE, FALSE, result_id, user_data, fp);
    if(gam->gas->show_target_gff)
        Alignment_display_gff(alignment, query, target, gam->translate,
             FALSE, Model_Type_has_genomic_target(gam->gas->type),
             result_id, user_data, fp);
    if(gam->gas->ryo)
        Alignment_display_ryo(alignment, query, target,
                              gam->gas->ryo, gam->translate, rank,
                              user_data, self_data, fp);
    fflush(fp);
    return;
    }

static void GAM_Result_display(GAM_Result *gam_result){
    register gint i;
    register Alignment *alignment;
    g_assert(gam_result);
    for(i = 0; i < gam_result->alignment_list->len; i++){
        alignment = gam_result->alignment_list->pdata[i];
        GAM_display_alignment(gam_result->gam, alignment,
                gam_result->query, gam_result->target,
                i+1, 0, gam_result->user_data, gam_result->self_data, stdout);
        }
    return;
    }

void GAM_Result_submit(GAM_Result *gam_result){
    register GAM_QueryResult *gqr;
    g_assert(gam_result);
    GAM_lock(gam_result->gam);
    if(gam_result->gam->gas->best_n){
        GAM_QueryResult gqr_lookup = {.query_id = gam_result->query->id};
        void *tree_node = tfind((void*)&gqr_lookup,
                                      &gam_result->gam->bestn_tree,
                                      GAM_compare_id);
        gqr = tree_node ? *(GAM_QueryResult **)tree_node : NULL;
        if(!gqr){
            gqr = GAM_QueryResult_create(gam_result->gam,
                                         gam_result->query->id);
            tsearch((void*)gqr, &gam_result->gam->bestn_tree, GAM_compare_id);
            }
        g_assert(gqr);
        g_assert(!strcmp(gqr->query_id, gam_result->query->id));
        GAM_QueryResult_submit(gqr, gam_result);
    } else {
        GAM_Result_display(gam_result);
        }
    GAM_unlock(gam_result->gam);
    return;
    }

void GAM_lock(GAM *gam){
#ifdef USE_PTHREADS
    pthread_mutex_lock(&gam->gam_lock);
#endif /* USE_PTHREADS */
    return;
    }

void GAM_unlock(GAM *gam){
#ifdef USE_PTHREADS
    pthread_mutex_unlock(&gam->gam_lock);
#endif /* USE_PTHREADS */
    return;
    }

/**/

