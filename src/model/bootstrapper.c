/****************************************************************\
*                                                                *
*  bootstrapper: required prior to static model linking          *
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

#include <stdio.h>  /* For remove() */
#include <stdlib.h> /* For system() */
#include <string.h> /* For strstr() */
#include <strings.h> /* For strcasecmp() */

#include "ungapped.h"
#include "affine.h"
#include "est2genome.h"
#include "protein2dna.h"
#include "protein2genome.h"
#include "ner.h"
#include "coding2coding.h"
#include "coding2genome.h"
#include "cdna2genome.h"
#include "genome2genome.h"
#include "modeltype.h"

#include "submat.h"
#include "splice.h"

#include "optimal.h"
#include "heuristic.h"
#include "codegen.h"
#include "argument.h"
#include "sdp.h"

typedef struct {
          gchar *archive_path;
          gchar *header_path;
          gchar *lookup_path;
          gchar *object_path;
           gint  archive_member_count;
          FILE  *lookup_fp;
          FILE  *header_fp;
         GTree  *model_dictionary;
   GStringChunk *model_name_chunk;
} Bootstrapper;

static void Bootstrapper_destroy(Bootstrapper *bs){
    g_tree_destroy(bs->model_dictionary);
    g_string_chunk_free(bs->model_name_chunk);
    g_free(bs->header_path);
    g_free(bs->lookup_path);
    g_free(bs->object_path);
    g_free(bs->archive_path);
    g_free(bs);
    return;
    }

static gint Bootstrapper_model_name_compare(gconstpointer a,
                                            gconstpointer b){
    register const gchar *id_a = a, *id_b = b;
    return strcmp(id_a, id_b);
    }

static Bootstrapper *Bootstrapper_create(void){
    register Bootstrapper *bs = g_new(Bootstrapper, 1);
    g_message("Starting building bootstrapping components");
    bs->header_path = g_strdup("c4_model_archive.h"),
    bs->lookup_path = g_strdup("c4_model_lookup.c");
    bs->object_path = g_strdup("c4_model_lookup.o");
    bs->archive_path = g_strdup("c4_model_archive.a");
    bs->archive_member_count = 0;
    bs->header_fp = fopen(bs->header_path, "w");
    if(!bs->header_fp)
        g_error("Could not open [%s] for writing", bs->header_path);
    bs->lookup_fp = fopen(bs->lookup_path, "w");
    if(!bs->lookup_fp)
        g_error("Could not open [%s] for writing", bs->lookup_path);
    fprintf(bs->lookup_fp,
            "#include <glib.h>\n"
            "#include \"%s\"\n"
            "\n"
            "gpointer Bootstrapper_lookup(gchar *name){\n",
            bs->header_path);
    fprintf(bs->header_fp,
            "#include <glib.h>\n"
            "\n"
            "#include \"c4.h\"\n"
            "#include \"viterbi.h\"\n"
            "#include \"alignment.h\"\n"
            "\n"
            "gpointer Bootstrapper_lookup(gchar *name);\n"
            "\n");
    bs->model_dictionary = g_tree_new(Bootstrapper_model_name_compare);
    bs->model_name_chunk = g_string_chunk_new(64);
    return bs;
    }

static void Bootstrapper_add_to_archive(Bootstrapper *bs,
                                        gchar *path){
    register gchar *option, *command;
    register gint ret_val;
    register gchar *key, *tmp;
    register gchar *ar = "ar",
                   *ar_init = "cSq", /* czr for OSF */
                   *ar_append = "Sq"; /* zr for OSF */
    tmp = (gchar*)g_getenv("C4_AR");
    if(tmp)
        ar = tmp;
    tmp = (gchar*)g_getenv("C4_AR_INIT");
    if(tmp)
        ar_init = tmp;
    tmp = (gchar*)g_getenv("C4_AR_APPEND");
    if(tmp)
        ar_append = tmp;
    if(g_tree_lookup(bs->model_dictionary, path)){
        g_error("Duplicate model in [%s]", path);
    } else {
        g_message("Adding: [%s]", path);
        key = g_string_chunk_insert(bs->model_name_chunk, path);
        g_tree_insert(bs->model_dictionary, key, key);
        }
    if(bs->archive_member_count){
        option = ar_append;
    } else { /* Is first */
        option = ar_init;
        if(Codegen_file_exists(bs->archive_path))
            remove(bs->archive_path);
        }
    command = g_strdup_printf("%s -%s %s %s",
                              ar, option, bs->archive_path, path);
    g_message("Adding object file [%s] to archive", path);
    g_print("%s\n", command);
    ret_val = system(command);
    if(ret_val)
        g_error("Problem adding to archive [%d][%s]", ret_val, command);
    bs->archive_member_count++;
    g_free(command);
    return;
    }

static void Bootstrapper_index_archive(Bootstrapper *bs){
    register gchar *command;
    register gint ret_val;
    command = g_strdup_printf("ranlib %s", bs->archive_path);
    g_message("Indexing archive [%s]", bs->archive_path);
    g_print("%s\n", command);
    ret_val = system(command);
    g_message("Done");
    if(ret_val)
        g_error("Problem indexing archive [%d][%s]", ret_val, command);
    bs->archive_member_count++;
    g_free(command);
    return;
    }

static void Bootstrapper_close(Bootstrapper *bs){
    register gchar *command;
    register gchar *cc_command = "gcc";
    register gchar *cc_flags = "-O2";
    register gchar *tmp;
    fprintf(bs->lookup_fp, "    return NULL;\n");
    fprintf(bs->lookup_fp, "    }\n");
    fclose(bs->header_fp);
    fclose(bs->lookup_fp);
    tmp = (gchar*)g_getenv("CC");
    if(tmp)
        cc_command = tmp;
    tmp = (gchar*)g_getenv("CFLAGS");
    if(tmp)
        cc_flags = tmp;
    command = g_strconcat(cc_command, " ", cc_flags,
                          " -I" SOURCE_ROOT_DIR "/src/c4"
                          " -I" SOURCE_ROOT_DIR "/src/sequence"
                          " -I" SOURCE_ROOT_DIR "/src/general"
                          " -I" SOURCE_ROOT_DIR "/src/struct"
                          " " GLIB_CFLAGS,
                          " -o ", bs->object_path,
                          " -c ", bs->lookup_path, NULL);
    g_message("Compiling with: [%s]", command);
    if(system(command))
        g_error("Problem compiling with: [%s]", command);
    g_free(command);
    /**/
    Bootstrapper_add_to_archive(bs, bs->object_path);
    Bootstrapper_index_archive(bs);
    remove(bs->lookup_path);
    remove(bs->object_path);
    g_message("Finished building bootstrapping components");
    return;
    }

static void Bootstrapper_add_codegen(Bootstrapper *bs,
                                     Codegen *codegen){
    g_assert(codegen);
    g_assert(codegen->name);
    Bootstrapper_add_to_archive(bs, codegen->object_path);
    g_message("Adding function [%s]", codegen->name);
    fprintf(bs->lookup_fp,
            "    if(!strcmp(name, \"%s\"))\n", codegen->name);
    fprintf(bs->lookup_fp,
            "        return %s;\n", codegen->name);
    fprintf(bs->header_fp, "C4_Score %s%s;\n",
            codegen->name, Viterbi_DP_Func_ARGS_STR);
    return;
    }

static void Bootstrapper_process_Optimal(Bootstrapper *bs,
                                         Optimal *optimal){
    register GPtrArray *codegen_list;
    register Codegen *codegen;
    register gint i;
    if(!optimal)
        return;
    g_message("Process optimal [%s]", optimal->name);
    codegen_list = Optimal_make_Codegen_list(optimal);
    for(i = 0; i < codegen_list->len; i++){
        codegen = codegen_list->pdata[i];
        Bootstrapper_add_codegen(bs, codegen);
        Codegen_destroy(codegen);
        }
    g_ptr_array_free(codegen_list, TRUE);
    return;
    }

static void Bootstrapper_process_C4_Model(Bootstrapper *bs,
                C4_Model *model, gboolean generate_bsdp_code,
                                 gboolean generate_sdp_code){
    register gint i;
    register Optimal *optimal = Optimal_create(model, NULL,
                       Optimal_Type_SCORE
                      |Optimal_Type_PATH
                      |Optimal_Type_REDUCED_SPACE,
                      TRUE);
    register Heuristic *heuristic;
    register SDP *sdp;
    register GPtrArray *codegen_list, *optimal_list;
    register Codegen *codegen;
    g_message("Bootstrapper processing model [%s]", model->name);
    Bootstrapper_process_Optimal(bs, optimal);
    if(generate_bsdp_code){
        heuristic = Heuristic_create(model);
        optimal_list = Heuristic_get_Optimal_list(heuristic);
        for(i = 0; i < optimal_list->len; i++)
            Bootstrapper_process_Optimal(bs, optimal_list->pdata[i]);
        Heuristic_destroy(heuristic);
        g_ptr_array_free(optimal_list, TRUE);
        }
    if(generate_sdp_code){
        sdp = SDP_create(model);
        codegen_list = SDP_get_codegen_list(sdp);
        for(i = 0; i < codegen_list->len; i++){
            codegen = codegen_list->pdata[i];
            Bootstrapper_add_codegen(bs, codegen);
            }
        g_ptr_array_free(codegen_list, TRUE);
        SDP_destroy(sdp);
        }
    Optimal_destroy(optimal);
    return;
    }

static void Bootstrapper_process(Bootstrapper *bs, gchar *models){
    register gchar **model_list = g_strsplit(models, " ", 1024);
    register gint i;
    register Model_Type type;
    register C4_Model *model;
    for(i = 0; model_list[i]; i++){
        g_message("Processing model [%s]", model_list[i]);
        type = Model_Type_from_string(model_list[i]);
        switch(type){
            case Model_Type_UNGAPPED:
                /* ungapped dna2dna */
                model = Ungapped_create(Match_Type_DNA2DNA);
                Bootstrapper_process_C4_Model(bs, model, FALSE, FALSE);
                C4_Model_destroy(model);
                /* ungapped protein2protein */
                model = Ungapped_create(Match_Type_PROTEIN2PROTEIN);
                Bootstrapper_process_C4_Model(bs, model, FALSE, FALSE);
                C4_Model_destroy(model);
                /* ungapped dna2protein */
                model = Ungapped_create(Match_Type_DNA2PROTEIN);
                Bootstrapper_process_C4_Model(bs, model, FALSE, FALSE);
                C4_Model_destroy(model);
                /* ungapped protein2dna */
                model = Ungapped_create(Match_Type_PROTEIN2DNA);
                Bootstrapper_process_C4_Model(bs, model, FALSE, FALSE);
                C4_Model_destroy(model);
                break;
            case Model_Type_UNGAPPED_TRANS:
                model = Ungapped_create(Match_Type_CODON2CODON);
                Bootstrapper_process_C4_Model(bs, model, FALSE, FALSE);
                C4_Model_destroy(model);
                break;
            case Model_Type_AFFINE_GLOBAL:
                model = Affine_create(Affine_Model_Type_GLOBAL,
                                      Alphabet_Type_DNA,
                                      Alphabet_Type_DNA, FALSE);
                Bootstrapper_process_C4_Model(bs, model, FALSE, FALSE);
                C4_Model_destroy(model);
                /**/
                model = Affine_create(Affine_Model_Type_GLOBAL,
                                      Alphabet_Type_PROTEIN,
                                      Alphabet_Type_PROTEIN, FALSE);
                Bootstrapper_process_C4_Model(bs, model, FALSE, FALSE);
                C4_Model_destroy(model);
                break;
            case Model_Type_AFFINE_BESTFIT:
                model = Affine_create(Affine_Model_Type_BESTFIT,
                                      Alphabet_Type_DNA,
                                      Alphabet_Type_DNA, FALSE);
                Bootstrapper_process_C4_Model(bs, model, FALSE, FALSE);
                C4_Model_destroy(model);
                /**/
                model = Affine_create(Affine_Model_Type_BESTFIT,
                                      Alphabet_Type_PROTEIN,
                                      Alphabet_Type_PROTEIN, FALSE);
                Bootstrapper_process_C4_Model(bs, model, FALSE, FALSE);
                C4_Model_destroy(model);
                break;
            case Model_Type_AFFINE_LOCAL:
                model = Affine_create(Affine_Model_Type_LOCAL,
                                      Alphabet_Type_DNA,
                                      Alphabet_Type_DNA, FALSE);
                Bootstrapper_process_C4_Model(bs, model, TRUE, TRUE);
                C4_Model_destroy(model);
                /**/
                model = Affine_create(Affine_Model_Type_LOCAL,
                                      Alphabet_Type_PROTEIN,
                                      Alphabet_Type_PROTEIN, FALSE);
                Bootstrapper_process_C4_Model(bs, model, TRUE, TRUE);
                C4_Model_destroy(model);
                break;
            case Model_Type_AFFINE_OVERLAP:
                model = Affine_create(Affine_Model_Type_OVERLAP,
                                      Alphabet_Type_DNA,
                                      Alphabet_Type_DNA, FALSE);
                Bootstrapper_process_C4_Model(bs, model, FALSE, FALSE);
                C4_Model_destroy(model);
                break;
            case Model_Type_EST2GENOME:
                model = EST2Genome_create();
                Bootstrapper_process_C4_Model(bs, model, TRUE, TRUE);
                C4_Model_destroy(model);
                break;
            case Model_Type_NER:
                model = NER_create(Alphabet_Type_DNA,
                                   Alphabet_Type_DNA);
                Bootstrapper_process_C4_Model(bs, model, TRUE, TRUE);
                C4_Model_destroy(model);
                /**/
                model = NER_create(Alphabet_Type_PROTEIN,
                                   Alphabet_Type_PROTEIN);
                Bootstrapper_process_C4_Model(bs, model, TRUE, TRUE);
                C4_Model_destroy(model);
                break;
            case Model_Type_PROTEIN2DNA:
                model = Protein2DNA_create(Affine_Model_Type_LOCAL);
                Bootstrapper_process_C4_Model(bs, model, TRUE, TRUE);
                C4_Model_destroy(model);
                break;
            case Model_Type_PROTEIN2DNA_BESTFIT:
                model = Protein2DNA_create(Affine_Model_Type_BESTFIT);
                Bootstrapper_process_C4_Model(bs, model, FALSE, FALSE);
                C4_Model_destroy(model);
                break;
            case Model_Type_PROTEIN2GENOME:
                model = Protein2Genome_create(Affine_Model_Type_LOCAL);
                Bootstrapper_process_C4_Model(bs, model, TRUE, TRUE);
                C4_Model_destroy(model);
                break;
            case Model_Type_PROTEIN2GENOME_BESTFIT:
                model = Protein2Genome_create(Affine_Model_Type_BESTFIT);
                Bootstrapper_process_C4_Model(bs, model, FALSE, FALSE);
                C4_Model_destroy(model);
                break;
            case Model_Type_CODING2CODING:
                model = Coding2Coding_create();
                Bootstrapper_process_C4_Model(bs, model, TRUE, TRUE);
                C4_Model_destroy(model);
                break;
            case Model_Type_CODING2GENOME:
                model = Coding2Genome_create();
                Bootstrapper_process_C4_Model(bs, model, FALSE, TRUE);
                C4_Model_destroy(model);
                break;
            case Model_Type_CDNA2GENOME:
                model = CDNA2Genome_create();
                Bootstrapper_process_C4_Model(bs, model, FALSE, TRUE);
                C4_Model_destroy(model);
                break;
            case Model_Type_GENOME2GENOME:
                model = Genome2Genome_create();
                Bootstrapper_process_C4_Model(bs, model, FALSE, TRUE);
                C4_Model_destroy(model);
                break;
            default:
                g_error("Unknown Model type [%s]", model_list[i]);
                break;
            }
        }
    g_strfreev(model_list);
    return;
    }

int Argument_main(Argument *arg){
    register Bootstrapper *bs = Bootstrapper_create();
    register gchar *models = (gchar*)g_getenv("C4_COMPILED_MODELS");
    Heuristic_ArgumentSet_create(arg);
    Match_ArgumentSet_create(arg);
    Affine_ArgumentSet_create(arg);
    NER_ArgumentSet_create(arg);
    Intron_ArgumentSet_create(arg);
    Codegen_ArgumentSet_create(arg);
    Argument_process(arg, "exonerate:bootstrapper", NULL, NULL);
    if((!models) || (!strcasecmp(models, "all")))
        models = "u u:t a:g a:b a:l a:o e2g ner p2d p2d:b p2g p2g:b c2c c2g cd2g";
    Bootstrapper_process(bs, models);
    Bootstrapper_close(bs);
    Bootstrapper_destroy(bs);
    return 0;
    }

