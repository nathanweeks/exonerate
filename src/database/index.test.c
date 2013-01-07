/****************************************************************\
*                                                                *
*  Library for manipulation of FASTA format databases            *
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

#include <string.h>
#include "index.h"

gint Argument_main(Argument *arg){
    register Index *index;
    register Alphabet *alphabet = Alphabet_create(Alphabet_Type_DNA, FALSE);
    register Sequence *query = Sequence_create("query", NULL,
            /* "ATGCAGCCCCGG", */
            /* "ACGTAGAGGCCGAGGATGCCTCCCATTCTCAGCT", */
            /* "ATGCAGCCCCGGGTACTCCTTGTTGTTGCCCTCCTGGCGCTCCTGGCCTCTGCCC", */
            "ATGCANCCCCNGGTACTCCTTGTTNTTNCCCTCNTGGCGNTCCTGGNCTCTNCCC",
            /*
               "CCGGNAAGCTCANCTTGGACCACCGACTCTCGANTGNNTCGCCGCGGGAGCCGGNTGGAN"
               "AACCTGAGCGGGACTGGNAGAAGGAGCAGAGGGAGGCAGCACCCGGCGTGACGGNAGTGT"
               "GTGGGGCACTCAGGCCTTCCGCAGTGTCATCTGCCACACGGAAGGCACGGCCACGGGCAG"
               "GGGGGTCTATGATCTTCTGCATGCCCAGCTGGCATGGCCCCACGTAGAGTGGNNTGGCGT"
               "CTCGGTGCTGGTCAGCGACACGTTGTCCTGGCTGGGCAGGTCCAGCTCCCGGAGGACCTG"
               "GGGCTTCAGCTTCCCGTAGCGCTGGCTGCAGTGACGGATGCTCTTGCGCTGCCATTTCTG"
               "GGTGCTGTCACTGTCCTTGCTCACTCCAAACCAGTTCGGCGGTCCCCCTGCGGATGGTCT"
               "GTGTTGATGGACGTTTGGGCTTTGCAGCACCGGCCGCCGAGTTCATGGTNGGGTNAAGAG"
               "ATTTGGGTTTTTTCN",
               */
            0, Sequence_Strand_FORWARD, alphabet);
    register Match *match;
    register HSP_Param *hsp_param;
    register ArgumentSet *as = ArgumentSet_create("Input options");
    register GPtrArray *index_hsp_set_list;
    register Index_HSPset *index_hsp_set;
    register gint i;
    gchar *path;
    ArgumentSet_add_option(as, 'i', "index", "path",
            "exonerate sequence index file (.esi)", "none",
            Argument_parse_string, &path);
    Argument_absorb_ArgumentSet(arg, as);
    HSPset_ArgumentSet_create(arg);
    Match_ArgumentSet_create(arg);
    Argument_process(arg, "index.test", NULL, NULL);
    if(!strcmp(path, "none")){ /* To ensure 'make check' does not fail */
        g_warning("No path set for test index file");
        return 0;
        }
    match = Match_find(Match_Type_DNA2DNA);
    hsp_param = HSP_Param_create(match, FALSE);
    /**/
    index = Index_open(path);
    index_hsp_set_list = Index_get_HSPsets(index, hsp_param, query, FALSE);
    if(index_hsp_set_list){
        for(i = 0; i < index_hsp_set_list->len; i++){
            index_hsp_set = index_hsp_set_list->pdata[i];
            HSPset_print(index_hsp_set->hsp_set);
            Index_HSPset_destroy(index_hsp_set);
            }
        g_ptr_array_free(index_hsp_set_list, TRUE);
        }
    Index_destroy(index);
    /**/
    HSP_Param_destroy(hsp_param);
    Alphabet_destroy(alphabet);
    Sequence_destroy(query);
    return 0;
    }

