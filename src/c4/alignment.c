/****************************************************************\
*                                                                *
*  C4 dynamic programming library - alignment code               *
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
#include <ctype.h>  /* For tolower() */
#include <time.h>   /* For time() */

#include "match.h"
#include "alignment.h"
#include "exonerate_util.h"

Alignment_ArgumentSet *Alignment_ArgumentSet_create(Argument *arg){
    register ArgumentSet *as;
    static Alignment_ArgumentSet aas = {80, TRUE};
    if(arg){
        as = ArgumentSet_create("Alignment options");
        ArgumentSet_add_option(as, '\0', "alignmentwidth", NULL,
            "Alignment display width", "80", Argument_parse_int,
            &aas.alignment_width);
        ArgumentSet_add_option(as, '\0', "forwardcoordinates", NULL,
            "Report all coordinates on the forward strand", "TRUE",
            Argument_parse_boolean, &aas.forward_strand_coords);
        Argument_absorb_ArgumentSet(arg, as);
        }
    return &aas;
    }

/**/

Alignment *Alignment_create(C4_Model *model, Region *region,
                            C4_Score score){
    register Alignment *alignment = g_new(Alignment, 1);
    g_assert(model);
    g_assert(region);
    alignment->ref_count = 1;
    alignment->region = Region_share(region);
    alignment->score = score;
    alignment->operation_list = g_ptr_array_new();
    alignment->model = C4_Model_share(model);
    return alignment;
    }

Alignment *Alignment_share(Alignment *alignment){
    g_assert(alignment);
    alignment->ref_count++;
    return alignment;
    }

void Alignment_destroy(Alignment *alignment){
    register gint i;
    g_assert(alignment);
    if(--alignment->ref_count)
        return;
    for(i = 0; i < alignment->operation_list->len; i++)
        g_free(alignment->operation_list->pdata[i]);
    g_ptr_array_free(alignment->operation_list, TRUE);
    Region_destroy(alignment->region);
    C4_Model_destroy(alignment->model);
    g_free(alignment);
    return;
    }

void Alignment_add(Alignment *alignment, C4_Transition *transition,
                   gint length){
    register AlignmentOperation *alignment_operation, *prev;
    g_assert(alignment);
    g_assert(transition);
    if(alignment->operation_list->len){
        prev = alignment->operation_list->pdata
              [alignment->operation_list->len-1];
        if(prev->transition == transition){
            prev->length += length;
            g_assert(prev->length >= 0); /* Cannot have -ve total */
            if(prev->length == 0){ /* Drop the transition */
                g_free(prev);
                g_ptr_array_set_size(alignment->operation_list,
                                     alignment->operation_list->len-1);
                }
            return;
            }
        g_assert(prev->transition->output == transition->input);
        }
    alignment_operation = g_new(AlignmentOperation, 1);
    alignment_operation->transition = transition;
    alignment_operation->length = length;
    g_ptr_array_add(alignment->operation_list, alignment_operation);
    return;
    }

/**/

static gchar Alignment_match_get_symbol(Sequence *seq,
                                        gint pos, gint advance,
                                        Translate *translate){
    register gchar symbol = '\0';
    switch(advance){
        case 1:
            symbol = Sequence_get_symbol(seq, pos);
            break;
        case 3:
            g_assert(translate);
            symbol = Translate_base(translate,
                               Sequence_get_symbol(seq, pos),
                               Sequence_get_symbol(seq, pos+1),
                               Sequence_get_symbol(seq, pos+2));
            break;
        default:
            g_error("Cannot process match advance of [%d]", advance);
            break;
        }
    return symbol;
    }

static gchar *Alignment_match_get_string(Sequence *seq,
                       gint pos, gint advance, gint max,
                       Translate *translate){
    register gint ch;
    switch(max){
        case 1:
            g_assert(advance == 1);
            return g_strdup_printf("%c", Sequence_get_symbol(seq, pos));
            break;
        case 3:
            switch(advance){
                case 1:
                    g_assert(seq->alphabet->type == Alphabet_Type_PROTEIN);
                    ch = Sequence_get_symbol(seq, pos);
                    return g_strdup(Alphabet_aa2tla(ch));
                    break;
                case 3:
                    g_assert(seq->alphabet->type == Alphabet_Type_DNA);
                    return g_strdup_printf("%c%c%c", /* codon */
                                    Sequence_get_symbol(seq, pos),
                                    Sequence_get_symbol(seq, pos+1),
                                    Sequence_get_symbol(seq, pos+2));
                    break;
                default:
                    g_error("Cannot handle advance [%d]", advance);
                    break;
                }
            break;
        default:
            g_error("Cannot find match display for max advance of [%d]",
                    advance);
            break;
        }
    return NULL;
    }

/**/

static gchar Alignment_get_gene_orientation(Alignment *alignment){
    register gint i;
    register AlignmentOperation *ao;
    for(i = 0; i < alignment->operation_list->len; i++){
        ao = alignment->operation_list->pdata[i];
        if(ao->transition->label == C4_Label_5SS)
            return '+';
        if(ao->transition->label == C4_Label_3SS)
            return '-';
        }
    return '.';
    }

static gint Alignment_get_coordinate(Alignment *alignment,
                                     Sequence *query,
                                     Sequence *target,
                                     gboolean on_query,
                                     gboolean report_start){
    register gint pos;
    register Alignment_ArgumentSet *aas
        = Alignment_ArgumentSet_create(NULL);
    if(on_query){
        if(report_start){
            pos = alignment->region->query_start;
        } else {
            pos = Region_query_end(alignment->region);
            }
        if(aas->forward_strand_coords
        && (query->strand == Sequence_Strand_REVCOMP)){
            pos = query->len - pos;
            }
    } else { /* on_target */
        if(report_start){
            pos = alignment->region->target_start;
        } else {
            pos = Region_target_end(alignment->region);
            }
        if(aas->forward_strand_coords
        && (target->strand == Sequence_Strand_REVCOMP)){
            pos = target->len - pos;
            }
        }
    return pos;
    }

static gint Alignment_convert_coordinate(Alignment *alignment,
                                         Sequence *query,
                                         Sequence *target,
                                         gint query_pos,
                                         gint target_pos,
                                         gboolean on_query){
    register gint pos;
    register Alignment_ArgumentSet *aas
        = Alignment_ArgumentSet_create(NULL);
    if(on_query){
        pos = query_pos;
        if(aas->forward_strand_coords
        && (query->strand == Sequence_Strand_REVCOMP)){
            pos = query->len - pos;
            }
    } else { /* on_target */
        pos = target_pos;
        if(aas->forward_strand_coords
        && (target->strand == Sequence_Strand_REVCOMP)){
            pos = target->len - pos;
            }
        }
    return pos;
    }

static gint Alignment_get_max_pos_len(Alignment *alignment,
                          Sequence *query, Sequence *target){
    register gint qmax = MAX(
      Alignment_get_coordinate(alignment, query, target, TRUE, TRUE),
      Alignment_get_coordinate(alignment, query, target, TRUE, FALSE));
    register gint tmax = MAX(
      Alignment_get_coordinate(alignment, query, target, FALSE, TRUE),
      Alignment_get_coordinate(alignment, query, target, FALSE, FALSE));
    register gint max = MAX(qmax, tmax);
    register gchar *tmpstr = g_strdup_printf("%d", max);
    register gint maxlen = strlen(tmpstr);
    g_free(tmpstr);
    return maxlen;
    }

/**/

typedef struct {
    gint query_pos;
    gint target_pos;
} AlignmentPosition;

typedef struct {
    gint query_separation;
    gint target_separation;
} AlignmentSeparation;

typedef struct {
      GString *outer_query;
      GString *inner_query;
      GString *middle;
      GString *inner_target;
      GString *outer_target;
    GPtrArray *row_marker; /* List containing AlignmentPosition */
         gint  max_pos_len;
         gint  width;
         gint  limit;
         gint  query_intron_count;
         gint  target_intron_count;
         gint  joint_intron_count;
         gint  intron_advance_query;
         gint  intron_advance_target;
        gchar  gene_orientation;
         gint  ner_count;
         gint  ner_advance_query;
         gint  ner_advance_target;
         gint  curr_split_codon_count; /* Add 2 for each split codon */
    GPtrArray *split_codon_separation_list;
               /* List containing AlignmentSeparation */
} AlignmentView;

static AlignmentView *AlignmentView_create(Alignment *alignment,
                                           Sequence *query,
                                           Sequence *target){
    register AlignmentView *av = g_new(AlignmentView, 1);
    register AlignmentOperation *ao;
    register AlignmentSeparation *curr_split_codon_separation = NULL;
    register gint i;
    register Alignment_ArgumentSet *aas
           = Alignment_ArgumentSet_create(NULL);
    av->outer_query  = g_string_sized_new(16);
    if(alignment->model->max_query_advance == 3)
        av->inner_query = g_string_sized_new(16);
    else
        av->inner_query = NULL;
    av->middle = g_string_sized_new(16);
    if(alignment->model->max_target_advance == 3)
        av->inner_target = g_string_sized_new(16);
    else
        av->inner_target = NULL;
    av->outer_target  = g_string_sized_new(16);
    av->row_marker = g_ptr_array_new();
    av->max_pos_len = Alignment_get_max_pos_len(alignment,
                                                query, target);
    av->width = aas->alignment_width-((av->max_pos_len + 5) << 1);
    g_assert(av->width > 0);
    av->limit = av->width;
    av->query_intron_count = 0;
    av->target_intron_count = 0;
    av->joint_intron_count = 0;
    av->intron_advance_query = 0;
    av->intron_advance_target = 0;
    av->gene_orientation = Alignment_get_gene_orientation(alignment);
    av->ner_count = 0;
    av->ner_advance_query = 0;
    av->ner_advance_target = 0;
    av->curr_split_codon_count = 0;
    av->split_codon_separation_list = g_ptr_array_new();
    for(i = 0; i < alignment->operation_list->len; i++){
        ao = alignment->operation_list->pdata[i];
        if(curr_split_codon_separation){
            if(ao->transition->label == C4_Label_SPLIT_CODON){
                g_ptr_array_add(av->split_codon_separation_list,
                                curr_split_codon_separation);
                curr_split_codon_separation = NULL;
            } else {
                curr_split_codon_separation->query_separation
                    += (ao->length * ao->transition->advance_query);
                curr_split_codon_separation->target_separation
                    += (ao->length * ao->transition->advance_target);
                }
        } else {
            if(ao->transition->label == C4_Label_SPLIT_CODON){
                curr_split_codon_separation
                        = g_new(AlignmentSeparation, 1);
                curr_split_codon_separation->query_separation
                    = (ao->length * ao->transition->advance_query);
                curr_split_codon_separation->target_separation
                    = (ao->length * ao->transition->advance_target);
                }
            }
        }
    /* Check we have an even number of split codons */
    g_assert(!curr_split_codon_separation);
    return av;
    }

static void AlignmentView_destroy(AlignmentView *av){
    register gint i;
    for(i = 0; i < av->split_codon_separation_list->len; i++)
        g_free(av->split_codon_separation_list->pdata[i]);
    g_ptr_array_free(av->split_codon_separation_list, TRUE);
    g_string_free(av->outer_query,  TRUE);
    g_string_free(av->middle, TRUE);
    g_string_free(av->outer_target,  TRUE);

    if(av->inner_query)
        g_string_free(av->inner_query, TRUE);
    if(av->inner_target)
        g_string_free(av->inner_target, TRUE);

    for(i = 0; i < av->row_marker->len; i++)
        g_free(av->row_marker->pdata[i]);
    g_ptr_array_free(av->row_marker, TRUE);
    g_free(av);
    return;
    }

static void AlignmentView_add(AlignmentView *av,
                              gchar *query_string,
                              gchar *inner_query_string,
                              gchar *match_string,
                              gchar *inner_target_string,
                              gchar *target_string,
                              gint query_pos, gint target_pos){
    register AlignmentPosition *apos;
    register gint i;
    g_assert(strlen(query_string) == strlen(match_string));
    g_assert(strlen(match_string) == strlen(target_string));
    if(av->inner_query){
        if(inner_query_string){
            g_assert(strlen(match_string)
                  == strlen(inner_query_string));
            g_string_append(av->inner_query, inner_query_string);
        } else {
            for(i = strlen(match_string)-1; i >=0; i--)
                g_string_append_c(av->inner_query, ' ');
            }
        }
    if(av->inner_target){
        if(inner_target_string){
            g_assert(strlen(match_string)
                  == strlen(inner_target_string));
            g_string_append(av->inner_target, inner_target_string);
        } else {
            for(i = strlen(match_string)-1; i >=0; i--)
                g_string_append_c(av->inner_target, ' ');
            }
        }
    g_string_append(av->outer_query,  query_string);
    g_string_append(av->middle, match_string);
    g_string_append(av->outer_target,  target_string);
    if(av->outer_query->len >= av->limit){
        apos = g_new(AlignmentPosition, 1);
        apos->query_pos  = query_pos;
        apos->target_pos = target_pos;
        g_ptr_array_add(av->row_marker, apos);
        av->limit += av->width;
        }
    return;
    }

typedef struct {
    gchar *match_string;
    gchar *codon;
} AlignmentMatchData;

static void Alignment_match_translate_reverse(gchar *dna, gint length,
                                              gpointer user_data){
    register AlignmentMatchData *amd = user_data;
    register gint i;
    for(i = 0; i < 3; i++)
        if(dna[i] == amd->codon[i])
            amd->match_string[i] = '!';
    return;
    }

static gchar Alignment_get_equiv_symbol(gchar symbol_a, gchar symbol_b,
                                        Submat *submat){
    register gint score;
    g_assert(symbol_a);
    g_assert(symbol_b);
    if(submat){
        score = Submat_lookup(submat, symbol_a, symbol_b);
        if(score == 0)
            return '.';
        if(score > 0){
            if(toupper(symbol_a) == toupper(symbol_b))
                return '|';
            else
                return ':';
            }
    } else {
        if(symbol_a == symbol_b)
            return '|';
        }
    return ' ';
    }

static gchar *Alignment_get_codon_match_string(gchar *codon, gchar aa,
                        Submat *protein_submat, Translate *translate){
    register gchar codon_aa = Translate_codon(translate, codon);
    register gchar match_symbol;
    register gchar *codon_match;
    gchar aa_seq[2];
    AlignmentMatchData amd;
    g_assert(translate);
    g_assert(protein_submat);
    match_symbol = Alignment_get_equiv_symbol(codon_aa, aa,
                                              protein_submat);
    codon_match = g_strnfill(3, match_symbol);
    if(match_symbol != '|'){
        amd.match_string = codon_match;
        amd.codon = codon;
        aa_seq[0] = aa;
        aa_seq[1] = '\0';
        Translate_reverse(translate, aa_seq, 1,
                          Alignment_match_translate_reverse, &amd);
        }
    return codon_match;
    }

static void AlignmentView_add_MATCH(AlignmentView *av,
                C4_Transition *transition,
                gint total_length, Sequence *query, Sequence *target,
                gint query_pos, gint target_pos,
                Submat *dna_submat, Submat *protein_submat,
                Translate *translate){
    register gchar query_symbol, target_symbol;
    register gchar *query_string, *target_string;
    register gchar *inner_query_string = NULL,
                   *inner_target_string = NULL;
    register gint max_advance;
    register gint i;
    register gint curr_query_pos = query_pos,
                  curr_target_pos = target_pos;
    register Match *match = (Match*)transition->label_data;
    gchar match_string[4]; /* FIXME: temp: should be max_advance long */
    for(i = 0; i < total_length; i++){
        max_advance = MAX(transition->advance_query,
                          transition->advance_target);
        query_string = Alignment_match_get_string(query,
                           curr_query_pos, transition->advance_query,
                           max_advance, translate);
        target_string = Alignment_match_get_string(target,
                           curr_target_pos, transition->advance_target,
                           max_advance, translate);
        query_symbol = Alignment_match_get_symbol(query,
                           curr_query_pos, transition->advance_query,
                           translate);
        target_symbol = Alignment_match_get_symbol(target,
                           curr_target_pos, transition->advance_target,
                           translate);
        if(transition->advance_query == 3)
            inner_query_string = Alphabet_aa2tla(query_symbol);
        if(transition->advance_target == 3)
            inner_target_string = Alphabet_aa2tla(target_symbol);
        if(match){
            match->display_func(match, query, target,
                                curr_query_pos, curr_target_pos,
                                match_string);
        } else { /* In the absense of Match */
            g_assert(transition->advance_query == 1);
            g_assert(transition->advance_target == 1);
            match_string[0] = (Sequence_get_symbol(query, query_pos)
                            == Sequence_get_symbol(target, target_pos))
                            ?'|':' ';
            match_string[1] = '\0';
            }
        AlignmentView_add(av, query_string, inner_query_string,
                              match_string,
                              inner_target_string, target_string,
                              curr_query_pos, curr_target_pos);
        g_free(query_string);
        g_free(target_string);
        curr_query_pos += transition->advance_query;
        curr_target_pos += transition->advance_target;
        }
    return;
    }

static void AlignmentView_add_GAP(AlignmentView *av,
                gint advance_query, gint advance_target,
                gint total_length, Sequence *query, Sequence *target,
                gint query_pos, gint target_pos, Translate *translate){
    register gint i, j;
    register gint curr_query_pos = query_pos,
                  curr_target_pos = target_pos;
    register gchar *seq_string = g_new(gchar, 4),
                   *match_string = g_new(gchar, 4),
                   *gap_string = g_new(gchar, 4);
    register Alphabet_Type emitted_alphabet_type;
    register gchar tr_codon, *tr_codon_name;
    register gboolean is_translating
        = ( ( (query->alphabet->type == Alphabet_Type_PROTEIN)
           && (target->alphabet->type == Alphabet_Type_DNA))
          ||( (query->alphabet->type == Alphabet_Type_DNA)
           && (target->alphabet->type == Alphabet_Type_PROTEIN))
          || ((advance_query|advance_target) == 3));
    g_assert(!(advance_query && advance_target));
    match_string[0] = ' ';
    gap_string[0] = '-';
    if(advance_query)
        emitted_alphabet_type = query->alphabet->type;
    else /* advance_target */
        emitted_alphabet_type = target->alphabet->type;
    for(i = 0; i < total_length; i++){
        for(j = 0; j < (advance_query|advance_target); j++){
            if(advance_query)
                seq_string[j] = Sequence_get_symbol(query, curr_query_pos+j);
            else /* advance_target */
                seq_string[j] = Sequence_get_symbol(target, curr_target_pos+j);
            match_string[j] = match_string[0];
            gap_string[j] = gap_string[0];
            }
        seq_string[j] = match_string[j] = gap_string[j] = '\0';
        tr_codon_name = NULL;
        if(is_translating){
            if(emitted_alphabet_type == Alphabet_Type_PROTEIN){
                strncpy(seq_string, Alphabet_aa2tla(seq_string[0]), 3);
                match_string[2] = match_string[1] = match_string[0];
                gap_string[2] = gap_string[1] = gap_string[0];
                seq_string[3] = match_string[3] = gap_string[3] = '\0';
                }
            if((advance_query|advance_target) == 3){
                gap_string[0] = '<';
                gap_string[1] = '-';
                gap_string[2] = '>';
                gap_string[3] = '\0';
                tr_codon = Translate_codon(translate, seq_string);
                tr_codon_name = Alphabet_aa2tla(tr_codon);
                }
            }
        if(advance_query)
            AlignmentView_add(av, seq_string,
                                  tr_codon_name,
                                  match_string,
                                  is_translating?gap_string:NULL,
                                  gap_string,
                                  curr_query_pos, curr_target_pos);
        else
            AlignmentView_add(av, gap_string,
                                  is_translating?gap_string:NULL,
                                  match_string,
                                  tr_codon_name,
                                  seq_string,
                                  curr_query_pos, curr_target_pos);
        curr_query_pos += advance_query;
        curr_target_pos += advance_target;
        }
    g_free(seq_string);
    g_free(match_string);
    g_free(gap_string);
    return;
    }
/* Display gap:[- N] codon[<->   NNN] frameshift [-*N]
 */

static void AlignmentView_set_consensus_ss_string(AlignmentView *av,
    gboolean is_5_prime, gchar *splice_site, gchar *consensus_string){
    register gchar cons_a, cons_b;
    if(av->gene_orientation == '+'){
        if(is_5_prime){ /* FWD: gt..ag */
            cons_a = 'G';
            cons_b = 'T';
        } else {
            cons_a = 'A';
            cons_b = 'G';
            }
    } else {
        g_assert(av->gene_orientation == '-');
        if(is_5_prime){ /* REV: ct..ac */
            cons_a = 'A';
            cons_b = 'C';
        } else {
            cons_a = 'C';
            cons_b = 'T';
            }
        }
    if(toupper(splice_site[0]) == cons_a)
        consensus_string[0] = '+';
    else
        consensus_string[0] = '-';
    if(toupper(splice_site[1]) == cons_b)
        consensus_string[1] = '+';
    else
        consensus_string[1] = '-';
    return;
    }

static void AlignmentView_add_SPLICE_SITE(AlignmentView *av,
                gint advance_query, gint advance_target,
                gint total_length, Sequence *query, Sequence *target,
                gint query_pos, gint target_pos, gboolean is_5_prime,
                C4_Transition *last_match){
    register gchar *gap_string = "  ";
    static gchar qy_seq_string[3], tg_seq_string[3],
                 qy_cons_string[3], tg_cons_string[3];
    qy_cons_string[0] = qy_cons_string[1] = ' ';
    tg_cons_string[0] = tg_cons_string[1] = ' ';
    qy_cons_string[2] = tg_cons_string[2] = '\0';
    qy_seq_string[2] = tg_seq_string[2] = '\0';
    if(advance_query == 2){
        qy_seq_string[0] = Sequence_get_symbol(query, query_pos);
        qy_seq_string[1] = Sequence_get_symbol(query, query_pos+1);
        AlignmentView_set_consensus_ss_string(av, is_5_prime,
                                qy_seq_string, qy_cons_string);
        strdown(qy_seq_string);
        }
    if(advance_target == 2){
        tg_seq_string[0] = Sequence_get_symbol(target, target_pos);
        tg_seq_string[1] = Sequence_get_symbol(target, target_pos+1);
        AlignmentView_set_consensus_ss_string(av, is_5_prime,
                                tg_seq_string, tg_cons_string);
        strdown(tg_seq_string);
        }
    if(advance_query == 2){
        if(advance_target == 2){ /* Joint intron */
            AlignmentView_add(av, qy_seq_string,
                          qy_cons_string, gap_string, tg_cons_string,
                          tg_seq_string,
                          query_pos, target_pos);
        } else { /* Query intron */
            g_assert(advance_target == 0);
            strdown(qy_seq_string);
            g_assert(last_match);
            if(last_match->advance_query == 3)
                AlignmentView_add(av, qy_seq_string, qy_cons_string, gap_string,
                                  gap_string, gap_string, query_pos, target_pos);
            else
                AlignmentView_add(av, qy_seq_string,
                                  NULL, qy_cons_string, NULL,
                                  gap_string, query_pos, target_pos);
            }
    } else { /* Target intron */
        g_assert(advance_query == 0);
        g_assert(advance_target == 2);
        strdown(tg_seq_string);
        g_assert(last_match);
        if(last_match->advance_target == 3)
            AlignmentView_add(av, gap_string, gap_string,
                              gap_string, tg_cons_string, tg_seq_string,
                              query_pos, target_pos);
        else
            AlignmentView_add(av, gap_string,
                              NULL, tg_cons_string, NULL,
                              tg_seq_string, query_pos, target_pos);
        }
    return;
    }

static void AlignmentView_add_INTRON(AlignmentView *av,
                gint advance_query, gint advance_target,
                Sequence *query, Sequence *target,
                gint query_pos, gint target_pos,
                C4_Transition *last_match){
    register gchar *dir_sign = "????", *label = NULL;
    register gint fill, intron_count = 0;
    register gchar *name_string, *middle_string, *gap_string,
                   *pad_string, *intron_name;
    g_assert(last_match);
    if(av->gene_orientation == '+'){
        dir_sign = ">>>>";
    } else if(av->gene_orientation == '-'){
        dir_sign = "<<<<";
        }
    if(advance_query){
        if(advance_target){
            intron_count = ++av->joint_intron_count;
            intron_name = "Joint";
            label = g_strdup_printf("%d bp // %d bp",
                        advance_query+4, advance_target+4);
        } else {
            intron_count = ++av->query_intron_count;
            intron_name = "Query";
            label = g_strdup_printf("%d bp", advance_query+4);
            }
    } else {
        intron_count = ++av->target_intron_count;
        intron_name = "Target";
        label = g_strdup_printf("%d bp", advance_target+4);
        }
    name_string = g_strdup_printf("%s %s Intron %d %s",
                      dir_sign, intron_name, intron_count, dir_sign);
    g_assert(strlen(name_string) > strlen(label));
    fill = (strlen(name_string)-strlen(label))+1;
    middle_string = g_strdup_printf("%*c%s%*c",
            ((fill|1)>>1), ' ', label, ((fill-1)>>1), ' ');
    gap_string = g_strnfill(strlen(name_string), '.');
    pad_string = g_strnfill(strlen(name_string), '^');
    g_assert(strlen(name_string) == strlen(middle_string));
    g_assert(strlen(middle_string) == strlen(gap_string));
    if(advance_query){
        if(advance_target){ /* joint intron */
            AlignmentView_add(av, name_string, NULL, middle_string, NULL,
                              name_string, query_pos, target_pos);
        } else { /* query intron */
            if(last_match->advance_query == 3)
                AlignmentView_add(av, gap_string, pad_string, middle_string,
                                  pad_string, name_string, query_pos, target_pos);
            else
                AlignmentView_add(av, gap_string, NULL, middle_string, NULL,
                                  name_string, query_pos, target_pos);
            }
    } else { /* target intron */
        if(last_match->advance_target == 3)
            AlignmentView_add(av, name_string, pad_string, middle_string,
                              pad_string, gap_string, query_pos, target_pos);
        else
            AlignmentView_add(av, name_string, NULL, middle_string, NULL,
                              gap_string, query_pos, target_pos);
        }
    g_free(label);
    g_free(name_string);
    g_free(gap_string);
    g_free(pad_string);
    g_free(middle_string);
    return;
    }

static void AlignmentView_add_NER(AlignmentView *av,
                          gint advance_query, gint advance_target,
                          gint query_pos, gint target_pos){
    register gchar
        *upper_string  = g_strdup_printf("%d", advance_query),
        *middle_string = g_strdup_printf("NER %d", ++av->ner_count),
        *lower_string  = g_strdup_printf("%d", advance_target);
    register gint upper_len = strlen(upper_string),
                  middle_len = strlen(middle_string),
                  lower_len = strlen(lower_string),
                  max_len;
    register gchar *upper_padded, *middle_padded, *lower_padded;
    max_len = upper_len;
    if(max_len < middle_len)
        max_len = middle_len;
    if(max_len < lower_len)
        max_len = lower_len;
    upper_padded = g_strdup_printf("--<%*c%s%*c>--",
            1+(((max_len-upper_len)+1)>>1), ' ',
            upper_string,
            1+((max_len-upper_len)>>1), ' ');
    middle_padded = g_strdup_printf("--<%*c%s%*c>--",
            1+(((max_len-middle_len)+1)>>1), ' ',
            middle_string,
            1+((max_len-middle_len)>>1), ' ');
    lower_padded = g_strdup_printf("--<%*c%s%*c>--",
            1+(((max_len-lower_len)+1)>>1), ' ',
            lower_string,
            1+((max_len-lower_len)>>1), ' ');
    g_free(upper_string);
    g_free(middle_string);
    g_free(lower_string);
    AlignmentView_add(av, upper_padded, NULL, middle_padded, NULL,
                      lower_padded, query_pos, target_pos);
    g_free(upper_padded);
    g_free(middle_padded);
    g_free(lower_padded);
    return;
    }

#define AlignmentView_combine_advance(advance_query, advance_target) \
        (((advance_query) << 8) | (advance_target))

static void AlignmentView_add_SPLIT_CODON(AlignmentView *av,
            Sequence *query, Sequence *target,
            gint advance_query, gint advance_target,
            gint query_pos, gint target_pos,
            Submat *protein_submat, Translate *translate){
    register gchar *query_string = NULL, *target_string = NULL,
                   *match_string = NULL, *codon_match;
    register gchar qy_aa_symbol = '\0', tg_aa_symbol = '\0',
                  *qy_aa_name = NULL, *tg_aa_name = NULL,
                   *tr_codon_name;
    register gchar *inner_query_string = NULL,
                   *inner_target_string = NULL;
    register gint pos_to_start = -1, num_to_print;
    register gint qp0 = 0, qp1 = 0, qp2 = 0,
                  tp0 = 0, tp1 = 0, tp2 = 0;
    register gboolean query_is_dna, target_is_dna, before_intron;
    gchar qy_codon[4], tg_codon[4];
    register AlignmentSeparation *as
        = av->split_codon_separation_list->pdata
         [av->curr_split_codon_count >> 1];
    g_assert(advance_query || advance_target);
    query_is_dna = (query->alphabet->type == Alphabet_Type_DNA);
    target_is_dna = (target->alphabet->type == Alphabet_Type_DNA);
    before_intron = (av->curr_split_codon_count & 1)?FALSE:TRUE;
    qy_codon[0] = qy_codon[1] = qy_codon[2] = qy_codon[3] = '\0';
    tg_codon[0] = tg_codon[1] = tg_codon[2] = tg_codon[3] = '\0';
    /**/
    if(query_is_dna && target_is_dna){
        switch(AlignmentView_combine_advance(
                    advance_query, advance_target)){
            case AlignmentView_combine_advance(1, 1):
                if(before_intron){
                    pos_to_start = 0;
                    qp0 = query_pos;
                    qp1 = query_pos + as->query_separation;
                    qp2 = query_pos + as->query_separation + 1;
                    tp0 = target_pos;
                    tp1 = target_pos + as->target_separation;
                    tp2 = target_pos + as->target_separation + 1;
                } else { /* after_intron */
                    pos_to_start = 2;
                    qp0 = query_pos - as->query_separation;
                    qp1 = query_pos - as->query_separation + 1;
                    qp2 = query_pos;
                    tp0 = target_pos - as->target_separation;
                    tp1 = target_pos - as->target_separation + 1;
                    tp2 = target_pos;
                    }
                break;
            case AlignmentView_combine_advance(2, 2):
                if(before_intron){
                    pos_to_start = 0;
                    qp0 = query_pos;
                    qp1 = query_pos + 1;
                    qp2 = query_pos + as->query_separation;
                    tp0 = target_pos;
                    tp1 = target_pos + 1;
                    tp2 = target_pos + as->target_separation;
                } else { /* after_intron */
                    pos_to_start = 1;
                    qp0 = query_pos - as->query_separation;
                    qp1 = query_pos;
                    qp2 = query_pos + 1;
                    tp0 = target_pos - as->target_separation;
                    tp1 = target_pos;
                    tp2 = target_pos + 1;
                    }
                break;
            default:
                g_error("Unexpected d2d split codon [%d,%d]",
                        advance_query, advance_target);
                break;
            }
        qy_codon[0] = Sequence_get_symbol(query, qp0);
        qy_codon[1] = Sequence_get_symbol(query, qp1);
        qy_codon[2] = Sequence_get_symbol(query, qp2);
        tg_codon[0] = Sequence_get_symbol(target, tp0);
        tg_codon[1] = Sequence_get_symbol(target, tp1);
        tg_codon[2] = Sequence_get_symbol(target, tp2);
    } else {
        if(query_is_dna){
            g_assert(target->alphabet->type == Alphabet_Type_PROTEIN);
            tg_aa_symbol = Sequence_get_symbol(target, target_pos);
            switch(AlignmentView_combine_advance(
                        advance_query, advance_target)){
                case AlignmentView_combine_advance(1, 0):
                    pos_to_start = 0;
                    qp0 = query_pos;
                    qp1 = query_pos + as->query_separation;
                    qp2 = query_pos + as->query_separation + 1;
                    break;
                case AlignmentView_combine_advance(2, 0):
                    pos_to_start = 0;
                    qp0 = query_pos;
                    qp1 = query_pos + 1;
                    qp2 = query_pos + as->query_separation;
                    break;
                case AlignmentView_combine_advance(2, 1):
                    pos_to_start = 1;
                    qp0 = query_pos - as->query_separation;
                    qp1 = query_pos;
                    qp2 = query_pos + 1;
                    break;
                case AlignmentView_combine_advance(1, 1):
                    pos_to_start = 2;
                    qp0 = query_pos - as->query_separation;
                    qp1 = query_pos - as->query_separation + 1;
                    qp2 = query_pos;
                    break;
                default:
                    g_error("Unexpected d2p split codon [%d,%d]",
                            advance_query, advance_target);
                    break;
                }
            qy_codon[0] = Sequence_get_symbol(query, qp0);
            qy_codon[1] = Sequence_get_symbol(query, qp1);
            qy_codon[2] = Sequence_get_symbol(query, qp2);
        } else {
            g_assert(query->alphabet->type == Alphabet_Type_PROTEIN);
            g_assert(target->alphabet->type == Alphabet_Type_DNA);
            qy_aa_symbol = Sequence_get_symbol(query, query_pos);
            switch(AlignmentView_combine_advance(
                        advance_query, advance_target)){
                case AlignmentView_combine_advance(0, 1):
                    pos_to_start = 0;
                    tp0 = target_pos;
                    tp1 = target_pos + as->target_separation;
                    tp2 = target_pos + as->target_separation + 1;
                    break;
                case AlignmentView_combine_advance(0, 2):
                    pos_to_start = 0;
                    tp0 = target_pos;
                    tp1 = target_pos + 1;
                    tp2 = target_pos + as->target_separation;
                    break;
                case AlignmentView_combine_advance(1, 2):
                    pos_to_start = 1;
                    tp0 = target_pos - as->target_separation;
                    tp1 = target_pos;
                    tp2 = target_pos + 1;
                    break;
                case AlignmentView_combine_advance(1, 1):
                    pos_to_start = 2;
                    tp0 = target_pos - as->target_separation;
                    tp1 = target_pos - as->target_separation + 1;
                    tp2 = target_pos;
                    break;
                default:
                    g_error("Unexpected p2d split codon [%d,%d]",
                            advance_query, advance_target);
                    break;
                }
            tg_codon[0] = Sequence_get_symbol(target, tp0);
            tg_codon[1] = Sequence_get_symbol(target, tp1);
            tg_codon[2] = Sequence_get_symbol(target, tp2);
            }
        }
    g_assert(pos_to_start != -1);
    av->curr_split_codon_count++;
/**/
    if(!query_is_dna)
        qy_aa_name = Alphabet_aa2tla(qy_aa_symbol);
    if(!target_is_dna)
        tg_aa_name = Alphabet_aa2tla(tg_aa_symbol);
    num_to_print = MAX(advance_query, advance_target);
    query_string = g_strdup_printf("{%.*s}",
                    num_to_print,
                    (query_is_dna?qy_codon:qy_aa_name)+pos_to_start);
    target_string = g_strdup_printf("{%.*s}",
                     num_to_print,
                     (target_is_dna?tg_codon:tg_aa_name)+pos_to_start);
    strup(qy_codon);
    strup(tg_codon);
    if(query_is_dna){
        g_assert(!qy_aa_symbol);
        qy_aa_symbol = Translate_codon(translate, qy_codon);
        tr_codon_name = Alphabet_aa2tla(qy_aa_symbol);
        inner_query_string = g_strdup_printf("{%.*s}",
                  num_to_print, tr_codon_name+pos_to_start);
        }
    if(target_is_dna){
        g_assert(!tg_aa_symbol);
        tg_aa_symbol = Translate_codon(translate, tg_codon);
        tr_codon_name = Alphabet_aa2tla(tg_aa_symbol);
        inner_target_string = g_strdup_printf("{%.*s}",
                      num_to_print, tr_codon_name+pos_to_start);
        }
    /**/
    if(query_is_dna){
        if(target_is_dna){ /* d,d */
            codon_match = g_strnfill(3, Alignment_get_equiv_symbol(
                            qy_aa_symbol, tg_aa_symbol,
                            protein_submat));
        } else { /* d,p */
            codon_match = Alignment_get_codon_match_string(
                            qy_codon, tg_aa_symbol,
                            protein_submat, translate);
            }
    } else { /* p,d */
        codon_match = Alignment_get_codon_match_string(
                        tg_codon, qy_aa_symbol,
                        protein_submat, translate);
        }
    match_string = g_strdup_printf("{%.*s}",
                         num_to_print, codon_match+pos_to_start);
    g_free(codon_match);
    /**/
    g_assert(query_string);
    g_assert(match_string);
    g_assert(target_string);
    AlignmentView_add(av, query_string, inner_query_string,
                      match_string, inner_target_string,
                      target_string, query_pos, target_pos);
    g_free(query_string);
    g_free(match_string);
    g_free(target_string);
    if(inner_query_string)
        g_free(inner_query_string);
    if(inner_target_string)
        g_free(inner_target_string);
    return;
    }

static void AlignmentView_add_FRAMESHIFT(AlignmentView *av,
                gint advance_query, gint advance_target,
                gint total_length, Sequence *query, Sequence *target,
                gint query_pos, gint target_pos, Translate *translate){
    register gint i, j;
    register gint curr_query_pos = query_pos,
                  curr_target_pos = target_pos;
    static gchar seq_string[4], match_string[4], gap_string[4];
    register Alphabet_Type emitted_alphabet_type;
    g_assert(!(advance_query && advance_target));
    match_string[0] = '#'; /* Frameshift */
    gap_string[0] = '-';
    if(advance_query)
        emitted_alphabet_type = query->alphabet->type;
    else /* advance_target */
        emitted_alphabet_type = target->alphabet->type;
    for(i = 0; i < total_length; i++){
        for(j = 0; j < (advance_query|advance_target); j++){
            if(advance_query)
                seq_string[j] = Sequence_get_symbol(query, curr_query_pos+j);
            else /* advance_target */
                seq_string[j] = Sequence_get_symbol(target, curr_target_pos+j);
            match_string[j] = match_string[0];
            gap_string[j] = gap_string[0];
            }
        seq_string[j] = match_string[j] = gap_string[j] = '\0';
        if(emitted_alphabet_type == Alphabet_Type_PROTEIN){
            strncpy(seq_string, Alphabet_aa2tla(seq_string[0]), 3);
            match_string[2] = match_string[1] = match_string[0];
            gap_string[2] = gap_string[1] = gap_string[0];
            seq_string[3] = match_string[3] = gap_string[3] = '\0';
            }
        if(advance_query)
            AlignmentView_add(av, seq_string,
                                  match_string,
                                  match_string,
                                  gap_string,
                                  gap_string,
                                  curr_query_pos, curr_target_pos);
        else
            AlignmentView_add(av, gap_string,
                                  gap_string,
                                  match_string,
                                  match_string,
                                  seq_string,
                                  curr_query_pos, curr_target_pos);
        curr_query_pos += advance_query;
        curr_target_pos += advance_target;
        }
    return;
    }

static void AlignmentView_add_label_operation(AlignmentView *av,
                C4_Transition *transition,
                gint total_length, Sequence *query, Sequence *target,
                gint query_pos, gint target_pos,
                Submat *dna_submat, Submat *protein_submat,
                Translate *translate,
                gboolean next_has_same_label, C4_Transition **last_match){
    switch(transition->label){
        case C4_Label_NONE:
            g_assert(!transition->advance_query);
            g_assert(!transition->advance_target);
            break;
        case C4_Label_MATCH:
            (*last_match) = transition;
            AlignmentView_add_MATCH(av, transition,
                    total_length, query, target,
                    query_pos, target_pos,
                    dna_submat, protein_submat, translate);
            break;
        case C4_Label_GAP:
            AlignmentView_add_GAP(av,
                    transition->advance_query,
                    transition->advance_target,
                    total_length, query, target,
                    query_pos, target_pos, translate);
            break;
        case C4_Label_5SS: /*fallthrough*/
            AlignmentView_add_SPLICE_SITE(av,
                    transition->advance_query,
                    transition->advance_target,
                    total_length, query, target,
                    query_pos, target_pos, TRUE, (*last_match));
            break;
        case C4_Label_3SS:
            AlignmentView_add_SPLICE_SITE(av,
                    transition->advance_query,
                    transition->advance_target,
                    total_length, query, target,
                    query_pos, target_pos, FALSE, (*last_match));
            break;
        case C4_Label_INTRON:
            av->intron_advance_query += (transition->advance_query
                                         * total_length);
            av->intron_advance_target += (transition->advance_target
                                         * total_length);
            if(!next_has_same_label){
                AlignmentView_add_INTRON(av,
                        av->intron_advance_query,
                        av->intron_advance_target,
                        query, target,
                        query_pos, target_pos, (*last_match));
                av->intron_advance_query = 0;
                av->intron_advance_target = 0;
                }
            break;
        case C4_Label_NER:
            av->ner_advance_query += (transition->advance_query
                                     * total_length);
            av->ner_advance_target += (transition->advance_target
                                     * total_length);
            if(!next_has_same_label){
                AlignmentView_add_NER(av, av->ner_advance_query,
                                          av->ner_advance_target,
                                          query_pos, target_pos);
                av->ner_advance_query = 0;
                av->ner_advance_target = 0;
                }
            break;
        case C4_Label_SPLIT_CODON:
            g_assert(total_length == 1);
            AlignmentView_add_SPLIT_CODON(av, query, target,
                              transition->advance_query,
                              transition->advance_target,
                              query_pos, target_pos,
                              protein_submat, translate);
            break;
        case C4_Label_FRAMESHIFT:
            AlignmentView_add_FRAMESHIFT(av,
                    transition->advance_query,
                    transition->advance_target,
                    total_length, query, target,
                    query_pos, target_pos, translate);
            break;
        default:
            g_error("AlignmentView cannot use label [%s]",
                    C4_Label_get_name(transition->label));
            break;
        }
    return;
    }

static void AlignmentView_prepare(AlignmentView *av,
                       Alignment *alignment,
                       Sequence *query, Sequence *target,
                       Submat *dna_submat, Submat *protein_submat,
                       Translate *translate){
    register gint i;
    register gint query_pos = alignment->region->query_start,
                  target_pos = alignment->region->target_start;
    register AlignmentPosition *apos;
    register AlignmentOperation *ao, *prev_ao;
    register gint total_length;
    C4_Transition *last_match = NULL;
    apos = g_new(AlignmentPosition, 1);
    apos->query_pos  = alignment->region->query_start-1;
    apos->target_pos = alignment->region->target_start-1;
    g_ptr_array_add(av->row_marker, apos);
    prev_ao = alignment->operation_list->pdata[0];
    total_length = prev_ao->length;
    for(i = 1; i < alignment->operation_list->len; i++){
        ao = alignment->operation_list->pdata[i];
        if(prev_ao->transition == ao->transition){
            /* group */
            total_length += ao->length;
        } else {
            /* report */
            AlignmentView_add_label_operation(av, prev_ao->transition,
                total_length, query, target, query_pos, target_pos,
                dna_submat, protein_submat, translate,
                (prev_ao->transition->label == ao->transition->label),
                &last_match);
            query_pos += (prev_ao->transition->advance_query
                        * total_length);
            target_pos += (prev_ao->transition->advance_target
                         * total_length);
            /* start new */
            prev_ao = ao;
            total_length = prev_ao->length;
            }
        }
    AlignmentView_add_label_operation(av, prev_ao->transition,
        total_length, query, target, query_pos, target_pos,
        dna_submat, protein_submat, translate, FALSE,
        &last_match);
    apos = g_new(AlignmentPosition, 1);
    apos->query_pos  = Region_query_end(alignment->region)-1;
    apos->target_pos = Region_target_end(alignment->region)-1;
    g_ptr_array_add(av->row_marker, apos);
    return;
    }

static gboolean AlignmentView_string_is_empty(GString *str, gint pos,
                                              gint width){
    register gint i;
    for(i = 0; i < width; i++)
        if(str->str[pos+i] != ' ')
            return FALSE;
    return TRUE;
    }

static void AlignmentView_prepare_seq(GString *outer, GString *inner,
                                      gint pos, gint width){
    register gint i;
    register gchar t;
    for(i = 0; i < width; i++){
        if(inner->str[pos+i] == ' '){ /* Swap empty */
             t = inner->str[pos+i];
             inner->str[pos+i] = outer->str[pos+i];
             outer->str[pos+i] = t;
             continue;
             }
        if(inner->str[pos+i] == '^') /* Replace padding */
            inner->str[pos+i] = ' ';
        }
    return;
    }

static void AlignmentView_replace_padding(GString *str, gint pos, gint width){
    register gint i;
    for(i = 0; i < width; i++){
        if(str->str[pos+i] == '^')
            str->str[pos+i] = ' ';
        }
    return;
    }

static void AlignmentView_display_row(AlignmentView *av,
            gint row, gint pos, gint width, gint maxposlen,
            Sequence *query, Sequence *target, FILE *fp){
    register AlignmentPosition *apos1 = av->row_marker->pdata[row],
                               *apos2 = av->row_marker->pdata[row+1];
    register gint p1q, p1t, p2q, p2t;
    register Alignment_ArgumentSet *aas
        = Alignment_ArgumentSet_create(NULL);
    register gboolean show_inner_query = FALSE,
                      show_inner_target = FALSE;
    p1q = apos1->query_pos+1;
    p2q = apos2->query_pos+1;
    p1t = apos1->target_pos+1;
    p2t = apos2->target_pos+1;
    if(aas->forward_strand_coords){
        if(query->strand == Sequence_Strand_REVCOMP){
            p1q = query->len - p1q - 1;
            p2q = query->len - p2q + 1;
            }
        if(target->strand == Sequence_Strand_REVCOMP){
            p1t = target->len - p1t - 1;
            p2t = target->len - p2t + 1;
            }
        }
    if(av->inner_query
    && (!AlignmentView_string_is_empty(av->inner_query, pos, width))){
        show_inner_query = TRUE;
        AlignmentView_prepare_seq(av->outer_query,
                                  av->inner_query, pos, width);
        }
    if(av->inner_target
    && (!AlignmentView_string_is_empty(av->inner_target, pos, width))){
        show_inner_target = TRUE;
        AlignmentView_prepare_seq(av->outer_target,
                                  av->inner_target, pos, width);
        }
    AlignmentView_replace_padding(av->outer_query, pos, width);
    AlignmentView_replace_padding(av->outer_target, pos, width);
    fprintf(fp, " %*d : %.*s : %*d\n", maxposlen, p1q+1,
            width, av->outer_query->str+pos,
            maxposlen, p2q);
    if(show_inner_query)
        fprintf(fp, " %*s   %.*s\n", maxposlen, " ", width,
                                 av->inner_query->str+pos);
    fprintf(fp, " %*s   %.*s\n", maxposlen, " ", width,
                             av->middle->str+pos);
    if(show_inner_target)
        fprintf(fp, " %*s   %.*s\n", maxposlen, " ", width,
                                 av->inner_target->str+pos);
    fprintf(fp, " %*d : %.*s : %*d\n", maxposlen, p1t+1,
            width, av->outer_target->str+pos,
            maxposlen, p2t);
    return;
    }

static void AlignmentView_display(AlignmentView *av,
                                  Sequence *query, Sequence *target,
                                  FILE *fp){
    register gint pos = 0, pause, row = 0;
    pause = av->outer_query->len-av->width;
    while(pos < pause){
        AlignmentView_display_row(av, row, pos, av->width,
                                  av->max_pos_len, query, target, fp);
        pos += av->width;
        row++;
        fprintf(fp, "\n");
        }
    AlignmentView_display_row(av, row, pos, av->outer_query->len-pos,
                              av->max_pos_len, query, target, fp);
    fprintf(fp, "\n");
    return;
    }

/**/

void Alignment_display(Alignment *alignment,
                       Sequence *query, Sequence *target,
                       Submat *dna_submat, Submat *protein_submat,
                       Translate *translate, FILE *fp){
    register AlignmentView *av = AlignmentView_create(alignment,
                                                      query, target);
    g_assert(alignment);
    g_assert(!alignment->model->is_open);
    /* Display header */
    fprintf(fp, "\n"
            "C4 Alignment:\n"
            "------------\n"
            "         Query: %s%s%s\n"
            "        Target: %s%s%s\n"
            "         Model: %s\n"
            "     Raw score: %d\n"
            "   Query range: %d -> %d\n"
            "  Target range: %d -> %d\n\n",
            query->id, query->def?" ":"", query->def?query->def:"",
            target->id, target->def?" ":"", target->def?target->def:"",
            alignment->model->name,
            alignment->score,
            Alignment_get_coordinate(alignment, query, target,
                                     TRUE, TRUE),
            Alignment_get_coordinate(alignment, query, target,
                                     TRUE, FALSE),
            Alignment_get_coordinate(alignment, query, target,
                                     FALSE, TRUE),
            Alignment_get_coordinate(alignment, query, target,
                                     FALSE, FALSE));
    AlignmentView_prepare(av, alignment, query, target,
                          dna_submat, protein_submat, translate);
    AlignmentView_display(av, query, target, fp);
    AlignmentView_destroy(av);
    return;
    }

/**/

static gint Alignment_get_equivalenced_matching(Alignment *alignment,
              Sequence *query, Sequence *target,
              Translate *translate, gboolean report_id,
              gpointer user_data){
    register gint i, j, match = 0;
    register AlignmentOperation *ao;
    register gint query_pos = alignment->region->query_start,
                  target_pos = alignment->region->target_start;
    register gchar qy_symbol, tg_symbol;
    for(i = 0; i < alignment->operation_list->len; i++){
        ao = alignment->operation_list->pdata[i];
        if(ao->transition->label == C4_Label_MATCH){
            for(j = 0; j < ao->length; j++){
                qy_symbol = Alignment_match_get_symbol(query, query_pos,
                                       ao->transition->advance_query,
                                       translate);
                tg_symbol = Alignment_match_get_symbol(target,
                                       target_pos,
                                       ao->transition->advance_target,
                                       translate);
                if(report_id){
                    if(toupper(qy_symbol) == toupper(tg_symbol))
                        match++;
                } else { /* report similarity */
                    if(C4_Calc_score(ao->transition->calc, query_pos,
                                  target_pos, user_data) > 0)
                        match++;
                    }
                query_pos += ao->transition->advance_query;
                target_pos += ao->transition->advance_target;
                }
        } else {
            query_pos += (ao->transition->advance_query
                          * ao->length);
            target_pos += (ao->transition->advance_target
                          * ao->length);
            }
        }
    return match;
    }
/* Reports the %id over the equivalenced regions of the alignment
 */

static gint Alignment_get_equivalenced_matching_region(Alignment *alignment,
              Sequence *query, Sequence *target,
              Translate *translate, gboolean report_id,
              gpointer user_data, gint exon_query_start, gint exon_query_end){
    register gint i, j, match = 0;
    register AlignmentOperation *ao;
    register gint query_pos = alignment->region->query_start,
                  target_pos = alignment->region->target_start;
    register gchar qy_symbol, tg_symbol;
    for(i = 0; i < alignment->operation_list->len; i++){
        ao = alignment->operation_list->pdata[i];
        if(ao->transition->label == C4_Label_MATCH){
            for(j = 0; j < ao->length; j++){
                if(query_pos > exon_query_end)
                    return match;
                if(query_pos >= exon_query_start){
                    qy_symbol = Alignment_match_get_symbol(query, query_pos,
                                           ao->transition->advance_query,
                                           translate);
                    tg_symbol = Alignment_match_get_symbol(target,
                                           target_pos,
                                           ao->transition->advance_target,
                                           translate);
                    if(report_id){
                        if(toupper(qy_symbol) == toupper(tg_symbol))
                            match++;
                    } else { /* report similarity */
                        if(C4_Calc_score(ao->transition->calc, query_pos,
                                      target_pos, user_data) > 0)
                            match++;
                        }
                    }
                query_pos += ao->transition->advance_query;
                target_pos += ao->transition->advance_target;
                }
        } else {
            query_pos += (ao->transition->advance_query
                          * ao->length);
            target_pos += (ao->transition->advance_target
                          * ao->length);
            }
        }
    return match;
    }
/* Reports the %id over the equivalenced regions of the alignment
 */

static gint Alignment_get_equivalenced_total(Alignment *alignment){
    register gint i, total = 0;
    register AlignmentOperation *ao;
    for(i = 0; i < alignment->operation_list->len; i++){
        ao = alignment->operation_list->pdata[i];
        if(ao->transition->label == C4_Label_MATCH)
            total += ao->length;
        }
    return total;
    }

/* report the number of gap bases found */
static int Alignment_get_gaps(Alignment *alignment){
    int i, total = 0;
    AlignmentOperation *ao;
    for(i = 0; i < alignment->operation_list->len; i++){
        ao = alignment->operation_list->pdata[i];
        if(ao->transition->label == C4_Label_GAP)
            total += ao->length;
        }
    return total;
    }

static gint Alignment_get_equivalenced_total_region(Alignment *alignment,
                          gint exon_query_start, gint exon_query_end){
    register gint i, j, total = 0;
    register gint query_pos = alignment->region->query_start,
                  target_pos = alignment->region->target_start;
    register AlignmentOperation *ao;
    for(i = 0; i < alignment->operation_list->len; i++){
        ao = alignment->operation_list->pdata[i];
        if(ao->transition->label == C4_Label_MATCH){
            for(j = 0; j < ao->length; j++){
                if(query_pos > exon_query_end)
                    return total;
                if(query_pos >= exon_query_start)
                    total++;
                query_pos += ao->transition->advance_query;
                target_pos += ao->transition->advance_target;
                }
        } else {
            query_pos += (ao->transition->advance_query
                          * ao->length);
            target_pos += (ao->transition->advance_target
                          * ao->length);
            }
        }
    return total;
    }
/* FIXME: optimisation: should use gene object to avoid full alignment pass
 */

static gfloat Alignment_get_percent_score_region(Alignment *alignment,
              Sequence *query, Sequence *target,
              Translate *translate, gboolean report_id, gpointer user_data,
              gint exon_query_start, gint exon_query_end){
    return (((gfloat)Alignment_get_equivalenced_matching_region(alignment,
                     query, target, translate, report_id, user_data,
                     exon_query_start, exon_query_end))
          / ((gfloat)Alignment_get_equivalenced_total_region(alignment,
                         exon_query_start, exon_query_end)))*100;
    }

static float Alignment_get_blast_percent_id(Alignment *alignment,
              Sequence *query, Sequence *target,
              Translate *translate, gpointer user_data){
    return ((float)Alignment_get_equivalenced_matching(alignment,
            query, target, translate, TRUE, user_data)
          / ((float)Alignment_get_equivalenced_total(alignment) +
             (float)Alignment_get_gaps(alignment)))*100;
    }

/* percent query coverage == equivalenced total / query length */
static float Alignment_get_percent_query_coverage(Alignment *alignment,
              Sequence *query){
    return ((float)Alignment_get_equivalenced_total(alignment)
          / (float)query->len)*100;
    }

/* FIXME: should also count split codons */

static gfloat Alignment_get_percent_score(Alignment *alignment,
              Sequence *query, Sequence *target,
              Translate *translate, gboolean report_id,
              gpointer user_data){
    return (((gfloat)Alignment_get_equivalenced_matching(alignment,
                     query, target, translate, report_id, user_data))
          / ((gfloat)Alignment_get_equivalenced_total(alignment)))*100;
    }
/* FIXME: should also count split codons */

static gint Alignment_get_match_score(Alignment *alignment,
                                      Sequence *query, Sequence *target,
                                      gpointer user_data){
    register gint i, j;
    register AlignmentOperation *ao;
    register C4_Score score = 0;
    register gint query_pos = alignment->region->query_start,
                  target_pos = alignment->region->target_start;
    for(i = 0;  i < alignment->operation_list->len; i++){
        ao = alignment->operation_list->pdata[i];
        for(j = 0; j < ao->length; j++){
            if(ao->transition->label == C4_Label_MATCH){
                score += C4_Calc_score(ao->transition->calc, query_pos,
                                       target_pos, user_data);
                }
            query_pos += ao->transition->advance_query;
            target_pos += ao->transition->advance_target;
            }
        }
    g_message("MATCH SCORE [%d]", score);
    return score;
    }

static gint Alignment_get_self_match_score(Alignment *alignment,
                                     Sequence *query, Sequence *target,
                                     gpointer self_data){
    register gint i, j;
    register AlignmentOperation *ao;
    register C4_Score score = 0;
    register gint query_pos = alignment->region->query_start;
    for(i = 0; i < alignment->operation_list->len; i++){
        ao = alignment->operation_list->pdata[i];
        for(j = 0; j < ao->length; j++){
            if(ao->transition->label == C4_Label_MATCH){
                /* Use self_data to find self-comparison scores */
                score += C4_Calc_score(ao->transition->calc, query_pos,
                                       query_pos, self_data);
                }
            query_pos += ao->transition->advance_query;
            }
        }
    g_message("SELF MATCH SCORE [%d]", score);
    return score;
    }

static gfloat Alignment_get_percent_self(Alignment *alignment,
              Sequence *query, Sequence *target,
              Translate *translate, gpointer user_data, gpointer self_data){
    return (((gfloat)Alignment_get_match_score(alignment,
                     query, target, user_data))
          / ((gfloat)Alignment_get_self_match_score(alignment, query, target,
                                                    self_data)))*100;
    }
/* Returns score as a percentage of maximum possible
 * over the equivalenced transitions
 */

/**/

static void Alignment_print_sugar_block(Alignment *alignment,
            Sequence *query, Sequence *target, FILE *fp){
    fprintf(fp, "%s %d %d %c %s %d %d %c %d",
             query->id,
             Alignment_get_coordinate(alignment, query, target,
                                      TRUE, TRUE),
             Alignment_get_coordinate(alignment, query, target,
                                      TRUE, FALSE),
             Sequence_get_strand_as_char(query),
             target->id,
             Alignment_get_coordinate(alignment, query, target,
                                      FALSE, TRUE),
             Alignment_get_coordinate(alignment, query, target,
                                      FALSE, FALSE),
             Sequence_get_strand_as_char(target),
             alignment->score);
    return;
    }

static gchar Alignment_get_cigar_type(AlignmentOperation *ao,
                                      gint *move){
    if(!ao->transition->advance_query){
        (*move) = ao->transition->advance_target * ao->length;
        return 'D';
        }
    if(!ao->transition->advance_target){
        (*move) = ao->transition->advance_query * ao->length;
        return 'I';
        }
    (*move) = MAX(ao->transition->advance_query,
                  ao->transition->advance_target) * ao->length;
    return 'M';
    }

static void Alignment_print_cigar_block(Alignment *alignment,
            Sequence *query, Sequence *target, FILE *fp){
    register gint i = 0;
    register gchar *gap = "";
    register AlignmentOperation *ao;
    register gchar type, next_type;
    gint move, next_move;
    ao = alignment->operation_list->pdata[i];
    type = Alignment_get_cigar_type(ao, &move);
    for(i = 1; i < alignment->operation_list->len; i++){
        ao = alignment->operation_list->pdata[i];
        next_type = Alignment_get_cigar_type(ao, &next_move);
        if(type == next_type){
            move += next_move;
        } else {
            if(move)
                fprintf(fp, "%s%c %d", gap, type, move);
            move = next_move;
            type = next_type;
            gap = " ";
            }
        }
    if(move)
        fprintf(fp, "%s%c %d", gap, type, move);
    return;
    }

static void Alignment_print_vulgar_block(Alignment *alignment,
            Sequence *query, Sequence *target, FILE *fp){
    register gint i;
    register gchar *gap = "";
    register C4_Label curr_label;
    register gint curr_advance_query, curr_advance_target;
    register AlignmentOperation *ao;
    register gboolean curr_is_codon = FALSE;
    ao = alignment->operation_list->pdata[0];
    curr_label = ao->transition->label;
    curr_advance_query = (ao->transition->advance_query * ao->length);
    curr_advance_target = (ao->transition->advance_target * ao->length);
    for(i = 1; i < alignment->operation_list->len; i++){
        ao = alignment->operation_list->pdata[i];
        if((ao->transition->label == curr_label)
         && (curr_advance_query || (!ao->transition->advance_query))
         && (curr_advance_target || (!ao->transition->advance_target))
         && (curr_is_codon == ((ao->transition->advance_query == 3)
                               && (ao->transition->advance_target == 3)))){
            curr_advance_query += (ao->transition->advance_query
                                  * ao->length);
            curr_advance_target += (ao->transition->advance_target
                                  * ao->length);
        } else {
            switch(curr_label){
                case C4_Label_NONE:
                    break;
                case C4_Label_MATCH:
                    g_assert(curr_advance_query && curr_advance_target);
                    fprintf(fp, "%s%c %d %d", gap, curr_is_codon?'C':'M',
                                              curr_advance_query,
                                              curr_advance_target);
                    gap = " ";
                    break;
                case C4_Label_GAP:
                    g_assert(curr_advance_query || curr_advance_target);
                    g_assert(!  (curr_advance_query
                              && curr_advance_target));
                    fprintf(fp, "%sG %d %d", gap, curr_advance_query,
                                                  curr_advance_target);
                    gap = " ";
                    break;
                case C4_Label_NER:
                    fprintf(fp, "%sN %d %d", gap, curr_advance_query,
                                                  curr_advance_target);
                    gap = " ";
                    break;
                case C4_Label_5SS:
                    fprintf(fp, "%s5 %d %d", gap, curr_advance_query,
                                                  curr_advance_target);
                    gap = " ";
                    break;
                case C4_Label_3SS:
                    fprintf(fp, "%s3 %d %d", gap, curr_advance_query,
                                                  curr_advance_target);
                    gap = " ";
                    break;
                case C4_Label_INTRON:
                    fprintf(fp, "%sI %d %d", gap, curr_advance_query,
                                                  curr_advance_target);
                    gap = " ";
                    break;
                case C4_Label_SPLIT_CODON:
                    fprintf(fp, "%sS %d %d", gap, curr_advance_query,
                                                  curr_advance_target);
                    gap = " ";
                    break;
                case C4_Label_FRAMESHIFT:
                    fprintf(fp, "%sF %d %d", gap, curr_advance_query,
                                                   curr_advance_target);
                    gap = " ";
                    break;
                default:
                    g_error("Unknown C4_Label [%d]", curr_label);
                    break;
                }
            curr_label = ao->transition->label;
            curr_is_codon = ((ao->transition->advance_query == 3)
                          && (ao->transition->advance_target == 3));
            curr_advance_query = (ao->transition->advance_query
                                  * ao->length);
            curr_advance_target = (ao->transition->advance_target
                                  * ao->length);
            }
        }
    return;
    }

typedef enum {
    Alignment_RYO_TOKEN_STRING,
    /**/
    Alignment_RYO_TOKEN_ID,           /* [QT] */
    Alignment_RYO_TOKEN_DEF,          /* [QT] */
    Alignment_RYO_TOKEN_LEN,          /* [QT] */
    Alignment_RYO_TOKEN_STRAND,       /* [QT] */
    Alignment_RYO_TOKEN_SEQ,          /* [QT] */
    Alignment_RYO_TOKEN_TYPE,         /* [QT] */
    /**/
    Alignment_RYO_TOKEN_ALIGN_BEGIN,  /* [QT] */
    Alignment_RYO_TOKEN_ALIGN_END,    /* [QT] */
    Alignment_RYO_TOKEN_ALIGN_LEN,    /* [QT] */
    Alignment_RYO_TOKEN_ALIGN_SEQ,    /* [QT] */
    /**/
    Alignment_RYO_TOKEN_CODING_BEGIN, /* [QT] */
    Alignment_RYO_TOKEN_CODING_END,   /* [QT] */
    Alignment_RYO_TOKEN_CODING_LEN,   /* [QT] */
    Alignment_RYO_TOKEN_CODING_SEQ,   /* [QT] */
    /**/
    Alignment_RYO_TOKEN_SCORE,
    Alignment_RYO_TOKEN_MODEL_NAME,
    Alignment_RYO_TOKEN_RANK,
    Alignment_RYO_TOKEN_BLAST_PERCENT_ID,
    Alignment_RYO_TOKEN_PERCENT_ID,
    Alignment_RYO_TOKEN_PERCENT_QUERY_COVERAGE,
    Alignment_RYO_TOKEN_PERCENT_SIMILARITY,
    Alignment_RYO_TOKEN_PERCENT_SELF,
    Alignment_RYO_TOKEN_GENE_ORIENTATION,
    Alignment_RYO_TOKEN_EQUIVALENCED_TOTAL,
    Alignment_RYO_TOKEN_EQUIVALENCED_IDENTCAL,
    Alignment_RYO_TOKEN_EQUIVALENCED_SIMILAR,
    Alignment_RYO_TOKEN_EQUIVALENCED_MISMATCHES,
    /**/
    Alignment_RYO_TOKEN_SUGAR_BLOCK,
    Alignment_RYO_TOKEN_CIGAR_BLOCK,
    Alignment_RYO_TOKEN_VULGAR_BLOCK,
    /**/
    Alignment_RYO_TOKEN_PTO_OPEN,
    Alignment_RYO_TOKEN_PTO_CLOSE,
    /**/
    Alignment_RYO_TOKEN_PTO_SEQ,     /* [QT] */
    Alignment_RYO_TOKEN_PTO_ADVANCE, /* [QT] */
    Alignment_RYO_TOKEN_PTO_BEGIN,   /* [QT] */
    Alignment_RYO_TOKEN_PTO_END,     /* [QT] */
    /**/
    Alignment_RYO_TOKEN_PTO_NAME,
    Alignment_RYO_TOKEN_PTO_SCORE,
    Alignment_RYO_TOKEN_PTO_LABEL
} Alignment_RYO_Token;

typedef struct {
    Alignment_RYO_Token  token;
                GString *str;      /* For TOKEN_STRING */
               gboolean  on_query; /* For [QT] */
} Alignment_RYO_ComplexToken;

static Alignment_RYO_ComplexToken *Alignment_RYO_ComplexToken_create(
       Alignment_RYO_Token token, gchar *str, gboolean on_query){
    register Alignment_RYO_ComplexToken *rct
     = g_new(Alignment_RYO_ComplexToken, 1);
    rct->token = token;
    if(str)
        rct->str = g_string_new(str);
    else
        rct->str = NULL;
    rct->on_query = on_query;
    return rct;
    }

static void Alignment_RYO_ComplexToken_destroy(
            Alignment_RYO_ComplexToken *rct){
    if(rct->str)
        g_string_free(rct->str, TRUE);
    g_free(rct);
    return;
    }

static void Alignment_RYO_add_string(GPtrArray *token_list, gchar *str){
    register Alignment_RYO_ComplexToken *rct = NULL;
    if(token_list->len){
        rct = token_list->pdata[token_list->len-1];
        if(rct->token != Alignment_RYO_TOKEN_STRING)
            rct = NULL;
        }
    if(rct){
        g_assert(rct->str);
        g_string_append(rct->str, str);
    } else {
        rct = Alignment_RYO_ComplexToken_create(
              Alignment_RYO_TOKEN_STRING, str, FALSE);
        g_ptr_array_add(token_list, rct);
        }
    return;
    }

static void Alignment_RYO_add_token(GPtrArray *token_list,
                                    Alignment_RYO_Token token,
                                    gboolean on_query){
    register Alignment_RYO_ComplexToken *rct
     = Alignment_RYO_ComplexToken_create(token, NULL, on_query);
    g_ptr_array_add(token_list, rct);
    return;
    }

static gint Alignment_RYO_tokenise_strand(GPtrArray *token_list,
                           gchar *format, gint pos, gboolean on_query){
    register gint move = 1;
    switch(format[pos]){
        case 'i':
            Alignment_RYO_add_token(token_list,
                      Alignment_RYO_TOKEN_ID, on_query);
            break;
        case 'd':
            Alignment_RYO_add_token(token_list,
                      Alignment_RYO_TOKEN_DEF, on_query);
            break;
        case 'l':
            Alignment_RYO_add_token(token_list,
                      Alignment_RYO_TOKEN_LEN, on_query);
            break;
        case 's':
            Alignment_RYO_add_token(token_list,
                      Alignment_RYO_TOKEN_SEQ, on_query);
            break;
        case 'S':
            Alignment_RYO_add_token(token_list,
                      Alignment_RYO_TOKEN_STRAND, on_query);
            break;
        case 't':
            Alignment_RYO_add_token(token_list,
                      Alignment_RYO_TOKEN_TYPE, on_query);
            break;
        case 'a':
            switch(format[pos+1]){
                case 'b':
                    Alignment_RYO_add_token(token_list,
                        Alignment_RYO_TOKEN_ALIGN_BEGIN, on_query);
                    break;
                case 'e':
                    Alignment_RYO_add_token(token_list,
                        Alignment_RYO_TOKEN_ALIGN_END, on_query);
                    break;
                case 'l':
                    Alignment_RYO_add_token(token_list,
                        Alignment_RYO_TOKEN_ALIGN_LEN, on_query);
                    break;
                case 's':
                    Alignment_RYO_add_token(token_list,
                        Alignment_RYO_TOKEN_ALIGN_SEQ, on_query);
                    break;
                default:
                    g_error("Unknown [%%%c%c%c] in format string [%s]",
                        format[pos-1], format[pos], format[pos+1],
                        format);
                    break;
                }
            move++;
            break;
        case 'c':
            switch(format[pos+1]){
                case 'b':
                    Alignment_RYO_add_token(token_list,
                        Alignment_RYO_TOKEN_CODING_BEGIN, on_query);
                    break;
                case 'e':
                    Alignment_RYO_add_token(token_list,
                        Alignment_RYO_TOKEN_CODING_END, on_query);
                    break;
                case 'l':
                    Alignment_RYO_add_token(token_list,
                        Alignment_RYO_TOKEN_CODING_LEN, on_query);
                    break;
                case 's':
                    Alignment_RYO_add_token(token_list,
                        Alignment_RYO_TOKEN_CODING_SEQ, on_query);
                    break;

                default:
                    g_error("Unknown [%%%c%c%c] in format string [%s]",
                        format[pos-1], format[pos], format[pos+1],
                        format);
                    break;
                }
            move++;
            break;
        default:
            g_error("Unknown [%%%c%c] in format string [%s]",
                    format[pos-1], format[pos], format);
            break;
        }
    return move;
    }

static gint Alignment_RYO_tokenise_PTO_strand(GPtrArray *token_list,
                           gchar *format, gint pos, gboolean on_query){
    register gint move = 1;
    switch(format[pos]){
        case 's':
            Alignment_RYO_add_token(token_list,
                        Alignment_RYO_TOKEN_PTO_SEQ, on_query);
            break;
        case 'a':
            Alignment_RYO_add_token(token_list,
                        Alignment_RYO_TOKEN_PTO_ADVANCE, on_query);
            break;
        case 'b':
            Alignment_RYO_add_token(token_list,
                        Alignment_RYO_TOKEN_PTO_BEGIN, on_query);
            break;
        case 'e':
            Alignment_RYO_add_token(token_list,
                        Alignment_RYO_TOKEN_PTO_END, on_query);
            break;
        default:
            g_error("Unknown [%%%c%c%c] in format string [%s]",
                    format[pos-2], format[pos-1], format[pos], format);
            break;
        }
    return move;
    }

static gint Alignment_RYO_tokenise_PTO(GPtrArray *token_list,
                           gchar *format, gint pos){
    register gint move = 1;
    switch(format[pos]){
        case 'q':
            move += Alignment_RYO_tokenise_PTO_strand(token_list,
                        format, pos+1, TRUE);
            break;
        case 't':
            move += Alignment_RYO_tokenise_PTO_strand(token_list,
                        format, pos+1, FALSE);
            break;
        case 'n':
            Alignment_RYO_add_token(token_list,
                        Alignment_RYO_TOKEN_PTO_NAME, FALSE);
            break;
        case 's':
            Alignment_RYO_add_token(token_list,
                        Alignment_RYO_TOKEN_PTO_SCORE, FALSE);
            break;
        case 'l':
            Alignment_RYO_add_token(token_list,
                        Alignment_RYO_TOKEN_PTO_LABEL, FALSE);
            break;
        default:
            g_error("Unknown [%%%c%c] in format string [%s]",
                    format[pos-1], format[pos], format);
            break;
        }
    return move;
    }

static GPtrArray *Alignment_RYO_tokenise(gchar *format){
    register gint i;
    register GPtrArray *token_list = g_ptr_array_new();
    gchar mini_str[2];
    mini_str[1] = '\0';
    for(i = 0; format[i]; i++){
        switch(format[i]){
            case '\\':
                switch(format[i+1]){
                    case '\\':
                        Alignment_RYO_add_string(token_list, "\\");
                        break;
                    case 'n':
                        Alignment_RYO_add_string(token_list, "\n");
                        break;
                    case 't':
                        Alignment_RYO_add_string(token_list, "\t");
                        break;
                    case '{':
                        Alignment_RYO_add_string(token_list, "{");
                        break;
                    case '}':
                        Alignment_RYO_add_string(token_list, "}");
                        break;
                    default:
                        g_error("Unknown [\\%c] in ryo string [%s]",
                                format[i+1], format);
                        break;
                    }
                i++;
                break;
            case '%':
                switch(format[i+1]){
                    case '%':
                        Alignment_RYO_add_string(token_list, "%");
                        break;
                    case 'q':
                        i += Alignment_RYO_tokenise_strand(token_list,
                                format, i+2, TRUE);
                        break;
                    case 't':
                        i += Alignment_RYO_tokenise_strand(token_list,
                                format, i+2, FALSE);
                        break;
                    case 's':
                        Alignment_RYO_add_token(token_list,
                            Alignment_RYO_TOKEN_SCORE, FALSE);
                        break;
                    case 'm':
                        Alignment_RYO_add_token(token_list,
                            Alignment_RYO_TOKEN_MODEL_NAME, FALSE);
                        break;
                    case 'r':
                        Alignment_RYO_add_token(token_list,
                            Alignment_RYO_TOKEN_RANK, FALSE);
                        break;
                    case 'p':
                        switch(format[i+2]){
                            case 'c':
                                Alignment_RYO_add_token(token_list,
                                 Alignment_RYO_TOKEN_PERCENT_QUERY_COVERAGE, FALSE);
                                break;
                            case 'I':
                                Alignment_RYO_add_token(token_list,
                                 Alignment_RYO_TOKEN_BLAST_PERCENT_ID, FALSE);
                                break;
                            case 'i':
                                Alignment_RYO_add_token(token_list,
                                 Alignment_RYO_TOKEN_PERCENT_ID, FALSE);
                                break;
                            case 's':
                                Alignment_RYO_add_token(token_list,
                                 Alignment_RYO_TOKEN_PERCENT_SIMILARITY,
                                 FALSE);
                                break;
                            case 'S':
                                Alignment_RYO_add_token(token_list,
                                 Alignment_RYO_TOKEN_PERCENT_SELF,
                                 FALSE);
                                break;
                            default:
                                g_error(
                                   "Unknown [%%%c%c] in format string",
                                   format[i+1], format[i+2]);
                                break;
                            }
                        i++;
                        break;
                    case 'e':
                        switch(format[i+2]){
                            case 't':
                                Alignment_RYO_add_token(token_list,
                                 Alignment_RYO_TOKEN_EQUIVALENCED_TOTAL,
                                 FALSE);
                                break;
                            case 'i':
                                Alignment_RYO_add_token(token_list,
                             Alignment_RYO_TOKEN_EQUIVALENCED_IDENTCAL,
                                 FALSE);
                                break;
                            case 's':
                                Alignment_RYO_add_token(token_list,
                             Alignment_RYO_TOKEN_EQUIVALENCED_SIMILAR,
                                 FALSE);
                                break;
                            case 'm':
                                Alignment_RYO_add_token(token_list,
                            Alignment_RYO_TOKEN_EQUIVALENCED_MISMATCHES,
                                 FALSE);
                                break;
                            default:
                                g_error(
                                   "Unknown [%%%c%c] in format string",
                                   format[i+1], format[i+2]);
                                break;
                            }
                        i++;
                        break;
                    case 'g':
                        Alignment_RYO_add_token(token_list,
                          Alignment_RYO_TOKEN_GENE_ORIENTATION, FALSE);
                        break;
                    case 'S':
                        Alignment_RYO_add_token(token_list,
                          Alignment_RYO_TOKEN_SUGAR_BLOCK, FALSE);
                        break;
                    case 'C':
                        Alignment_RYO_add_token(token_list,
                          Alignment_RYO_TOKEN_CIGAR_BLOCK, FALSE);
                        break;
                    case 'V':
                        Alignment_RYO_add_token(token_list,
                          Alignment_RYO_TOKEN_VULGAR_BLOCK, FALSE);
                        break;
                    case 'P':
                        i += Alignment_RYO_tokenise_PTO(token_list,
                                format, i+2);
                        break;
                    default:
                        g_error("Unknown [%%%c] in format string [%s]",
                                format[i+1], format);
                        break;
                    }
                i++;
                break;
            case '{':
                Alignment_RYO_add_token(token_list,
                  Alignment_RYO_TOKEN_PTO_OPEN, FALSE);
                break;
            case '}':
                Alignment_RYO_add_token(token_list,
                  Alignment_RYO_TOKEN_PTO_CLOSE, FALSE);
                break;
            default:
                mini_str[0] = format[i];
                Alignment_RYO_add_string(token_list, mini_str);
                break;
            }
        }
    return token_list;
    }
/*  %[qt][idlsSt] {query,target}{id,def,len,seq,Strand,type]
 *  %[qt]a[bels] {query,target}align{begin,end,len,seq}
 *  %s score
 *  %m model_name
 *  %r rank
 *  %p[isS] percent {id,similarity,self}
 *  %e[tid] equivalenced {total,identical,similar}
 *  %g gene_orientation
 *  %S sugar block
 *  %C cigar block
 *  %V vulgar block
 *  %% %
 *  \n newline
 *  \t tab
 *  \\ \
 *
 *  \{ {
 *  \} }
 *  { start per-transition section
 *  }   end per-transition section
 *  %P[qt][sabe] per-transition-output{query,target}
 *                                    {seq,advance,begin,end}
 *  %P[nsl] per-transition-output{name,score,label}
 */

static void Alignment_RYO_token_list_destroy(GPtrArray *token_list){
    register gint i;
    register Alignment_RYO_ComplexToken *rct;
    for(i = 0; i < token_list->len; i++){
        rct = token_list->pdata[i];
        Alignment_RYO_ComplexToken_destroy(rct);
        }
    g_ptr_array_free(token_list, TRUE);
    return;
    }

/**/

typedef struct {
             Alignment *alignment;
    AlignmentOperation *ao;
                  gint  operation_id;
                  gint  operation_pos;
                  gint  query_pos;
                  gint  target_pos;
              C4_Score *cell;
              C4_Score *start_cell;
              C4_Score  score;
              gpointer  user_data; /* No ref stored */
} Alignment_Position;

static Alignment_Position *Alignment_Position_create(
                    Alignment *alignment, gpointer user_data){
    register Alignment_Position *ap = g_new(Alignment_Position, 1);
    register gint i;
    register gint shadow_total
                    = 1 + alignment->model->total_shadow_designations;
    register C4_Calc *calc;
    ap->alignment = Alignment_share(alignment);
    ap->ao = alignment->operation_list->pdata[0];
    ap->operation_id = 0;
    ap->operation_pos = 0;
    ap->query_pos = alignment->region->query_start;
    ap->target_pos = alignment->region->target_start;
    ap->cell = g_new0(C4_Score, shadow_total);
    ap->user_data = user_data;
    if(alignment->model->start_state->cell_start_func){
        ap->start_cell = alignment->model->start_state->cell_start_func
            (ap->query_pos, ap->target_pos, user_data);
        for(i = 0; i < shadow_total; i++)
            ap->cell[i] = ap->start_cell[i];
        }
    if(alignment->model->init_func)
        alignment->model->init_func(alignment->region, user_data);
    for(i = 0; i < alignment->model->calc_list->len; i++){
        calc = alignment->model->calc_list->pdata[i];
        if(calc && calc->init_func)
            calc->init_func(alignment->region, user_data);
        }
    return ap;
    }

static void Alignment_Position_destroy(Alignment_Position *ap){
    register gint i;
    register C4_Calc *calc;
    for(i = 0; i < ap->alignment->model->calc_list->len; i++){
        calc = ap->alignment->model->calc_list->pdata[i];
        if(calc && calc->exit_func)
            calc->exit_func(ap->alignment->region, ap->user_data);
        }
    if(ap->alignment->model->exit_func)
        ap->alignment->model->exit_func(ap->alignment->region,
                                        ap->user_data);
    Alignment_destroy(ap->alignment);
    g_free(ap);
    return;
    }

static void Alignment_Position_set_shadows(Alignment_Position *ap){
    register gint i;
    register C4_State *state = ap->ao->transition->input;
    register C4_Shadow *shadow;
    for(i = 0; i < state->src_shadow_list->len; i++){
        shadow = state->src_shadow_list->pdata[i];
        ap->cell[1 + shadow->designation] = shadow->start_func(
                ap->query_pos, ap->target_pos,
                ap->user_data);
        }
    for(i = 0; i < ap->ao->transition->dst_shadow_list->len; i++){
        shadow = ap->ao->transition->dst_shadow_list->pdata[i];
        shadow->end_func(ap->cell[1 + shadow->designation],
                ap->query_pos + ap->ao->transition->advance_query,
                ap->target_pos + ap->ao->transition->advance_target,
                ap->user_data);
        }
    return;
    }

static gboolean Alignment_Position_next(Alignment_Position *ap){
    ap->query_pos += ap->ao->transition->advance_query;
    ap->target_pos += ap->ao->transition->advance_target;
    /* Use next transition */
    if(++ap->operation_pos < ap->ao->length){
        Alignment_Position_set_shadows(ap);
        return TRUE;
        }
    /* Use next operation */
    if(++ap->operation_id < ap->alignment->operation_list->len){
        ap->operation_pos = 0;
        ap->ao = ap->alignment->operation_list->pdata[ap->operation_id];
        Alignment_Position_set_shadows(ap);
        return TRUE;
        }
    return FALSE;
    }

/**/

typedef struct {
        gint  begin;
        gint  end;
    Sequence *seq;
} Alignment_Coding;

static Alignment_Coding *Alignment_Coding_create(Alignment *alignment,
                                                 Sequence *query,
                                                 Sequence *target,
                                                 gpointer user_data,
                                                 gboolean on_query){
    register Alignment_Coding *ac = g_new(Alignment_Coding, 1);
    register Alignment_Position *ap
        = Alignment_Position_create(alignment, user_data);
    register Sequence *src = NULL;
    register GString *seq = g_string_sized_new(1024);
    register gint i, advance, pos;
    if(on_query){
        g_assert(query->alphabet->type == Alphabet_Type_DNA);
        src = query;
    } else {
        g_assert(target->alphabet->type == Alphabet_Type_DNA);
        src = target;
        }
    while(Alignment_Position_next(ap)){
        if(on_query){
            advance = ap->ao->transition->advance_query;
            pos = ap->query_pos;
        } else {
            advance = ap->ao->transition->advance_target;
            pos = ap->target_pos;
            }
        switch(ap->ao->transition->label){
            case C4_Label_MATCH:
                if(advance == 3){
                    if(!seq->len)
                        ac->begin = pos;
                    g_string_append_c(seq, Sequence_get_symbol(src, pos));
                    g_string_append_c(seq, Sequence_get_symbol(src, pos+1));
                    g_string_append_c(seq, Sequence_get_symbol(src, pos+2));
                    ac->end = pos;
                    }
                break;
            case C4_Label_SPLIT_CODON:
                g_assert(advance);
                g_assert(seq);
                g_assert(seq->len);
                for(i = 0; i < advance; i++)
                    g_string_append_c(seq, Sequence_get_symbol(src, pos+i));
                break;
            case C4_Label_GAP:
                if(advance == 3){
                    g_assert(seq);
                    g_assert(seq->len);
                    g_string_append_c(seq, Sequence_get_symbol(src, pos));
                    g_string_append_c(seq, Sequence_get_symbol(src, pos+1));
                    g_string_append_c(seq, Sequence_get_symbol(src, pos+2));
                    }
                break;
            default: /* ignored */
                break;
            }
        }
    Alignment_Position_destroy(ap);
    ac->seq = Sequence_create("coding", NULL, seq->str, seq->len,
                              src->strand, src->alphabet);
    g_string_free(seq, TRUE);
    return ac;
    }

static void Alignment_Coding_destroy(Alignment_Coding *ac){
    Sequence_destroy(ac->seq);
    g_free(ac);
    return;
    }

/**/

static void Alignment_RYO_token_list_print(GPtrArray *token_list,
            Alignment *alignment, Sequence *query, Sequence *target,
            Translate *translate, gint rank,
            gpointer user_data, gpointer self_data, FILE *fp){
    register gint i, j, pto_start = -1;
    register Alignment_RYO_ComplexToken *rct;
    register Sequence *seq, *subseq;
    register Alignment_Position *ap = NULL;
    register Alignment_Coding *qy_ac = NULL, *tg_ac = NULL, *ac;
    for(i = 0; i < token_list->len; i++){
        rct = token_list->pdata[i];
        seq = rct->on_query?query:target;
        switch(rct->token){
            case Alignment_RYO_TOKEN_STRING:
                fprintf(fp, "%s", rct->str->str);
                break;
            /**/
            case Alignment_RYO_TOKEN_ID:
                fprintf(fp, "%s", seq->id);
                break;
            case Alignment_RYO_TOKEN_DEF:
                if(seq->def)
                    fprintf(fp, "%s", seq->def);
                break;
            case Alignment_RYO_TOKEN_LEN:
                fprintf(fp, "%d", seq->len);
                break;
            case Alignment_RYO_TOKEN_STRAND:
                fprintf(fp, "%c", Sequence_get_strand_as_char(seq));
                break;
            case Alignment_RYO_TOKEN_SEQ:
                Sequence_print_fasta_block(seq, fp);
                break;
            case Alignment_RYO_TOKEN_TYPE:
                fprintf(fp, "%s",
                        Alphabet_Type_get_name(seq->alphabet->type));
                break;
            /**/
            case Alignment_RYO_TOKEN_ALIGN_BEGIN:
                fprintf(fp, "%d",
                    Alignment_get_coordinate(alignment, query, target,
                        rct->on_query, TRUE));
                break;
            case Alignment_RYO_TOKEN_ALIGN_END:
                fprintf(fp, "%d",
                    Alignment_get_coordinate(alignment, query, target,
                        rct->on_query, FALSE));
                break;
            case Alignment_RYO_TOKEN_ALIGN_LEN:
                fprintf(fp, "%d",
                        rct->on_query?alignment->region->query_length
                                     :alignment->region->target_length);
                break;
            case Alignment_RYO_TOKEN_ALIGN_SEQ:
                if(rct->on_query){
                    subseq = Sequence_subseq(seq,
                                alignment->region->query_start,
                                alignment->region->query_length);
                } else {
                    subseq = Sequence_subseq(seq,
                                alignment->region->target_start,
                                alignment->region->target_length);
                    }
                Sequence_print_fasta_block(subseq, fp);
                Sequence_destroy(subseq);
                break;
            /**/
            case Alignment_RYO_TOKEN_CODING_BEGIN:
                ac = rct->on_query?qy_ac:tg_ac;
                if(!ac)
                    ac = Alignment_Coding_create(alignment,
                            query, target, user_data, rct->on_query);
                fprintf(fp, "%d", Alignment_convert_coordinate(alignment,
                              query, target,
                              ac->begin, ac->begin, rct->on_query));
                break;
            case Alignment_RYO_TOKEN_CODING_END:
                ac = rct->on_query?qy_ac:tg_ac;
                if(!ac)
                    ac = Alignment_Coding_create(alignment,
                            query, target, user_data, rct->on_query);
                fprintf(fp, "%d", Alignment_convert_coordinate(alignment,
                              query, target,
                              ac->end, ac->end, rct->on_query));
                break;
            case Alignment_RYO_TOKEN_CODING_LEN:
                ac = rct->on_query?qy_ac:tg_ac;
                if(!ac)
                    ac = Alignment_Coding_create(alignment,
                            query, target, user_data, rct->on_query);
                fprintf(fp, "%d", ac->seq->len);
                break;
            case Alignment_RYO_TOKEN_CODING_SEQ:
                ac = rct->on_query?qy_ac:tg_ac;
                if(!ac)
                    ac = Alignment_Coding_create(alignment,
                            query, target, user_data, rct->on_query);
                Sequence_print_fasta_block(ac->seq, fp);
                break;
            /**/
            case Alignment_RYO_TOKEN_SCORE:
                fprintf(fp, "%d", alignment->score);
                break;
            case Alignment_RYO_TOKEN_MODEL_NAME:
                fprintf(fp, "%s", alignment->model->name);
                break;
            case Alignment_RYO_TOKEN_RANK:
                if(rank == -1)
                    fprintf(fp, "%%_EXONERATE_BESTN_RANK_%%");
                else
                    fprintf(fp, "%d", rank);
                break;
            /**/
            case Alignment_RYO_TOKEN_PERCENT_QUERY_COVERAGE:
                fprintf(fp, "%2.2f", Alignment_get_percent_query_coverage(
                                 alignment, query));
                break;
            case Alignment_RYO_TOKEN_BLAST_PERCENT_ID:
                fprintf(fp, "%2.2f", Alignment_get_blast_percent_id(
                                 alignment, query, target,
                                 translate, user_data));
                break;
            case Alignment_RYO_TOKEN_PERCENT_ID:
                fprintf(fp, "%2.2f", Alignment_get_percent_score(
                                 alignment, query, target,
                                 translate, TRUE, user_data));
                break;
            case Alignment_RYO_TOKEN_PERCENT_SIMILARITY:
                fprintf(fp, "%2.2f", Alignment_get_percent_score(
                                 alignment, query, target,
                                 translate, FALSE, user_data));
                break;
            case Alignment_RYO_TOKEN_PERCENT_SELF:
                fprintf(fp, "%2.2f", Alignment_get_percent_self(
                                 alignment, query, target,
                                 translate, user_data, self_data));
                break;
            case Alignment_RYO_TOKEN_GENE_ORIENTATION:
                fprintf(fp, "%c",
                     Alignment_get_gene_orientation(alignment));
                break;
            /**/
            case Alignment_RYO_TOKEN_EQUIVALENCED_TOTAL:
                fprintf(fp, "%d", Alignment_get_equivalenced_total(
                                 alignment));
                break;
            case Alignment_RYO_TOKEN_EQUIVALENCED_IDENTCAL:
                fprintf(fp, "%d", Alignment_get_equivalenced_matching(
                                 alignment, query, target,
                                 translate, TRUE, user_data));
                break;
            case Alignment_RYO_TOKEN_EQUIVALENCED_SIMILAR:
                fprintf(fp, "%d", Alignment_get_equivalenced_matching(
                                 alignment, query, target,
                                 translate, FALSE, user_data));
                break;
            case Alignment_RYO_TOKEN_EQUIVALENCED_MISMATCHES:
                fprintf(fp, "%d",
                        Alignment_get_equivalenced_total(alignment)
                      - Alignment_get_equivalenced_matching(
                                 alignment, query, target,
                                 translate, TRUE, user_data));
                break;
            /**/
            case Alignment_RYO_TOKEN_SUGAR_BLOCK:
                Alignment_print_sugar_block(alignment, query, target, fp);
                break;
            case Alignment_RYO_TOKEN_CIGAR_BLOCK:
                Alignment_print_cigar_block(alignment, query, target, fp);
                break;
            case Alignment_RYO_TOKEN_VULGAR_BLOCK:
                Alignment_print_vulgar_block(alignment, query, target, fp);
                break;
            /**/
            case Alignment_RYO_TOKEN_PTO_OPEN:
                if(pto_start != -1)
                    g_error("Cannot nest PTO brackets");
                pto_start = i;
                ap = Alignment_Position_create(alignment, user_data);
                break;
            case Alignment_RYO_TOKEN_PTO_CLOSE:
                if(pto_start == -1)
                    g_error("No opening PTO bracket in --ryo string");
                if(Alignment_Position_next(ap)){
                    i = pto_start;
                } else {
                    pto_start = -1;
                    Alignment_Position_destroy(ap);
                    ap = NULL;
                    }
                break;
            /**/
            case Alignment_RYO_TOKEN_PTO_SEQ:
                g_assert(pto_start != -1);
                if(rct->on_query){
                    if(ap->ao->transition->advance_query){
                        for(j = 0; j < ap->ao->transition->advance_query; j++)
                            fprintf(fp, "%c",
                                Sequence_get_symbol(query, ap->query_pos+j));
                    } else {
                        fprintf(fp, "-");
                        }
                } else {
                    if(ap->ao->transition->advance_target){
                        for(j = 0; j < ap->ao->transition->advance_target; j++)
                            fprintf(fp, "%c",
                                Sequence_get_symbol(target, ap->target_pos+j));
                    } else {
                        fprintf(fp, "-");
                        }
                    }
                break;
            case Alignment_RYO_TOKEN_PTO_ADVANCE:
                fprintf(fp, "%d",
                   rct->on_query?ap->ao->transition->advance_query
                                :ap->ao->transition->advance_target);
                break;
            case Alignment_RYO_TOKEN_PTO_BEGIN:
                fprintf(fp, "%d", Alignment_convert_coordinate(alignment,
                              query, target,
                              ap->query_pos, ap->target_pos, rct->on_query));
                break;
            case Alignment_RYO_TOKEN_PTO_END:
                fprintf(fp, "%d", Alignment_convert_coordinate(alignment,
                              query, target,
                              ap->query_pos+ap->ao->transition->advance_query,
                              ap->target_pos+ap->ao->transition->advance_target,
                              rct->on_query));
                break;
            /**/
            case Alignment_RYO_TOKEN_PTO_NAME:
                fprintf(fp, "%s", ap->ao->transition->name);
                break;
            case Alignment_RYO_TOKEN_PTO_SCORE:
                fprintf(fp, "%d",
                        C4_Calc_score(ap->ao->transition->calc,
                                      ap->query_pos,
                                      ap->target_pos, user_data));
                break;
            case Alignment_RYO_TOKEN_PTO_LABEL:
                fprintf(fp, "%s",
                    C4_Label_get_name(ap->ao->transition->label));
                break;
            default:
                g_error("Unknown token [%d]", rct->token);
                break;
            }
        }
    if(pto_start != -1)
        g_error("No closing PTO bracket in --ryo string");
    if(qy_ac)
        Alignment_Coding_destroy(qy_ac);
    if(tg_ac)
        Alignment_Coding_destroy(tg_ac);
    return;
    }

void Alignment_display_ryo(Alignment *alignment,
        Sequence *query, Sequence *target, gchar *format,
        Translate *translate, gint rank,
        gpointer user_data, gpointer self_data, FILE *fp){
    register GPtrArray *token_list = Alignment_RYO_tokenise(format);
    Alignment_RYO_token_list_print(token_list, alignment,
                                   query, target, translate, rank,
                                   user_data, self_data, fp);
    Alignment_RYO_token_list_destroy(token_list);
    return;
    }

void Alignment_display_sugar(Alignment *alignment,
                             Sequence *query, Sequence *target, FILE *fp){
    fprintf(fp, "sugar: ");
    Alignment_print_sugar_block(alignment, query, target, fp);
    fprintf(fp, "\n");
    return;
    }
/* sugar: simple ungapped alignment report
 */

void Alignment_display_cigar(Alignment *alignment,
                             Sequence *query, Sequence *target, FILE *fp){
    fprintf(fp, "cigar: ");
    Alignment_print_sugar_block(alignment, query, target, fp);
    fprintf(fp, " ");
    Alignment_print_cigar_block(alignment, query, target, fp);
    fprintf(fp, "\n");
    return;
    }
/* cigar: concise idiosyncratic gapped alignment report
 * query_id query_start query_end query_strand \
 * target_id target_start target_end target_strand \
 * score \
 * <[MID] length> (m=Match i=Insert d=Delete)
 *
 */

void Alignment_display_vulgar(Alignment *alignment,
                              Sequence *query, Sequence *target, FILE *fp){
    fprintf(fp, "vulgar: ");
    Alignment_print_sugar_block(alignment, query, target, fp);
    fprintf(fp, " ");
    Alignment_print_vulgar_block(alignment, query, target, fp);
    fprintf(fp, "\n");
    return;
    }

/**/

static void Alignment_display_gff_header(Alignment *alignment,
        Sequence *query, Sequence *target, gboolean report_on_query, FILE *fp){
    time_t timenow = time(NULL);
    register struct tm *tm_ptr = localtime(&timenow);
    gchar curr_time_str[12];
    strftime(curr_time_str, 30, "%Y-%m-%d", tm_ptr);
    fprintf(fp, "#\n"
                "##gff-version 2\n"
                "##source-version %s:%s %s\n"
                "##date %s\n"
                "##type %s\n"
                "#\n",
                PACKAGE, alignment->model->name, VERSION,
                curr_time_str,
                Alphabet_Type_get_name(report_on_query
                                      ? query->alphabet->type
                                      : target->alphabet->type));
    fprintf(fp, "#\n# seqname source feature start end"
                " score strand frame attributes\n#\n");
    return;
    }

static void Alignment_display_gff_line(Alignment *alignment,
                                       Sequence *query,
                                       Sequence *target,
                                       gboolean report_on_query,
                                       gchar *feature,
                                       gint query_start,
                                       gint target_start,
                                       gint query_end,
                                       gint target_end,
                                       gboolean show_score, gint score,
                                       gboolean show_frame, gint frame,
                                       GPtrArray *attribute_list, FILE *fp){
    /* <seqname> <source> <feature> <start> <end> <score> <strand> <frame>
     * [attributes] [comments]
     */
    register Sequence *seq;
    register gint i, start, end, t;
    register gchar *attribute;
    if(report_on_query){
        seq = query;
        start = query_start;
        end = query_end;
    } else {
        seq = target;
        start = target_start;
        end = target_end;
        }
    if(seq->strand == Sequence_Strand_REVCOMP){
        start = seq->len - start;
        end = seq->len - end;
        /* swap coords for gff output */
        t = start;
        start = end;
        end = t;
        }
    /* Sanity checks for valid GFF output */
    g_assert(start >= 0);
    g_assert(end <= seq->len);
    g_assert(start < end);
    g_assert((!show_frame) || ((frame >= 0) && (frame <= 2)));
    fprintf(fp, "%s\t%s:%s\t%s\t%d\t%d\t",
            seq->id, PACKAGE, alignment->model->name, feature, start+1, end);
    if(show_score)
        fprintf(fp, "%d", score);
    else
        fprintf(fp, ".");
    fprintf(fp, "\t%c\t", Sequence_get_strand_as_char(seq));
    if(show_frame)
        fprintf(fp, "%d", frame);
    else
        fprintf(fp, ".");
    fprintf(fp, "\t");
    if(attribute_list){
        for(i = 0; i < attribute_list->len; i++){
            attribute = attribute_list->pdata[i];
            fprintf(fp, "%s", attribute);
            if((i+1) < attribute_list->len)
                fprintf(fp, " ; ");
            }
        }
    fprintf(fp, "\n");
    return;
    }

static void Alignment_free_attribute_list(GPtrArray *attribute_list){
    register gint i;
    for(i = 0; i < attribute_list->len; i++)
        g_free(attribute_list->pdata[i]);
    g_ptr_array_free(attribute_list, TRUE);
    return;
    }

static void Alignment_display_gff_exon(Alignment *alignment,
                                       Sequence *query,
                                       Sequence *target,
                                       Translate *translate,
                                       gboolean report_on_query,
                                       gint query_pos,
                                       gint target_pos,
                                       gint exon_query_start,
                                       gint exon_target_start,
                                       gint exon_query_gap,
                                       gint exon_target_gap,
                                       gint exon_query_frameshift,
                                       gint exon_target_frameshift,
                                       gpointer user_data, FILE *fp){
    register GPtrArray *attribute_list = g_ptr_array_new();
    g_ptr_array_add(attribute_list,
                    g_strdup_printf("insertions %d",
                                    report_on_query?exon_query_gap
                                                   :exon_target_gap));
    g_ptr_array_add(attribute_list,
                    g_strdup_printf("deletions %d",
                                    report_on_query?exon_target_gap
                                                   :exon_query_gap));
    g_ptr_array_add(attribute_list,
                    g_strdup_printf("identity %2.2f",
                                 Alignment_get_percent_score_region(alignment,
                                 query, target, translate, TRUE, user_data,
                                 exon_query_start, query_pos)));
    g_ptr_array_add(attribute_list,
                    g_strdup_printf("similarity %2.2f",
                                 Alignment_get_percent_score_region(alignment,
                                 query, target, translate, FALSE, user_data,
                                 exon_query_start, query_pos)));
    if(report_on_query){
        if(exon_query_frameshift)
            g_ptr_array_add(attribute_list,
                            g_strdup_printf("frameshifts %d",
                                            exon_query_frameshift));
    } else {
        if(exon_target_frameshift)
            g_ptr_array_add(attribute_list,
                            g_strdup_printf("frameshifts %d",
                                            exon_target_frameshift));
        }
    Alignment_display_gff_line(alignment, query, target,
                               report_on_query, "exon",
                               exon_query_start, exon_target_start,
                               query_pos, target_pos,
                               FALSE, 0, FALSE, 0, attribute_list, fp);
    Alignment_free_attribute_list(attribute_list);
    return;
    }

static void Alignment_display_gff_utr(Alignment *alignment,
                            Sequence *query, Sequence *target,
                            gboolean report_on_query, gboolean post_cds,
                            gint cds_query_start, gint cds_target_start,
                            gint cds_query_end, gint cds_target_end,
                            gint exon_query_start, gint exon_target_start,
                            gint query_pos, gint target_pos, FILE *fp){
    register gint curr_cds_query_start, curr_cds_target_start,
                  curr_utr_query_start, curr_utr_target_start;
    if(post_cds){
        curr_utr_query_start = MAX(exon_query_start, cds_query_end);
        curr_utr_target_start = MAX(exon_target_start, cds_target_end);
        Alignment_display_gff_line(alignment, query, target,
                               report_on_query, "utr3",
                               curr_utr_query_start, curr_utr_target_start,
                               query_pos, target_pos,
                               FALSE, 0, FALSE, 0, NULL, fp);
    } else if(cds_query_start == -1){
        Alignment_display_gff_line(alignment, query, target,
                               report_on_query, "utr5",
                               exon_query_start, exon_target_start,
                               query_pos, target_pos,
                               FALSE, 0, FALSE, 0, NULL, fp);
    } else {
        curr_cds_query_start = MAX(cds_query_start,
                                   exon_query_start);
        curr_cds_target_start = MAX(cds_target_start,
                                    exon_target_start);
        Alignment_display_gff_line(alignment, query, target,
                               report_on_query, "cds",
                               curr_cds_query_start, curr_cds_target_start,
                               query_pos, target_pos,
                               FALSE, 0, FALSE, 0, NULL, fp);
        }
    return;
    }

static void Alignment_display_gff_gene(Alignment *alignment,
        Sequence *query, Sequence *target, Translate *translate,
        gboolean report_on_query,
        gint result_id, gpointer user_data, FILE *fp){
    register gint i;
    register gint query_pos = alignment->region->query_start,
                  target_pos = alignment->region->target_start;
    register gint intron_id = 0, intron_length = 0;
    register gint exon_query_start = 0, exon_target_start = 0;
    register gint exon_query_gap = 0, exon_target_gap = 0;
    register gint exon_query_frameshift = 0, exon_target_frameshift = 0;
    register gint cds_query_start = -1, cds_target_start = -1,
                  cds_query_end = -1, cds_target_end = -1;
    register gboolean in_exon = FALSE, post_cds = FALSE;
    register AlignmentOperation *ao;
    register gchar gene_orientation
               = Alignment_get_gene_orientation(alignment);
    register GPtrArray *attribute_list = g_ptr_array_new();
    register gint curr_utr_query_start, curr_utr_target_start;
    /**/
    g_ptr_array_add(attribute_list,
                    g_strdup_printf("gene_id %d", result_id));
    g_ptr_array_add(attribute_list,
                    g_strdup_printf("sequence %s",
                        report_on_query?target->id:query->id));
    g_ptr_array_add(attribute_list,
                    g_strdup_printf("gene_orientation %c", gene_orientation));
    g_ptr_array_add(attribute_list,
                    g_strdup_printf("identity %2.2f",
                                 Alignment_get_percent_score(alignment,
                                 query, target, translate, TRUE, user_data)));
    g_ptr_array_add(attribute_list,
                    g_strdup_printf("similarity %2.2f",
                                 Alignment_get_percent_score(alignment,
                                 query, target, translate, FALSE, user_data)));
    Alignment_display_gff_line(alignment, query, target, report_on_query,
                               "gene",
                               alignment->region->query_start,
                               alignment->region->target_start,
                               Region_query_end(alignment->region),
                               Region_target_end(alignment->region),
                               TRUE, alignment->score,
                               FALSE, 0, attribute_list, fp);
    Alignment_free_attribute_list(attribute_list);
    for(i = 1; i < alignment->operation_list->len; i++){
        ao = alignment->operation_list->pdata[i];
        switch(ao->transition->label){
            case C4_Label_MATCH:
                if((ao->transition->advance_query == 1)
                && (ao->transition->advance_target == 1)){
                    if((cds_query_start != -1) && (!post_cds)){
                        Alignment_display_gff_line(alignment, query, target,
                               report_on_query, "cds",
                               exon_query_start, exon_target_start,
                               query_pos, target_pos,
                               FALSE, 0, FALSE, 0, NULL, fp);
                        post_cds = TRUE;
                        }
                } else {
                    if(cds_query_start == -1){ /* First coding exon */
                        if(in_exon){ /* Have UTR5 */
                            Alignment_display_gff_line(alignment, query, target,
                               report_on_query, "utr5",
                               exon_query_start, exon_target_start,
                               query_pos, target_pos,
                               FALSE, 0, FALSE, 0, NULL, fp);
                            }
                        cds_query_start = query_pos;
                        cds_target_start = target_pos;
                        }
                    cds_query_end = query_pos
                                  + (ao->transition->advance_query
                                   * ao->length);
                    cds_target_end = target_pos
                                   + (ao->transition->advance_target
                                    * ao->length);
                    }
                /*fallthrough*/
            case C4_Label_SPLIT_CODON:
                if(!in_exon){
                    exon_query_start = query_pos;
                    exon_target_start = target_pos;
                    exon_query_gap = 0;
                    exon_target_gap = 0;
                    exon_query_frameshift = 0;
                    exon_target_frameshift = 0;
                    in_exon = TRUE;
                    }
                break;
            case C4_Label_NONE:
                break;
            case C4_Label_GAP:
                exon_query_gap += (ao->transition->advance_query
                                   * ao->length);
                exon_target_gap += (ao->transition->advance_target
                                   * ao->length);
                break;
            case C4_Label_5SS:
                if(report_on_query){
                    g_assert(ao->transition->advance_query == 2);
                    g_assert(ao->transition->advance_target == 0);
                } else {
                    g_assert(ao->transition->advance_query == 0);
                    g_assert(ao->transition->advance_target == 2);
                    }
                if(in_exon){
                    Alignment_display_gff_utr(alignment,
                            query, target, report_on_query, post_cds,
                            cds_query_start, cds_target_start,
                            cds_query_end, cds_target_end,
                            exon_query_start, exon_target_start,
                            query_pos, target_pos, fp);
                    Alignment_display_gff_exon(alignment,
                          query, target, translate, report_on_query,
                          query_pos, target_pos,
                          exon_query_start, exon_target_start,
                          exon_query_gap, exon_target_gap,
                          exon_query_frameshift,
                          exon_target_frameshift, user_data, fp);
                    in_exon = FALSE;
                    }
                attribute_list = g_ptr_array_new();
                g_ptr_array_add(attribute_list,
                                g_strdup_printf("intron_id %d", intron_id+1));
                g_ptr_array_add(attribute_list,
                                g_strdup_printf("splice_site \"%c%c\"",
                   report_on_query?Sequence_get_symbol(query, query_pos)
                                  :Sequence_get_symbol(target, target_pos),
                   report_on_query?Sequence_get_symbol(query, query_pos+1)
                                  :Sequence_get_symbol(target, target_pos+1)));
                Alignment_display_gff_line(alignment, query, target,
                               report_on_query, "splice5",
                               query_pos, target_pos,
                               query_pos+2, target_pos+2,
                               FALSE, 0, FALSE, 0, attribute_list, fp);
                Alignment_free_attribute_list(attribute_list);
                intron_length = 0;
                break;
            case C4_Label_3SS:
                if(report_on_query){
                    g_assert(ao->transition->advance_query == 2);
                    g_assert(ao->transition->advance_target == 0);
                } else {
                    g_assert(ao->transition->advance_query == 0);
                    g_assert(ao->transition->advance_target == 2);
                    }
                if(in_exon){
                    Alignment_display_gff_utr(alignment,
                            query, target, report_on_query, post_cds,
                            cds_query_start, cds_target_start,
                            cds_query_end, cds_target_end,
                            exon_query_start, exon_target_start,
                            query_pos, target_pos, fp);
                    Alignment_display_gff_exon(alignment,
                          query, target, translate, report_on_query,
                          query_pos, target_pos,
                          exon_query_start, exon_target_start,
                          exon_query_gap, exon_target_gap,
                          exon_query_frameshift,
                          exon_target_frameshift, user_data, fp);
                    in_exon = FALSE;
                    }
                if(gene_orientation == '+'){
                    attribute_list = g_ptr_array_new();
                    g_ptr_array_add(attribute_list,
                                    g_strdup_printf("intron_id %d", ++intron_id));
                    Alignment_display_gff_line(alignment, query, target,
                                   report_on_query, "intron",
                                   query_pos-intron_length-2,
                                   target_pos-intron_length-2,
                                   query_pos+2, target_pos+2,
                                   FALSE, 0, FALSE, 0, attribute_list, fp);
                    Alignment_free_attribute_list(attribute_list);
                    }
                attribute_list = g_ptr_array_new();
                g_ptr_array_add(attribute_list,
                                g_strdup_printf("intron_id %d", intron_id-1));
                g_ptr_array_add(attribute_list,
                                g_strdup_printf("splice_site \"%c%c\"",
                   report_on_query?Sequence_get_symbol(query, query_pos)
                                  :Sequence_get_symbol(target, target_pos),
                   report_on_query?Sequence_get_symbol(query, query_pos+1)
                                  :Sequence_get_symbol(target, target_pos+1)));
                Alignment_display_gff_line(alignment, query, target,
                               report_on_query, "splice3",
                               query_pos, target_pos,
                               query_pos+2, target_pos+2,
                               FALSE, 0, FALSE, 0, attribute_list, fp);
                Alignment_free_attribute_list(attribute_list);
                intron_length = 0;
                break;
            case C4_Label_INTRON:
                if(report_on_query){
                    g_assert(ao->transition->advance_query == 1);
                    g_assert(ao->transition->advance_target == 0);
                } else {
                    g_assert(ao->transition->advance_query == 0);
                    g_assert(ao->transition->advance_target == 1);
                    }
                intron_length += ao->length;
                break;
            case C4_Label_FRAMESHIFT:
                exon_query_frameshift
                    += (ao->transition->advance_query * ao->length);
                exon_target_frameshift
                    += (ao->transition->advance_target * ao->length);
                break;
            case C4_Label_NER:
                g_error("Unexpected NER for gff gene output");
                break;
            default:
                g_error("Unknown C4_Label [%s]",
                        C4_Label_get_name(ao->transition->label));
                break;
            }
        query_pos += (ao->transition->advance_query * ao->length);
        target_pos += (ao->transition->advance_target * ao->length);
        }
    if(in_exon){
        if(cds_query_end != -1){
            if(cds_query_end != query_pos){ /* Have 3' UTR */
                curr_utr_query_start = MAX(exon_query_start, cds_query_end);
                curr_utr_target_start = MAX(exon_target_start, cds_target_end);
                Alignment_display_gff_line(alignment, query, target,
                               report_on_query, "utr3b",
                               curr_utr_query_start,
                               curr_utr_target_start,
                               query_pos, target_pos,
                               FALSE, 0, FALSE, 0, NULL, fp);
            } else {
                Alignment_display_gff_line(alignment, query, target,
                               report_on_query, "cds",
                               exon_query_start, exon_target_start,
                               query_pos, target_pos,
                               FALSE, 0, FALSE, 0, NULL, fp);
                }
            }
        Alignment_display_gff_exon(alignment,
              query, target, translate, report_on_query,
              query_pos, target_pos,
              exon_query_start, exon_target_start,
              exon_query_gap, exon_target_gap,
              exon_query_frameshift, exon_target_frameshift, user_data, fp);
        }
    return;
    }

static void Alignment_display_gff_similarity(Alignment *alignment,
        Sequence *query, Sequence *target, gboolean report_on_query,
        gint result_id, FILE *fp){
    register gint i, qp, tp;
    register gint query_pos = alignment->region->query_start,
                  target_pos = alignment->region->target_start;
    register AlignmentOperation *ao;
    register GPtrArray *attribute_list = g_ptr_array_new();
    g_ptr_array_add(attribute_list,
                    g_strdup_printf("alignment_id %d", result_id));
    if(report_on_query)
        g_ptr_array_add(attribute_list,
                    g_strdup_printf("Target %s", target->id));
    else
        g_ptr_array_add(attribute_list,
                    g_strdup_printf("Query %s", query->id));
    /* Align <seq_start> <target_start> [<length>] ; */
    for(i = 1; i < alignment->operation_list->len; i++){
        ao = alignment->operation_list->pdata[i];
        switch(ao->transition->label){
            case C4_Label_MATCH:
                qp = query_pos;
                tp = target_pos;
                if(query->strand == Sequence_Strand_REVCOMP)
                    qp = query->len - qp;
                if(target->strand == Sequence_Strand_REVCOMP)
                    tp = target->len - tp;
                if(report_on_query){
                    g_ptr_array_add(attribute_list,
                     g_strdup_printf("Align %d %d %d", qp+1, tp+1,
                                     ao->length*ao->transition->advance_query));
                } else {
                    g_ptr_array_add(attribute_list,
                     g_strdup_printf("Align %d %d %d", tp+1, qp+1,
                                     ao->length*ao->transition->advance_target));
                    }
                break;
            case C4_Label_NONE:
            case C4_Label_GAP:
            case C4_Label_NER:
            case C4_Label_5SS:
            case C4_Label_3SS:
            case C4_Label_INTRON:
            case C4_Label_SPLIT_CODON:
            case C4_Label_FRAMESHIFT:
                break;
            default:
                g_error("Unknown C4_Label [%s]",
                        C4_Label_get_name(ao->transition->label));
                break;
            }
        query_pos += (ao->transition->advance_query * ao->length);
        target_pos += (ao->transition->advance_target * ao->length);
        }
    Alignment_display_gff_line(alignment, query, target, report_on_query,
                               "similarity",
                               alignment->region->query_start,
                               alignment->region->target_start,
                               Region_query_end(alignment->region),
                               Region_target_end(alignment->region),
                               TRUE, alignment->score,
                               FALSE, 0, attribute_list, fp);
    Alignment_free_attribute_list(attribute_list);
    return;
    }

/**/

void Alignment_display_gff(Alignment *alignment,
                           Sequence *query, Sequence *target,
                           Translate *translate,
                           gboolean report_on_query,
                           gboolean report_on_genomic,
                           gint result_id, gpointer user_data, FILE *fp){
    fprintf(fp, "# --- START OF GFF DUMP ---\n#\n");
    Alignment_display_gff_header(alignment, query, target,
                                 report_on_query, fp);
    if(report_on_genomic){
        Alignment_display_gff_gene(alignment, query, target, translate,
                                   report_on_query, result_id, user_data, fp);
        }
    Alignment_display_gff_similarity(alignment, query, target,
                                     report_on_query, result_id, fp);
    fprintf(fp, "# --- END OF GFF DUMP ---\n#\n");
    return;
    }

/**/

/*
features: gene,intron,exon,splice[35],similarity
*/

/**/

#ifndef G_DISABLE_ASSERT
static gboolean Alignment_has_valid_alignment(Alignment *alignment,
                                              gpointer user_data){
    register gint query_pos = alignment->region->query_start,
                  target_pos = alignment->region->target_start;
    register C4_Score score;
    register C4_Calc *calc;
    register AlignmentOperation *alignment_operation;
    register C4_State *state;
    register C4_Transition *transition;
    register C4_Shadow *shadow;
    register gint i, j, k, t;
    register C4_Score *start_cell, *cell
     = g_new0(C4_Score, 1+alignment->model->total_shadow_designations);
/* FIXME: temp */
/* #define DEBUG */
#ifdef DEBUG
    g_message("Check model [%s] align score [%d]", alignment->model->name,
            alignment->score);
    for(i = 0; i < alignment->model->transition_list->len; i++){
        transition = alignment->model->transition_list->pdata[i];
        g_message("Check transition (%d) [%s] has dst_shad[%d] (%s:%s)",
                i, transition->name, transition->dst_shadow_list->len,
                transition->input->name, transition->output->name);
        }
    Region_print(alignment->region, "Alignment_has_valid_alignment");
#endif /* DEBUG */
    if(alignment->model->init_func)
        alignment->model->init_func(alignment->region, user_data);
    for(i = 0; i < alignment->model->calc_list->len; i++){
        calc = alignment->model->calc_list->pdata[i];
        if(calc && calc->init_func)
            calc->init_func(alignment->region, user_data);
        }
    if(alignment->model->start_state->cell_start_func){
        start_cell = alignment->model->start_state->cell_start_func
            (query_pos, target_pos, user_data);
#ifdef DEBUG
        g_print("USING start cell from (%d,%d): ",
                query_pos, target_pos);
#endif /* DEBUG */
        for(i = 0;
            i < (1 + alignment->model->total_shadow_designations); i++){
            cell[i] = start_cell[i];
#ifdef DEBUG
            g_print("[%d]", cell[i]);
#endif /* DEBUG */
            }
#ifdef DEBUG
        g_print("\n");
#endif /* DEBUG */
        }
    score = cell[0];
#ifdef DEBUG
    g_message("start with [%d] (%s)", score,
            alignment->model->name);
#endif /* DEBUG */
    for(i = 0; i < alignment->operation_list->len; i++){
        alignment_operation = alignment->operation_list->pdata[i];
        for(j = 0; j < alignment_operation->length; j++){
            state = alignment_operation->transition->input;
            for(k = 0; k < state->src_shadow_list->len; k++){
                shadow = state->src_shadow_list->pdata[k];
/* #define DEBUG_SHADOWS */
#ifdef DEBUG_SHADOWS
                g_message("start_func on [%s:%s:%s] (%d,%d) shadow(%d)",
                          state->name,
                          alignment_operation->transition->name,
                          shadow->name,
                          query_pos, target_pos, k);
#endif /* DEBUG_SHADOWS */
                cell[1 + shadow->designation] = shadow->start_func(
                  query_pos, target_pos, user_data);
                }
            transition = alignment_operation->transition;
            for(k = 0; k < transition->dst_shadow_list->len; k++){
                shadow = transition->dst_shadow_list->pdata[k];
#ifdef DEBUG_SHADOWS
                g_message("end_func on [%s:%s] (%d,%d) shadow(%d)=(%d)",
                          transition->name,
                          shadow->name,
                          query_pos, target_pos, k,
                          cell[1 + shadow->designation]);
#endif /* DEBUG_SHADOWS */
                shadow->end_func(cell[1 + shadow->designation],
                  query_pos, target_pos, user_data);
                }
            /**/
            t = C4_Calc_score(alignment_operation->transition->calc,
                              query_pos, target_pos, user_data);
            score += t;
#ifdef DEBUG
                g_message("add on [%d] (%s) => [%d] at (%d,%d) {%d,%d}",
                    t, alignment_operation->transition->name,
                    score, query_pos, target_pos, i, j);
#endif /* DEBUG */
            query_pos
                += alignment_operation->transition->advance_query;
            target_pos
                += alignment_operation->transition->advance_target;
            }
        }
    for(i = 0; i < alignment->model->calc_list->len; i++){
        calc = alignment->model->calc_list->pdata[i];
        if(calc && calc->exit_func)
            calc->exit_func(alignment->region, user_data);
        }
    if(alignment->model->exit_func)
        alignment->model->exit_func(alignment->region, user_data);
#ifdef DEBUG
    g_message("end with [%d]", score);
#endif /* DEBUG */
#ifdef DEBUG
    g_message("CHECKIT: [%d-%d=%d],[%d] [%d-%d=%d],[%d]",
             query_pos,
             alignment->region->query_start,
             (query_pos-alignment->region->query_start),
             alignment->region->query_length,
              target_pos,
              alignment->region->target_start,
             (target_pos-alignment->region->target_start),
             alignment->region->target_length);
#endif /* DEBUG */
    if(score != alignment->score) /* FIXME: temp */
        g_warning("Score difference DP[%d] TB[%d]",
                alignment->score, score);
    g_assert((query_pos-alignment->region->query_start)
             == alignment->region->query_length);
    g_assert((target_pos-alignment->region->target_start)
             == alignment->region->target_length);
    g_assert(score == alignment->score);
    g_free(cell);
    return TRUE;
    }
/* FIXME: Should change equality test in case the score is a float
 */
#endif /* G_DISABLE_ASSERT */

static gboolean Alignment_has_valid_path(Alignment *alignment){
    register C4_State *first_state, *last_state;
    register AlignmentOperation *alignment_operation, *prev;
    register gint i;
    g_assert(alignment);
    g_assert(alignment->operation_list);
    g_assert(alignment->operation_list->len > 0);
    alignment_operation = alignment->operation_list->pdata[0];
    first_state = alignment_operation->transition->input;
    g_assert(first_state == alignment->model->start_state->state);
    alignment_operation = alignment->operation_list->pdata
                         [alignment->operation_list->len-1];
    last_state = alignment_operation->transition->output;
    g_assert(last_state == alignment->model->end_state->state);
    for(i = 1; i < alignment->operation_list->len; i++){
        alignment_operation = alignment->operation_list->pdata[i];
        prev = alignment->operation_list->pdata[i-1];
        g_assert(alignment_operation->transition->input
              == prev->transition->output);
        }
    return TRUE;
    }

static gboolean Alignment_is_within_scope(Alignment *alignment,
                                          Region *seq_region){
    register gboolean start_at_query_edge,
                      start_at_target_edge,
                      end_at_query_edge,
                      end_at_target_edge;
    start_at_query_edge
        = (seq_region->query_start == alignment->region->query_start);
    start_at_target_edge
        = (seq_region->target_start == alignment->region->target_start);
    end_at_query_edge = (Region_query_end(alignment->region)
                         == Region_query_end(seq_region));
    end_at_target_edge = (Region_target_end(alignment->region)
                         == Region_target_end(seq_region));
    switch(alignment->model->start_state->scope){
        case C4_Scope_ANYWHERE:
            break;
        case C4_Scope_EDGE:
            g_assert(start_at_query_edge || start_at_target_edge);
            break;
        case C4_Scope_QUERY:
            g_assert(start_at_query_edge);
            break;
        case C4_Scope_TARGET:
            g_assert(start_at_target_edge);
            break;
        case C4_Scope_CORNER:
            g_assert(start_at_query_edge && start_at_target_edge);
            break;
        }
    switch(alignment->model->end_state->scope){
        case C4_Scope_ANYWHERE:
            break;
        case C4_Scope_EDGE:
            g_assert(end_at_query_edge || end_at_target_edge);
            break;
        case C4_Scope_QUERY:
            g_assert(end_at_query_edge);
            break;
        case C4_Scope_TARGET:
            g_assert(end_at_target_edge);
            break;
        case C4_Scope_CORNER:
            g_assert(end_at_query_edge && end_at_target_edge);
            break;
        }
    return TRUE;
    }

gboolean Alignment_is_valid(Alignment *alignment, Region *seq_region,
                            gpointer user_data){
    g_assert(alignment);
    g_assert(seq_region);
    g_assert(Region_is_within(seq_region, alignment->region));
    g_assert(Alignment_is_within_scope(alignment, seq_region));
    g_assert(Alignment_has_valid_path(alignment));
    g_assert(Alignment_has_valid_alignment(alignment, user_data));
    return TRUE;
    }

/**/

void Alignment_import_derived(Alignment *alignment,
                              Alignment *to_add,
                              C4_DerivedModel *derived_model){
    register gint i;
    register AlignmentOperation *alignment_operation;
    register C4_Transition *transition;
    g_assert(alignment->model == derived_model->original);
    g_assert(to_add->model == derived_model->derived);
    for(i = 0; i < to_add->operation_list->len; i++){
        alignment_operation = to_add->operation_list->pdata[i];
        transition = derived_model->transition_map
                               [alignment_operation->transition->id];
        Alignment_add(alignment, transition,
                      alignment_operation->length);
        }
    return;
    }

