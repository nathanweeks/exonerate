/****************************************************************\
*                                                                *
*  Match : A module for pairwise symbol comparison               *
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

#include <ctype.h> /* For toupper() */

#include "match.h"

static gchar *Match_Argument_parse_submat(gchar *arg_string,
                                          gpointer data){
    register Submat **submat = (Submat**)data;
    gchar *submat_str;
    register gchar *ret_val = Argument_parse_string(arg_string,
                                                   &submat_str);
    if(ret_val)
        return ret_val;
    (*submat) = Submat_create(submat_str);
    return NULL;
    }

static void Match_Argument_cleanup(gpointer user_data){
    Match_destroy_all();
    return;
    }

Match_ArgumentSet *Match_ArgumentSet_create(Argument *arg){
    register ArgumentSet *as;
    static Match_ArgumentSet mas = {0};
    if(arg){
        as = ArgumentSet_create("Symbol Comparison Options");
        ArgumentSet_add_option(as, 0, "softmaskquery", NULL,
            "Allow softmasking on the query sequence", "FALSE",
            Argument_parse_boolean, &mas.softmasked_query);
        ArgumentSet_add_option(as, 0, "softmasktarget", NULL,
            "Allow softmasking on the target sequence", "FALSE",
            Argument_parse_boolean, &mas.softmasked_target);
        /**/
        ArgumentSet_add_option(as, 'd', "dnasubmat", "name",
            "DNA substitution matrix", "nucleic",
            Match_Argument_parse_submat, &mas.dna_submat);
        ArgumentSet_add_option(as, 'p', "proteinsubmat", "name",
            "Protein substitution matrix", "blosum62",
            Match_Argument_parse_submat, &mas.protein_submat);
        /**/
        Argument_absorb_ArgumentSet(arg, as);
        Argument_add_cleanup(arg, Match_Argument_cleanup, NULL);
        /**/
    } else {
        if(!mas.translate){
            mas.translate = Translate_create(FALSE);
            /* FIXME this should be freed somewhere */
            }
        }
    return &mas;
    }

/**/

#define Match_Type_pair(qt, tt) ((qt) | ((tt) << 8))

Match_Type Match_Type_find(Alphabet_Type query_type,
                           Alphabet_Type target_type,
                           gboolean translate_both){
    register Match_Type type = 0;
    switch(Match_Type_pair(query_type, target_type)){
        case Match_Type_pair(Alphabet_Type_DNA, Alphabet_Type_DNA):
            type = translate_both?Match_Type_CODON2CODON
                                 :Match_Type_DNA2DNA;
            break;
        case Match_Type_pair(Alphabet_Type_PROTEIN,
                             Alphabet_Type_PROTEIN):
            type = Match_Type_PROTEIN2PROTEIN;
            break;
        case Match_Type_pair(Alphabet_Type_DNA,
                             Alphabet_Type_PROTEIN):
            type = Match_Type_DNA2PROTEIN;
            break;
        case Match_Type_pair(Alphabet_Type_PROTEIN,
                             Alphabet_Type_DNA):
            type = Match_Type_PROTEIN2DNA;
            break;
        default:
            g_error("Do not have a match type for [%s] [%s]",
                    Alphabet_Type_get_name(query_type),
                    Alphabet_Type_get_name(target_type));
            break;
        }
    return type;
    }

gchar *Match_Type_get_name(Match_Type type){
    register gchar *name = NULL;
    switch(type){
        case Match_Type_DNA2DNA:
            name = "dna2dna";
            break;
        case Match_Type_PROTEIN2PROTEIN:
            name = "protein2protein";
            break;
        case Match_Type_DNA2PROTEIN:
            name = "dna2protein";
            break;
        case Match_Type_PROTEIN2DNA:
            name = "protein2dna";
            break;
        case Match_Type_CODON2CODON:
            name = "codon";
            break;
        default:
            g_error("Bad Match_Type [%d]", type);
            break;
        }
    return name;
    }

static Match_Type Match_Type_mirror(Match_Type type){
    register Match_Type mirror = 0;
    switch(type){
        case Match_Type_DNA2DNA:         /*fallthrough*/
        case Match_Type_PROTEIN2PROTEIN: /*fallthrough*/
        case Match_Type_CODON2CODON:     /*fallthrough*/
            mirror = type;
            break;
        case Match_Type_DNA2PROTEIN:
            mirror = Match_Type_PROTEIN2DNA;
            break;
        case Match_Type_PROTEIN2DNA:
            mirror = Match_Type_DNA2PROTEIN;
            break;
        default:
            g_error("Bad Match_Type for mirroring [%d]", type);
        }
    return mirror;
    }

/**/

static gboolean Match_Strand_pos_is_valid(Match_Strand *strand,
                                          gint pos, gint length){
    g_assert(pos >= 0);
    g_assert((pos+strand->advance) <= length);
    return TRUE;
    }

#define Match_seqpos_is_masked(sequence, pos)              \
    Alphabet_is_masked((sequence)->alphabet,               \
                       Sequence_get_symbol(sequence, pos))

/**/

static Match_Score Match_1_dna_self_func(Match_Strand *strand,
                                         Sequence *seq, guint pos){
    register gchar symbol;
    g_assert(Match_Strand_pos_is_valid(strand, pos, seq->len));
    symbol = Sequence_get_symbol(seq, pos);
    return Submat_lookup(strand->match->mas->dna_submat, symbol, symbol);
    }

static Match_Score Match_1_protein_self_func(Match_Strand *strand,
                                             Sequence *seq, guint pos){
    register gchar symbol;
    g_assert(Match_Strand_pos_is_valid(strand, pos, seq->len));
    symbol = Sequence_get_symbol(seq, pos);
    return Submat_lookup(strand->match->mas->protein_submat, symbol, symbol);
    }

static gboolean Match_1_mask_func(Match_Strand *strand,
                                  Sequence *seq, guint pos){
    g_assert(Match_Strand_pos_is_valid(strand, pos, seq->len));
    return Match_seqpos_is_masked(seq, pos);
    }

/**/

static Match_Score Match_3_translate_self_func(Match_Strand *strand,
                                            Sequence *seq, guint pos){
    register gchar symbol;
    g_assert(Match_Strand_pos_is_valid(strand, pos, seq->len));
    symbol = Translate_base(strand->match->mas->translate,
                            Sequence_get_symbol(seq, pos),
                            Sequence_get_symbol(seq, pos+1),
                            Sequence_get_symbol(seq, pos+2));
    return Submat_lookup(strand->match->mas->protein_submat,
                         symbol, symbol);
    }

#if 0
static Match_Score Match_3_codon_self_func(Match_Strand *strand,
                                           Sequence *seq, guint pos){
    register gint codon;
    g_assert(Match_Strand_pos_is_valid(strand, pos, seq->len));
    codon = CodonSubmat_get_codon_base(strand->match->mas->codon_submat,
                                       Sequence_get_symbol(seq, pos),
                                       Sequence_get_symbol(seq, pos+1),
                                       Sequence_get_symbol(seq, pos+2));
    return CodonSubmat_lookup_codon(strand->match->mas->codon_submat,
                                    codon, codon);
    }
#endif /* 0 */

static gboolean Match_3_mask_func(Match_Strand *strand,
                                  Sequence *seq, guint pos){
    g_assert(Match_Strand_pos_is_valid(strand, pos, seq->len));
    if((Match_seqpos_is_masked(seq, pos))
    || (Match_seqpos_is_masked(seq, pos+1))
    || (Match_seqpos_is_masked(seq, pos+2)))
        return TRUE;
    return FALSE;
    }

/**/

static gchar Match_get_display_symbol(Submat *submat,
                                      gchar query_symbol,
                                      gchar target_symbol){
    register gint score;
    if(toupper(query_symbol) == toupper(target_symbol))
        return '|';
    score = Submat_lookup(submat, query_symbol, target_symbol);
    if(score == 0)
        return '.';
    if(score > 0)
        return ':';
    return ' ';
    }

static gboolean Match_pos_pair_is_valid(Match *match,
                             Sequence *query, Sequence *target,
                             guint query_pos, guint target_pos){
    g_assert(Match_Strand_pos_is_valid(match->query, query_pos,
                                       query->len));
    g_assert(Match_Strand_pos_is_valid(match->target, target_pos,
                                       target->len));
    return TRUE;
    }

/**/

static Match_Score Match_invalid_split_score_func(Match *match,
                              Sequence *query, Sequence *target,
                              guint qp1, guint qp2, guint qp3,
                              guint tp1, guint tp2, guint tp3){
    g_error("Match [%s] cannot use a split score func",
            Match_Type_get_name(match->type));
    return 0;
    }

static void Match_invalid_split_display_func(Match *match,
                              Sequence *query, Sequence *target,
                              guint qp1, guint qp2, guint qp3,
                              guint tp1, guint tp2, guint tp3,
                              gchar *display_str){
    g_error("Match [%s] cannot use a split display func",
            Match_Type_get_name(match->type));
    return;
    }

/**/

static Match_Score Match_1_1_dna_score_func(Match *match,
                              Sequence *query, Sequence *target,
                              guint query_pos, guint target_pos){
    g_assert(Match_pos_pair_is_valid(match, query, target,
                                     query_pos, target_pos));
    if(query->annotation)
        if(query->alphabet->type == Alphabet_Type_DNA)
            if((query_pos >= query->annotation->cds_start)
            && (query_pos < (query->annotation->cds_start
                            +query->annotation->cds_length)))
                return MATCH_IMPOSSIBLY_LOW_SCORE;
    return Submat_lookup(match->mas->dna_submat,
                         Sequence_get_symbol(query, query_pos),
                         Sequence_get_symbol(target, target_pos));
    }

static Match_Score Match_1_1_protein_score_func(Match *match,
                              Sequence *query, Sequence *target,
                              guint query_pos, guint target_pos){
    g_assert(Match_pos_pair_is_valid(match, query, target,
                                     query_pos, target_pos));
    return Submat_lookup(match->mas->protein_submat,
                         Sequence_get_symbol(query, query_pos),
                         Sequence_get_symbol(target, target_pos));
    }

static gchar *Match_1_1_dna_score_macro(void){
    return "((ud->query->annotation)\n"
           "&& (%QP >= ud->query->annotation->cds_start)\n"
           "&& (%QP < (ud->query->annotation->cds_start\n"
           "          +ud->query->annotation->cds_length)))\n"
           "?MATCH_IMPOSSIBLY_LOW_SCORE\n"
           ":Submat_lookup(ud->mas->dna_submat,\n"
           "     Sequence_get_symbol(ud->query, %QP),\n"
           "     Sequence_get_symbol(ud->target, %TP))";
    }

static gchar *Match_1_1_protein_score_macro(void){
    return "Submat_lookup(ud->mas->protein_submat,\n"
           "     Sequence_get_symbol(ud->query, %QP),\n"
           "     Sequence_get_symbol(ud->target, %TP))";
    }

static void Match_1_1_display_func(Match *match,
                                   Sequence *query, Sequence *target,
                                   guint query_pos, guint target_pos,
                                   gchar *display_str){
    register Submat *submat = (match->type == Match_Type_DNA2DNA)
                            ? match->mas->dna_submat
                            : match->mas->protein_submat;
    g_assert(Match_pos_pair_is_valid(match, query, target,
                                     query_pos, target_pos));
    display_str[0] = Match_get_display_symbol(submat,
                           Sequence_get_symbol(query, query_pos),
                           Sequence_get_symbol(target, target_pos));
    display_str[1] = '\0';
    return;
    }

/**/

static Match_Score Match_1_3_split_score_func(Match *match,
                              Sequence *query, Sequence *target,
                              guint qp1, guint qp2, guint qp3,
                              guint tp1, guint tp2, guint tp3){
    register gchar query_symbol, target_symbol;
    /* FIXME: need CDS checks here */
    query_symbol = Sequence_get_symbol(query, qp1);
    target_symbol = Translate_base(match->mas->translate,
                              Sequence_get_symbol(target, tp1),
                              Sequence_get_symbol(target, tp2),
                              Sequence_get_symbol(target, tp3));
    return Submat_lookup(match->mas->protein_submat, query_symbol,
                                                     target_symbol);
    }

static Match_Score Match_1_3_score_func(Match *match,
                              Sequence *query, Sequence *target,
                              guint query_pos, guint target_pos){
    g_assert(Match_pos_pair_is_valid(match, query, target,
                                     query_pos, target_pos));
    return Match_1_3_split_score_func(match, query, target,
            query_pos, -1, -1,
            target_pos, target_pos+1, target_pos+2);
    }

static gchar *Match_1_3_score_macro(void){
    return "Submat_lookup(ud->mas->protein_submat,\n"
           "              Sequence_get_symbol(ud->query, %QP),\n"
           "              Translate_base(ud->mas->translate,\n"
           "                  Sequence_get_symbol(ud->target, %TP),\n"
           "                  Sequence_get_symbol(ud->target, %TP+1),\n"
           "                  Sequence_get_symbol(ud->target, %TP+2)))\n";
    }

typedef struct {
    gchar *display_string;
    gchar  codon_a;
    gchar  codon_b;
    gchar  codon_c;
} Match_RT_Data;

static void Match_1_3_display_reverse_translate(gchar *dna, gint length,
                                                gpointer user_data){
    register Match_RT_Data *mrtd = user_data;
    if(dna[0] == mrtd->codon_a)
        mrtd->display_string[0] = '!';
    if(dna[1] == mrtd->codon_b)
        mrtd->display_string[1] = '!';
    if(dna[2] == mrtd->codon_c)
        mrtd->display_string[2] = '!';
    return;
    }

static void Match_1_3_split_display_func(Match *match,
                              Sequence *query, Sequence *target,
                              guint qp1, guint qp2, guint qp3,
                              guint tp1, guint tp2, guint tp3,
                              gchar *display_str){
    register gchar query_symbol, target_symbol, display_symbol;
    Match_RT_Data mrtd;
    gchar aaseq[2];
    query_symbol = Sequence_get_symbol(query, qp1);
    target_symbol = Translate_base(match->mas->translate,
                                   Sequence_get_symbol(target, tp1),
                                   Sequence_get_symbol(target, tp2),
                                   Sequence_get_symbol(target, tp3));
    display_symbol = Match_get_display_symbol(match->mas->protein_submat,
                                              query_symbol, target_symbol);
    display_str[0] = display_symbol;
    display_str[1] = display_symbol;
    display_str[2] = display_symbol;
    display_str[3] = '\0';
    if(query_symbol != target_symbol){
        /* Label with reverse translation */
        mrtd.display_string = display_str;
        mrtd.codon_a = Sequence_get_symbol(target, tp1);
        mrtd.codon_b = Sequence_get_symbol(target, tp2);
        mrtd.codon_c = Sequence_get_symbol(target, tp3);
        aaseq[0] = query_symbol;
        aaseq[1] = '\0';
        Translate_reverse(match->mas->translate,
                          aaseq, 1, Match_1_3_display_reverse_translate,
                          &mrtd);
        }
    return;
    }

static void Match_1_3_display_func(Match *match,
                                   Sequence *query, Sequence *target,
                                   guint query_pos, guint target_pos,
                                   gchar *display_str){
    g_assert(Match_pos_pair_is_valid(match, query, target,
                                     query_pos, target_pos));
    Match_1_3_split_display_func(match, query, target,
                                 query_pos, -1, -1,
                                 target_pos, target_pos+1, target_pos+2,
                                 display_str);
    return;
    }

/**/

static Match_Score Match_3_1_split_score_func(Match *match,
                              Sequence *query, Sequence *target,
                              guint qp1, guint qp2, guint qp3,
                              guint tp1, guint tp2, guint tp3){
    return Match_1_3_split_score_func(Match_swap(match),
                                      target, query,
                                      tp1, tp2, tp3,
                                      qp1, qp2, qp3);
    }

static Match_Score Match_3_1_score_func(Match *match,
                                    Sequence *query, Sequence *target,
                                    guint query_pos, guint target_pos){
    return Match_1_3_score_func(Match_swap(match),
                                target, query, target_pos, query_pos);
    }

static gchar *Match_3_1_score_macro(void){
    return "Submat_lookup(ud->mas->protein_submat,\n"
           "              Translate_base(ud->mas->translate,\n"
           "                  Sequence_get_symbol(ud->query, %QP),\n"
           "                  Sequence_get_symbol(ud->query, %QP+1),\n"
           "                  Sequence_get_symbol(ud->query, %QP+2)),\n"
           "              Sequence_get_symbol(ud->target, %TP))\n";
    }

static void Match_3_1_split_display_func(Match *match,
                              Sequence *query, Sequence *target,
                              guint qp1, guint qp2, guint qp3,
                              guint tp1, guint tp2, guint tp3,
                              gchar *display_str){
    Match_1_3_split_display_func(Match_swap(match),
                                 target, query,
                                 tp1, tp2, tp3,
                                 qp1, qp2, qp3,
                                 display_str);
    return;
    }

static void Match_3_1_display_func(Match *match,
                                   Sequence *query, Sequence *target,
                                   guint query_pos, guint target_pos,
                                   gchar *display_str){
    Match_1_3_display_func(Match_swap(match),
                           target, query, target_pos, query_pos,
                           display_str);
    return;
    }

/**/

#if 0
static Match_Score Match_3_3_split_score_func(Match *match,
                              Sequence *query, Sequence *target,
                              guint qp1, guint qp2, guint qp3,
                              guint tp1, guint tp2, guint tp3){
    if(query->annotation)
        if(query->alphabet->type == Alphabet_Type_DNA)
            if((qp1 < query->annotation->cds_start)
            || (qp1 >= (query->annotation->cds_start
                       +query->annotation->cds_length))
            || ((qp1 % 3) != (query->annotation->cds_start % 3)))
                return MATCH_IMPOSSIBLY_LOW_SCORE;
    return CodonSubmat_lookup_base(match->mas->codon_submat,
                Sequence_get_symbol(query, qp1),
                Sequence_get_symbol(query, qp2),
                Sequence_get_symbol(query, qp3),
                Sequence_get_symbol(target, tp1),
                Sequence_get_symbol(target, tp2),
                Sequence_get_symbol(target, tp3));
    }
#endif /* 0 */

/* FIXME: temp */
static Match_Score Match_3_3_split_score_func(Match *match,
                              Sequence *query, Sequence *target,
                              guint qp1, guint qp2, guint qp3,
                              guint tp1, guint tp2, guint tp3){
    register gchar query_symbol, target_symbol;
    if(query->annotation)
        if(query->alphabet->type == Alphabet_Type_DNA)
            if((qp1 < query->annotation->cds_start)
            || (qp1 >= (query->annotation->cds_start
                       +query->annotation->cds_length))
            || ((qp1 % 3) != (query->annotation->cds_start % 3)))
                return MATCH_IMPOSSIBLY_LOW_SCORE;
    query_symbol = Translate_base(match->mas->translate,
                       Sequence_get_symbol(query, qp1),
                       Sequence_get_symbol(query, qp2),
                       Sequence_get_symbol(query, qp3));
    target_symbol = Translate_base(match->mas->translate,
                       Sequence_get_symbol(target, tp1),
                       Sequence_get_symbol(target, tp2),
                       Sequence_get_symbol(target, tp3));
    return Submat_lookup(match->mas->protein_submat,
                query_symbol, target_symbol);
    }

static Match_Score Match_3_3_codon_score_func(Match *match,
                           Sequence *query, Sequence *target,
                           guint query_pos, guint target_pos){
    g_assert(Match_pos_pair_is_valid(match, query, target,
                                     query_pos, target_pos));
    return Match_3_3_split_score_func(match, query, target,
            query_pos, query_pos+1, query_pos+2,
            target_pos, target_pos+1, target_pos+2);
    }

#if 0
static gchar *Match_3_3_score_macro(void){
    return "((ud->query->annotation)\n"
           " && (ud->query->alphabet->type == Alphabet_Type_DNA)\n"
           " && (  ((%QP) < ud->query->annotation->cds_start)\n"
           "    || ((%QP) >= (ud->query->annotation->cds_start\n"
           "               +ud->query->annotation->cds_length))\n"
           "    || (((%QP) %% 3) != (ud->query->annotation->cds_start %% 3))))\n"
           "?MATCH_IMPOSSIBLY_LOW_SCORE\n"
           ":CodonSubmat_lookup_base(ud->mas->codon_submat,\n"
           "                  Sequence_get_symbol(ud->query, %QP),\n"
           "                  Sequence_get_symbol(ud->query, %QP+1),\n"
           "                  Sequence_get_symbol(ud->query, %QP+2),\n"
           "                  Sequence_get_symbol(ud->target, %TP),\n"
           "                  Sequence_get_symbol(ud->target, %TP+1),\n"
           "                  Sequence_get_symbol(ud->target, %TP+2))\n";
    }
/* FIXME: optimisations:
 *        remove seq type check
 *        precompute cds_start % 3
 */
#endif /* 0 */

/* FIXME: temp */
static gchar *Match_3_3_score_macro(void){
    return "((ud->query->annotation)\n"
           " && (ud->query->alphabet->type == Alphabet_Type_DNA)\n"
           " && (  ((%QP) < ud->query->annotation->cds_start)\n"
           "    || ((%QP) >= (ud->query->annotation->cds_start\n"
           "               +ud->query->annotation->cds_length))\n"
           "    || (((%QP) %% 3) != (ud->query->annotation->cds_start %% 3))))\n"
           "?MATCH_IMPOSSIBLY_LOW_SCORE\n"
           ":Submat_lookup(ud->mas->protein_submat,\n"
           "               Translate_base(ud->mas->translate, \n"
           "                  Sequence_get_symbol(ud->query, %QP),\n"
           "                  Sequence_get_symbol(ud->query, %QP+1),\n"
           "                  Sequence_get_symbol(ud->query, %QP+2)),\n"
           "               Translate_base(ud->mas->translate, \n"
           "                  Sequence_get_symbol(ud->target, %TP),\n"
           "                  Sequence_get_symbol(ud->target, %TP+1),\n"
           "                  Sequence_get_symbol(ud->target, %TP+2)))\n";
    }
/* FIXME: optimisations:
 *        remove seq type check
 *        precompute cds_start % 3
 */

static void Match_3_3_split_display_func(Match *match,
                              Sequence *query, Sequence *target,
                              guint qp1, guint qp2, guint qp3,
                              guint tp1, guint tp2, guint tp3,
                              gchar *display_str){
    register gchar query_symbol, target_symbol, display_symbol;
    query_symbol = Translate_base(match->mas->translate,
                       Sequence_get_symbol(query, qp1),
                       Sequence_get_symbol(query, qp2),
                       Sequence_get_symbol(query, qp3));
    target_symbol = Translate_base(match->mas->translate,
                       Sequence_get_symbol(target, tp1),
                       Sequence_get_symbol(target, tp2),
                       Sequence_get_symbol(target, tp3));
    display_symbol = Match_get_display_symbol(
                           match->mas->protein_submat,
                           query_symbol, target_symbol);
    /**/
    display_str[0] = display_symbol;
    if(query_symbol == target_symbol){
        if(toupper(Sequence_get_symbol(query, qp1))
        != toupper(Sequence_get_symbol(target, tp1)))
            display_str[0] = '+'; /* Mark Silent */
    } else {
        if(toupper(Sequence_get_symbol(query, qp1))
        == toupper(Sequence_get_symbol(target, tp1)))
            display_str[0] = '!'; /* Mark Conserved */
        }
    /**/
    display_str[1] = display_symbol;
    if(query_symbol == target_symbol){
        if(toupper(Sequence_get_symbol(query, qp2))
        != toupper(Sequence_get_symbol(target, tp2)))
            display_str[1] = '+'; /* Mark Silent */
    } else {
        if(toupper(Sequence_get_symbol(query, qp2))
        == toupper(Sequence_get_symbol(target, tp2)))
            display_str[1] = '!'; /* Mark Conserved */
        }
    /**/
    display_str[2] = display_symbol;
    if(query_symbol == target_symbol){
        if(toupper(Sequence_get_symbol(query, qp3))
        != toupper(Sequence_get_symbol(target, tp3)))
            display_str[2] = '+'; /* Mark Silent */
    } else {
        if(toupper(Sequence_get_symbol(query, qp3))
        == toupper(Sequence_get_symbol(target, tp3)))
            display_str[2] = '!'; /* Mark Conserved */
        }
    /**/
    display_str[3] = '\0';
    return;
    }

static void Match_3_3_display_func(Match *match,
                                   Sequence *query, Sequence *target,
                                   guint query_pos, guint target_pos,
                                   gchar *display_str){
    g_assert(Match_pos_pair_is_valid(match, query, target,
                                     query_pos, target_pos));
    Match_3_3_split_display_func(match, query, target,
            query_pos, query_pos+1, query_pos+2,
            target_pos, target_pos+1, target_pos+2, display_str);
    return;
    }

/**/

static Match_Strand *Match_Strand_create(Alphabet_Type alphabet_type,
       gboolean is_softmasked, gboolean is_translated,
       guint advance, Match *match){
    register Match_Strand *strand = g_new(Match_Strand, 1);
    strand->alphabet = Alphabet_create(alphabet_type, is_softmasked);
    strand->advance = advance;
    strand->is_translated = is_translated;
    switch(advance){
        case 1:
            strand->self_func = (alphabet_type == Alphabet_Type_DNA)
                              ? Match_1_dna_self_func
                              : Match_1_protein_self_func;
            strand->mask_func = Match_1_mask_func;
            break;
        case 3:
            /*
            strand->self_func = is_translated
                              ? Match_3_translate_self_func
                              : Match_3_codon_self_func;
                              */
            strand->self_func = Match_3_translate_self_func; /* FIXME: temp */
            strand->mask_func = Match_3_mask_func;
            break;
        default:
            g_error("Bad advance [%d]", advance);
            break;
        }
    strand->match = match;
    return strand;
    }

static void Match_Strand_destroy(Match_Strand *strand){
    Alphabet_destroy(strand->alphabet);
    g_free(strand);
    return;
    }

void Match_Strand_get_raw(Match_Strand *strand, Sequence *sequence,
                          guint pos, gchar *result){
    register gint i;
    g_assert(result);
    for(i = 0; i < strand->advance; i++)
        result[i] = Sequence_get_symbol(sequence, pos+i);
    result[i] = '\0';
    return;
    }

/**/

gchar *Match_Type_get_score_macro(Match_Type type){
    register gchar *macro = NULL;
    switch(type){
        case Match_Type_DNA2DNA:
            macro = Match_1_1_dna_score_macro();
            break;
        case Match_Type_PROTEIN2PROTEIN:
            macro = Match_1_1_protein_score_macro();
            break;
        case Match_Type_DNA2PROTEIN:
            macro = Match_3_1_score_macro();
            break;
        case Match_Type_PROTEIN2DNA:
            macro = Match_1_3_score_macro();
            break;
        case Match_Type_CODON2CODON:
            macro = Match_3_3_score_macro();
            break;
        default:
            g_error("Bad Match_Type [%d]", type);
            break;
        }
    g_assert(macro);
    return macro;
    }

/**/

static Match *Match_create_without_mirror(Match_Type type){
    register Match *match = g_new(Match, 1);
    register gboolean softmask_comparison;
    match->mas = Match_ArgumentSet_create(NULL);
    match->type = type;
    g_assert(match->mas->dna_submat);
    g_assert(match->mas->protein_submat);
    softmask_comparison = match->mas->softmasked_query
                       || match->mas->softmasked_target;
    switch(type){
        case Match_Type_DNA2DNA:
            match->query = Match_Strand_create(Alphabet_Type_DNA,
                   match->mas->softmasked_query, FALSE, 1, match);
            match->target = Match_Strand_create(Alphabet_Type_DNA,
                   match->mas->softmasked_target, FALSE, 1, match);
            match->comparison_alphabet
                      = Alphabet_create(Alphabet_Type_DNA,
                                        softmask_comparison);
            match->score_func = Match_1_1_dna_score_func;
            match->display_func = Match_1_1_display_func;
            match->split_score_func = Match_invalid_split_score_func;
            match->split_display_func = Match_invalid_split_display_func;
            break;
        case Match_Type_CODON2CODON:
            /* FIXME: should do this for mixed models only ?? */
#if 0
            g_assert(match->mas->dna_submat);
            g_assert(match->mas->codon_submat);
            CodonSubmat_add_nucleic(match->mas->codon_submat,
                                    match->mas->dna_submat);
#endif /* 0 */
            /**/
            match->query = Match_Strand_create(Alphabet_Type_DNA,
                     match->mas->softmasked_query, FALSE, 3, match);
            match->target = Match_Strand_create(Alphabet_Type_DNA,
                     match->mas->softmasked_target, FALSE, 3, match);
            match->comparison_alphabet
                = Alphabet_create(Alphabet_Type_PROTEIN,
                                  softmask_comparison);
            match->score_func = Match_3_3_codon_score_func;
            match->display_func = Match_3_3_display_func;
            match->split_score_func = Match_3_3_split_score_func;
            match->split_display_func = Match_3_3_split_display_func;
            break;
        case Match_Type_PROTEIN2PROTEIN:
            match->query = Match_Strand_create(Alphabet_Type_PROTEIN,
                     match->mas->softmasked_query, FALSE, 1, match);
            match->target = Match_Strand_create(Alphabet_Type_PROTEIN,
                     match->mas->softmasked_target, FALSE, 1, match);
            match->comparison_alphabet
                = Alphabet_create(Alphabet_Type_PROTEIN,
                                  softmask_comparison);
            match->score_func = Match_1_1_protein_score_func;
            match->display_func = Match_1_1_display_func;
            match->split_score_func = Match_invalid_split_score_func;
            match->split_display_func = Match_invalid_split_display_func;
            break;
        case Match_Type_DNA2PROTEIN:
            match->query = Match_Strand_create(Alphabet_Type_DNA,
                     match->mas->softmasked_query, TRUE, 3, match);
            match->target = Match_Strand_create(Alphabet_Type_PROTEIN,
                     match->mas->softmasked_target, FALSE, 1, match);
            match->comparison_alphabet
                = Alphabet_create(Alphabet_Type_PROTEIN,
                                  softmask_comparison);
            match->score_func = Match_3_1_score_func;
            match->display_func = Match_3_1_display_func;
            match->split_score_func = Match_3_1_split_score_func;
            match->split_display_func = Match_3_1_split_display_func;
            break;
        case Match_Type_PROTEIN2DNA:
            match->query = Match_Strand_create(Alphabet_Type_PROTEIN,
                   match->mas->softmasked_query, FALSE, 1, match);
            match->target = Match_Strand_create(Alphabet_Type_DNA,
                   match->mas->softmasked_target, TRUE, 3, match);
            /* protein_submat is used in target strand
             * because the comparison is done at the protein level
             */
            match->comparison_alphabet
                = Alphabet_create(Alphabet_Type_PROTEIN,
                                  softmask_comparison);
            match->score_func = Match_1_3_score_func;
            match->display_func = Match_1_3_display_func;
            match->split_score_func = Match_1_3_split_score_func;
            match->split_display_func = Match_1_3_split_display_func;
            break;
        default:
            g_error("Unknown Match Type [%d]", type);
            break;
        }
    match->mirror = NULL;
    return match;
    }

static Match *local_match_cache[Match_Type_TOTAL] = {0};

static Match *Match_create(Match_Type type){
    register Match *match;
    register Match_Type mirror_type = Match_Type_mirror(type);
    match = Match_create_without_mirror(type);
    match->mirror = Match_create_without_mirror(mirror_type);
    match->mirror->mirror = match;
    return match;
    }

Match *Match_find(Match_Type type){
    if(!local_match_cache[type]){
        local_match_cache[type] = Match_create(type);
        local_match_cache[Match_Type_mirror(type)]
            = local_match_cache[type]->mirror;
        }
    return local_match_cache[type];
    }

static void Match_destroy_without_mirror(Match *match){
    g_assert(match);
    match->mirror = NULL;
    Match_Strand_destroy(match->query);
    Match_Strand_destroy(match->target);
    Alphabet_destroy(match->comparison_alphabet);
    local_match_cache[match->type] = NULL;
    g_free(match);
    return;
    }

static void Match_destroy(Match *match){
    g_assert(match);
    Match_destroy_without_mirror(match->mirror);
    Match_destroy_without_mirror(match);
    return;
    }

void Match_destroy_all(void){
    register gint i;
    register Match *match;
    for(i = 0; i < Match_Type_TOTAL; i++){
        match = local_match_cache[i];
        if(match){
            local_match_cache[i] = NULL;
            local_match_cache[match->mirror->type] = NULL;
            Match_destroy(match);
            }
        }
    return;
    }

Match *Match_swap(Match *match){
    return match->mirror;
    }

Match_Score Match_max_score(Match *match){
    switch(match->type){
        case Match_Type_DNA2DNA:
            return Submat_max_score(match->mas->dna_submat);
        case Match_Type_PROTEIN2PROTEIN:
        case Match_Type_PROTEIN2DNA:
        case Match_Type_DNA2PROTEIN:
        case Match_Type_CODON2CODON:
            return Submat_max_score(match->mas->protein_submat);
        default:
            g_error("Unknown match type [%d]", match->type);
        }
    return 0;
    }

/**/

