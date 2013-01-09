/****************************************************************\
*                                                                *
*  Simple Alphabet Object                                        *
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
#include <strings.h> /* For strcasecmp() */
#include <ctype.h>  /* For toupper() */

#include "alphabet.h"

gchar *Alphabet_Argument_parse_alphabet_type(gchar *arg_string, gpointer data){
    register Alphabet_Type *dst_type = (Alphabet_Type*)data;
    gchar *type_str;
    register gchar *ret_val = Argument_parse_string(arg_string, &type_str);
    if(ret_val)
        return ret_val;
    (*dst_type) = Alphabet_name_get_type(type_str);
    return NULL;
    }

Alphabet_Type Alphabet_Type_guess(gchar *str){
    static const guchar *index = (guchar*)
    "----------------------------------------------------------------"
    "-A-C---G------N-----T------------A-C---G------N-----T-----------"
    "----------------------------------------------------------------"
    "----------------------------------------------------------------";
    register gint count = 0, len;
    register Alphabet_Type type;
    for(len = 0; str[len]; len++)
        if(index[(guchar)str[len]] != '-')
            count++;
    if((20.0*(gdouble)count) >= (17.0*(gdouble)len))
        type = Alphabet_Type_DNA;
    else
        type = Alphabet_Type_PROTEIN;
    return type;
    }

Alphabet_ArgumentSet *Alphabet_ArgumentSet_create(Argument *arg){
    register ArgumentSet *as;
    static Alphabet_ArgumentSet aas = {TRUE};
    if(arg){
        as = ArgumentSet_create("Alphabet Options");
        ArgumentSet_add_option(as, '\0', "useaatla", NULL,
           "Use three-letter abbreviation for AA names", "TRUE",
           Argument_parse_boolean, &aas.use_aa_tla);
        Argument_absorb_ArgumentSet(arg, as);
        }
    return &aas;
    }

/**/

#define Alphabet_Filter_macro(asterisk, upper_case, lower_case)        \
    (guchar*) ("------------------------------------------"            \
    asterisk "---------------------"                                   \
    "-" upper_case "------" lower_case "-----"                         \
    "----------------------------------------------------------------" \
    "----------------------------------------------------------------")

typedef enum {
    Alphabet_Filter_VALID_SEQ_SYMBOL,
    Alphabet_Filter_VALID_DNA_IUPAC,
    Alphabet_Filter_VALID_DNA_ACGTN,
    Alphabet_Filter_VALID_PROTEIN,
    Alphabet_Filter_CLEAN_DNA_ACGTN,
    Alphabet_Filter_CLEAN_DNA_IUPAC,
    Alphabet_Filter_CLEAN_PROTEIN,
    Alphabet_Filter_TO_UPPER,
    Alphabet_Filter_TO_LOWER,
    Alphabet_Filter_COMPLEMENT,
    Alphabet_Filter_SOFTMASK_DNA,
    Alphabet_Filter_SOFTMASK_PROTEIN,
    Alphabet_Filter_DNA_NOAMBIG,
    Alphabet_Filter_PROTEIN_NOAMBIG_SOFTMASK,
    Alphabet_Filter_PROTEIN_NOAMBIG,
    Alphabet_Filter_DNA_NOAMBIG_SOFTMASK,
    Alphabet_Filter_NUMBER_OF_FILTERS
} Alphabet_Filter;

static guchar *Alphabet_Filter_get(Alphabet_Filter type){
    static guchar *filter_set[Alphabet_Filter_NUMBER_OF_FILTERS] = {
    /* Alphabet_Filter_VALID_SEQ_SYMBOL: */
    Alphabet_Filter_macro("*", "ABCDEFGHI-KLMN-PQRSTUVWXYZ",
                               "abcdefghi-klmn-pqrstuvwxyz"),
    /* Alphabet_Filter_VALID_DNA_IUPAC: */
    Alphabet_Filter_macro("-", "ABCD--GH--K-MN---RST-VW-Y-",
                               "abcd--gh--k-mn---rst-vw-y-"),
    /* Alphabet_Filter_DNA_ACGTN: */
    Alphabet_Filter_macro("-", "A-C---G------N-----T------",
                               "a-c---g------n-----t------"),
    /* Alphabet_Filter_VALID_PROTEIN: */
    Alphabet_Filter_macro("*", "ABCDEFGHI-KLMN-PQRSTUVWXYZ",
                               "abcdefghi-klmn-pqrstuvwxyz"),
    /* Alphabet_Filter_CLEAN_DNA_ACGTN: */
    Alphabet_Filter_macro("-", "ANCNNNGNNNNNNNNNNNNTNNNNNN",
                               "ancnnngnnnnnnnnnnnntnnnnnn"),
    /* Alphabet_Filter_CLEAN_DNA_IUPAC: */
    Alphabet_Filter_macro("-", "ABCDNNGHNNKNMNNNNRSTNVWNYN",
                               "abcdnnghnnknmnnnnrstnvwnyn"),
    /* Alphabet_Filter_CLEAN_PROTEIN: */
    Alphabet_Filter_macro("*", "AXCDEFGHIXKLMNXPQRSTUVWXYX",
                               "axcdefghixklmnxpqrstuvwxyx"),
    /* Alphabet_Filter_TO_UPPER: */
    Alphabet_Filter_macro("*", "ABCDEFGHI-KLMN-PQRSTUVWXY-",
                               "ABCDEFGHI-KLMN-PQRSTUVWXY-"),
    /* Alphabet_Filter_TO_LOWER: */
    Alphabet_Filter_macro("*", "abcdefghi-klmn-pqrstuvwxy-",
                               "abcdefghi-klmn-pqrstuvwxy-"),
    /* Alphabet_Filter_COMPLEMENT: */
    Alphabet_Filter_macro("-", "TVGH--CD--M-KN---YSA-BW-R-",
                               "tvgh--cd--m-kn---ysa-bw-r-"),
    /* Alphabet_Filter_SOFTMASK_DNA: */
    Alphabet_Filter_macro("-", "ABCD--GH--K-MN---RST-VW-Y-",
                               "NNNN--NN--N-NN---NNN-NN-N-"),
    /* Alphabet_Filter_SOFTMASK_PROTEIN: */
    Alphabet_Filter_macro("*", "ABCDEFGHI-KLMN-PQRSTUVWXYZ",
                               "XXXXXXXXX-XXXX-XXXXXXXXXXX"),
    /* Alphabet_Filter_DNA_NOAMBIG: */
    Alphabet_Filter_macro("-", "A-C---G------------T------",
                               "A-C---G------------T------"),
    /* Alphabet_Filter_PROTEIN_NOAMBIG_SOFTMASK: */
    Alphabet_Filter_macro("*", "A-CDEFGHI-KLMN-PQRSTUVW-Y-",
                               "--------------------------"),
    /* Alphabet_Filter_PROTEIN_NOAMBIG: */
    Alphabet_Filter_macro("*", "A-CDEFGHI-KLMN-PQRSTUVW-Y-",
                               "A-CDEFGHI-KLMN-PQRSTUVW-Y-"),
    /* Alphabet_Filter_DNA_NOAMBIG_SOFTMASK: */
    Alphabet_Filter_macro("-", "A-C---G------------T------",
                               "--------------------------")
    };
    /* The protein alphabets include U for selenocysteine */
    g_assert(type >= 0);
    g_assert(type < Alphabet_Filter_NUMBER_OF_FILTERS);
    return filter_set[type];
    }

/**/

Alphabet *Alphabet_create(Alphabet_Type type, gboolean is_soft_masked){
    register Alphabet *a = g_new0(Alphabet, 1);
    a->ref_count = 1;
    a->type = type;
    a->is_soft_masked = is_soft_masked;
    /* Set look-up tables */
    a->unmasked = Alphabet_Filter_get(Alphabet_Filter_TO_UPPER);
    if(is_soft_masked){
        if(type == Alphabet_Type_PROTEIN){
            a->masked = Alphabet_Filter_get
                       (Alphabet_Filter_SOFTMASK_PROTEIN);
            a->non_ambig = Alphabet_Filter_get
                       (Alphabet_Filter_PROTEIN_NOAMBIG_SOFTMASK);
        } else {
            a->masked = Alphabet_Filter_get
                       (Alphabet_Filter_SOFTMASK_DNA);
            a->non_ambig = Alphabet_Filter_get
                       (Alphabet_Filter_DNA_NOAMBIG_SOFTMASK);
            }
    } else {
        a->masked = Alphabet_Filter_get(Alphabet_Filter_TO_UPPER);
        if(type == Alphabet_Type_PROTEIN){
            a->non_ambig = Alphabet_Filter_get
                (Alphabet_Filter_PROTEIN_NOAMBIG);
        } else {
            a->non_ambig = Alphabet_Filter_get
                (Alphabet_Filter_DNA_NOAMBIG);
            }
        }
    switch(type){
        case Alphabet_Type_DNA:
            a->member = (guchar*)"ACGT";
            a->is_valid
                = Alphabet_Filter_get(Alphabet_Filter_VALID_DNA_IUPAC);
            a->complement
                = Alphabet_Filter_get(Alphabet_Filter_COMPLEMENT);
            a->clean
                = Alphabet_Filter_get(Alphabet_Filter_CLEAN_DNA_IUPAC);
            break;
        case Alphabet_Type_PROTEIN:
            a->member = (guchar*)"ARNDCQEGHILKMFPSTWYUV*";
            a->is_valid
                = Alphabet_Filter_get(Alphabet_Filter_VALID_PROTEIN);
            a->complement = NULL;
            a->clean
                = Alphabet_Filter_get(Alphabet_Filter_CLEAN_PROTEIN);
            break;
        case Alphabet_Type_UNKNOWN:
            a->member = NULL;
            a->is_valid
                = Alphabet_Filter_get(Alphabet_Filter_VALID_SEQ_SYMBOL);
            a->complement = NULL;
            a->clean = NULL; /* Cannot clean unknown sequence type */
            break;
        default:
            g_error("Unknown alphabet type [%d]", type);
            break;
        }
    return a;
    }

void Alphabet_destroy(Alphabet *a){
    g_assert(a);
    if(--a->ref_count)
        return;
    g_free(a);
    return;
    }

Alphabet *Alphabet_share(Alphabet *a){
    g_assert(a);
    g_assert(a->ref_count);
    a->ref_count++;
    return a;
    }

/**/

guchar *Alphabet_get_filter_by_type(Alphabet *alphabet,
                                    Alphabet_Filter_Type filter_type){
    register guchar *filter = NULL;
    g_assert(filter_type >= 0);
    g_assert(filter_type < Alphabet_Filter_Type_NUMBER_OF_TYPES);
    switch(filter_type){
        case Alphabet_Filter_Type_UNMASKED:
            filter = alphabet->unmasked;
            break;
        case Alphabet_Filter_Type_MASKED:
            filter = alphabet->masked;
            break;
        case Alphabet_Filter_Type_IS_VALID:
            filter = alphabet->is_valid;
            break;
        case Alphabet_Filter_Type_COMPLEMENT:
            filter = alphabet->complement;
            break;
        case Alphabet_Filter_Type_CLEAN:
            filter = alphabet->clean;
            break;
        case Alphabet_Filter_Type_CLEAN_ACGTN:
            filter = Alphabet_Filter_get(
                     Alphabet_Filter_CLEAN_DNA_ACGTN);
            break;
        case Alphabet_Filter_Type_NON_AMBIG:
            filter = alphabet->non_ambig;
            break;
        default:
            g_error("Alphabet_Filter_Type unknown [%d]", filter_type);
            break;
        }
    return filter;
    }

gchar *Alphabet_Filter_Type_get_name(Alphabet_Filter_Type filter_type){
    register gchar *name = NULL;
    g_assert(filter_type >= 0);
    g_assert(filter_type < Alphabet_Filter_Type_NUMBER_OF_TYPES);
    switch(filter_type){
        case Alphabet_Filter_Type_UNMASKED:
            name = "unmasked";
            break;
        case Alphabet_Filter_Type_MASKED:
            name = "masked";
            break;
        case Alphabet_Filter_Type_IS_VALID:
            name = "valid";
            break;
        case Alphabet_Filter_Type_COMPLEMENT:
            name = "complement";
            break;
        case Alphabet_Filter_Type_CLEAN:
            name = "clean";
            break;
        case Alphabet_Filter_Type_CLEAN_ACGTN:
            name = "clean_acgtn";
            break;
        case Alphabet_Filter_Type_NON_AMBIG:
            name = "non_ambig";
            break;
        default:
            g_error("Alphabet_Filter_Type unknown [%d]", filter_type);
            break;
        }
    return name;
    }

/**/

gchar *Alphabet_Type_get_name(Alphabet_Type type){
    switch(type){
        case Alphabet_Type_DNA:
            return "DNA";
        case Alphabet_Type_PROTEIN:
            return "Protein";
        case Alphabet_Type_UNKNOWN:
            return "Unknown";
        }
    g_error("Unknown Alphabet Type [%d]", type);
    return NULL;
    }

Alphabet_Type Alphabet_name_get_type(gchar *name){
    register gint i;
    gchar *dna_string[5] = {"d", "dna", "n", "nt", "nucleotide"},
          *protein_string[4] = {"p", "protein", "aa", "aminoacid"},
          *unknown_string[3] = {"?", "unknown", "-"};
    for(i = 0; i < 5; i++)
        if(!strcasecmp(name, dna_string[i]))
            return Alphabet_Type_DNA;
    for(i = 0; i < 4; i++)
        if(!strcasecmp(name, protein_string[i]))
            return Alphabet_Type_PROTEIN;
    for(i = 0; i < 3; i++)
        if(!strcasecmp(name, unknown_string[i]))
            return Alphabet_Type_UNKNOWN;
    g_error("Unknown sequence type [%s]", name);
    return Alphabet_Type_UNKNOWN; /* not reached */
    }

gchar *Alphabet_aa2tla(gchar aa){
    register Alphabet_ArgumentSet *aas
           = Alphabet_ArgumentSet_create(NULL);
    static gchar *tla_names[25] = {
        "Ala", "Arg", "Asn", "Asp", "Cys", "Gln",
        "Glu", "Gly", "His", "Ile", "Leu", "Lys",
        "Met", "Phe", "Pro", "Ser", "Thr", "Trp",
        "Tyr", "Val", "Asx", "Zed", "Unk", "***",
        "Sec"
         };
    static gchar *short_names[25] = {
        "^A^", "^R^", "^N^", "^D^", "^C^", "^Q^",
        "^E^", "^G^", "^H^", "^I^", "^L^", "^K^",
        "^M^", "^F^", "^P^", "^S^", "^T^", "^W^",
        "^Y^", "^V^", "^B^", "^Z^", "^X^", "^*^",
        "^U^"};
    static guchar index[(1<<8)] = /* ARNDCQEGHILKMFPSTWYV */
     { 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25,
       25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25,
/* 32 */
       25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 23, 25, 25, 25, 25, 25,
       25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25,
/* 64 */
       25,  0, 20,  4,  3,  6, 13,  7,  8,  9, 25, 11, 10, 12,  2, 25,
       14,  5,  1, 15, 16, 25, 19, 17, 22, 18, 21, 25, 25, 25, 25, 25,
       25,  0, 20,  4,  3,  6, 13,  7,  8,  9, 25, 11, 10, 12,  2, 25,
       14,  5,  1, 15, 16, 25, 19, 17, 22, 18, 21, 25, 25, 25, 25, 25,
/* 128 */
       25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25,
       25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25,
       25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25,
       25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25,
/* 192 */
       25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25,
       25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25,
       25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25,
       25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25
/* 256 */
       };
    if(index[(guchar)aa] == 25)
        g_error("Unknown amino acid [%c]", aa);
    g_assert(short_names[index[(guchar)aa]][1] == toupper(aa));
    if(aas->use_aa_tla)
        return tla_names[index[(guchar)aa]];
    return short_names[index[(guchar)aa]];
    }

gchar *Alphabet_nt2ambig(gchar nt){
    static gchar *ambig[16] = {
        /* - */ NULL,
        /* G */ "G",
        /* A */ "A",
        /* R */ "AG",
        /* T */ "T",
        /* K */ "GT",
        /* W */ "AT",
        /* D */ "AGT",
        /* C */ "C",
        /* S */ "CG",
        /* M */ "AC",
        /* V */ "ACG",
        /* Y */ "CT",
        /* B */ "CGT",
        /* H */ "ACT",
        /* N */ "ACGT"
         };
    static guchar index[(1<<8)] = { /* GARTKWDCSMVYBHN */
         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
         0,  2, 13,  8,  7,  0,  0,  1, 14,  0,  0,  5,  0, 10, 15,  0,
         0,  0,  3,  9,  4,  0, 11,  6,  0, 12,  0,  0,  0,  0,  0,  0,
         0,  2, 13,  8,  7,  0,  0,  1, 14,  0,  0,  5,  0, 10, 15,  0,
         0,  0,  3,  9,  4,  0, 11,  6,  0, 12,  0,  0,  0,  0,  0,  0,
         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        };
    g_assert(index[(guchar)nt]);
    return ambig[index[(guchar)nt]];
    }

