/****************************************************************\
*                                                                *
*  Nucleotide Translation Code                                   *
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

#include <ctype.h>
#include <string.h>
#include <stdlib.h> /* For atoi() */

#include "translate.h"

/* ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
     CTAG
 0 - 0000 - blank          Representation of Nucleotides
 1 G 0001 C                -----------------------------
 2 A 0010 T
 3 R 0011 Y (AG)       The standard (IUB) nucleotides
 4 T 0100 A            are used in the array "-GARTKWDCSMVYBHN"
 5 K 0101 M (GT)       so that:
 6 W 0110 W (AT)
 7 D 0111 H (AGT)       1: Four bits per base are used.
 8 C 1000 G             2: A bit is set for each base represented.
 9 S 1001 S (CG)        3: Any base reversed is it's complement.
10 M 1010 K (AC)
11 V 1011 B (ACG)
12 Y 1100 R (CT)
13 B 1101 V (CGT)
14 H 1110 D (ACT)
15 N 1111 N (ATGC)
 :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: */

/**/

Translate_ArgumentSet *Translate_ArgumentSet_create(Argument *arg){
    register ArgumentSet *as;
    static Translate_ArgumentSet tas = {NULL};
    if(arg){
        as = ArgumentSet_create("Translation Options");
        ArgumentSet_add_option(as, '\0', "geneticcode", NULL,
            "Use built-in or custom genetic code", "1",
            Argument_parse_string, &tas.genetic_code);
        Argument_absorb_ArgumentSet(arg, as);
        }
    return &tas;
    }

/**/

static void Translate_initialise_nucleotide_data(Translate *t){
    register gint i;
    for(i = 0; i < Translate_NT_SET_SIZE; i++)
        t->nt2d[t->nt[i]] = t->nt2d[tolower(t->nt[i])] = i;
    t->nt2d['X'] = t->nt2d['x'] = t->nt2d['N'];
    t->nt2d['U'] = t->nt2d['u'] = t->nt2d['T'];
    return;
    }

static void Translate_initialise_peptide_data(Translate *t){
    register gint i, j;
    guchar pimagrp[Translate_PIMA_SET_SIZE][6] = {
     "aIV", "bLM",  "dFWY", "lND",   "kDE",  "oEQ",
     "nKR", "iST",  "hAG",  "cab",   "edH",  "mlk",
     "pon", "jihP", "fCcd", "rHmpi", "xfrj", "Xx*" };
    for(i = 0; i < Translate_AA_SET_SIZE; i++)
        t->aa2d[t->aa[i]] = i;
    for(i = 1; i < 23; i++)
        t->aamask[i] = (1L<<(i-1));  /* First is zero  */
    for(i = 0; i < Translate_PIMA_SET_SIZE; i++){
        t->aamask[t->aa2d[pimagrp[i][0]]]
      = t->aamask[t->aa2d[pimagrp[i][1]]];
        for(j = 2; pimagrp[i][j]; j++)
            t->aamask[t->aa2d[pimagrp[i][0]]]
         |= t->aamask[t->aa2d[pimagrp[i][j]]];
        }
    return;
    }

static void Translate_initialise_translation_data(Translate *t){
    register gchar a, b, c, x, y, z;
    register gint i, tmp;
    for(x = 0; x < 16; x++){
        for(y = 0; y < 16; y++){
            for(z = 0; z < 16; z++){
                tmp = 0;
                for(a = 0; a < 4; a++){
                    if(x != (x|1<<a))
                        continue;
                    for(b = 0; b < 4; b++){
                        if(y != (y|1<<b))
                            continue;
                        for(c = 0; c < 4; c++){
                            if(z != (z|1<<c))
                                continue;
                            tmp = (tmp|t->aamask[
                                t->aa2d[t->code[((a<<4)|(b<<2)|c)]]]);
                            }
                        }
                    }
                for(i = 0; t->aamask[i] != (tmp|t->aamask[i]); i++);
                t->trans[x|(y<<4)|(z<<8)] = i;
                }
            }
       }
    return;
    }

static void Translate_initialise_reverse_translate_data(Translate *t){
    register gint i;
    for(i = 0; t->code[i]; i++){
        if(!t->revtrans[t->code[i]])
            t->revtrans[t->code[i]] = g_ptr_array_new();
        g_ptr_array_add(t->revtrans[t->code[i]], GINT_TO_POINTER(i));
        }
    return;
    }

/**/

static gchar *Translate_convert_genetic_code(gchar *code){
    register gint a, b, c, n = 0;
    gint table[4] = {3, 2, 0, 1};
    register gchar *result;
    g_assert(code && (strlen(code) == 64));
    result = g_new(gchar, 65);
    result[64] = '\0';
    for(a = 0; a < 4; a++)
        for(b = 0; b < 4; b++)
            for(c = 0; c < 4; c++)
                result[n++] = code[(table[a] << 4)
                                 | (table[b] << 2)
                                 |  table[c]];
    return result;
    }
/* Converts genetic code from NCBI format (TCAG order)
 * to the format internally by exonerate (GATC order)
 */

typedef struct {
     gint  id;
    gchar *name;
    gchar *code;
} Translate_GeneticCode;

static gchar *Translate_get_builtin_genetic_code(gint id){
    register gint i;
    Translate_GeneticCode code[17] = {
        { 1, "The Standard Code",
        "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"},
        { 2, "The Vertebrate Mitochondrial Code",
        "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG"},
        { 3, "The Yeast Mitochondrial Code",
        "FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG"},
        { 4, "The Mold, Protozoan, and Coelenterate Mitochondrial Code"
             " and the Mycoplasma/Spiroplasma Code",
        "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"},
        { 5, "The Invertebrate Mitochondrial Code",
        "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG"},
        { 6, "The Ciliate, Dasycladacean and Hexamita Nuclear Code",
        "FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"},
        { 9, "The Echinoderm and Flatworm Mitochondrial Code",
        "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG"},
        {10, "The Euplotid Nuclear Code",
        "FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"},
        {11, "The Bacterial and Plant Plastid Code",
        "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"},
        {12, "The Alternative Yeast Nuclear Code",
        "FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"},
        {13, "The Ascidian Mitochondrial Code",
        "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG"},
        {14, "The Alternative Flatworm Mitochondrial Code",
        "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG"},
        {15, "Blepharisma Nuclear Code",
        "FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"},
        {16, "Chlorophycean Mitochondrial Code",
        "FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"},
        {21, "Trematode Mitochondrial Code",
        "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG"},
        {22, "Scenedesmus obliquus mitochondrial Code",
        "FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"},
        {23, "Thraustochytrium Mitochondrial Code",
        "FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"}
        };
    for(i = 0; i < 17; i++){
        if(code[i].id == id)
            return Translate_convert_genetic_code(code[i].code);
        }
    g_error("No built in genetic code corresponding to id [%d]", id);
    return NULL;
    }
/* Built in genetic codes taken from:
 * http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
 */

static gchar *Translate_get_genetic_code(gchar *code){
    register gint len, id;
    if(!code)
        return g_strdup(
         "GGGGEEDDVVVVAAAARRSSKKNNMIIITTTTW*CC**YYLLFFSSSSRRRRQQHHLLLLPPPP");
    len = strlen(code);
    if(len == 64)
        return Translate_convert_genetic_code(code);
    if((len == 1) || (len == 2)){
        id = atoi(code);
        return Translate_get_builtin_genetic_code(id);
        }
    g_error("Could not use genetic code [%s]", code);
    return NULL;
    }

/**/

Translate *Translate_create(gboolean use_pima){
    register Translate *t = g_new0(Translate, 1);
    register Translate_ArgumentSet *tas = Translate_ArgumentSet_create(NULL);
    t->ref_count = 1;
    t->nt = (guchar*)"-GARTKWDCSMVYBHN";
    t->aa = (guchar*)"-ARNDCQEGHILKMFPSTWYV*ablkonihdmcepjfrxX";
    t->code = (guchar*)Translate_get_genetic_code(tas->genetic_code);
    Translate_initialise_nucleotide_data(t);
    Translate_initialise_peptide_data(t);
    Translate_initialise_translation_data(t);
    Translate_initialise_reverse_translate_data(t);
    if(!use_pima)
        t->aa = (guchar*)"-ARNDCQEGHILKMFPSTWYV*XXXXXXXXXXXXXXXXXX";
    return t;
    }

Translate *Translate_share(Translate *t){
    g_assert(t);
    t->ref_count++;
    return t;
    }

void Translate_destroy(Translate *t){
    register gint i;
    g_assert(t);
    if(--t->ref_count)
        return;
    for(i = 0; t->code[i]; i++){
        if(t->revtrans[t->code[i]]){
            g_ptr_array_free(t->revtrans[t->code[i]], TRUE);
            t->revtrans[t->code[i]] = NULL;
            }
        }
    g_free(t->code);
    g_free(t);
    return;
    }

gint Translate_sequence(Translate *t, gchar *dna, gint dna_length,
                        gint frame, gchar *aaseq, guchar *filter){
    register gchar *dp, *ap = aaseq, *end;
    if((frame > 0) && (frame < 4)){
        end = dna + dna_length-2;
        if(filter){
            for(dp = dna + frame - 1; dp < end; dp += 3)
                *ap++ = Translate_base(t, filter[(guchar)dp[0]],
                                          filter[(guchar)dp[1]],
                                          filter[(guchar)dp[2]]);
        } else {
            for(dp = dna + frame - 1; dp < end; dp += 3)
                *ap++ = Translate_codon(t, dp);
            }
        *ap = '\0';
        return ap - aaseq;
        }
    if((frame < 0) && (frame > -4)){
        if(filter){
            for(dp = dna + dna_length + frame - 2; dp >= dna; dp -= 3)
                *ap++ = Translate_base(t, filter[(guchar)dp[0]],
                                          filter[(guchar)dp[1]],
                                          filter[(guchar)dp[2]]);
        } else {
            for(dp = dna + dna_length + frame - 2; dp >= dna; dp -= 3)
                *ap++ = Translate_codon(t, dp);
            }
        *ap = '\0';
        return ap - aaseq;
        }
    g_error("Invalid reading frame [%d]", frame);
    return 0; /* Never reached */
    }
/* Returns length of peptide generated
 */

void Translate_reverse(Translate *t, gchar *aaseq, gint length,
                       Translate_reverse_func trf, gpointer user_data){
    register gint i, j, total, *prod, id, codon;
    register gchar *nucleotide;
    g_assert(length > 0);
    g_assert(aaseq);
    g_assert(aaseq[0]);
    g_assert(trf);
    prod = g_new0(gint, length);
    nucleotide = g_new(gchar, (length*3)+1);
    nucleotide[length*3] = '\0';
    total = t->revtrans[(guchar)aaseq[0]]
          ? t->revtrans[(guchar)aaseq[0]]->len
          : 1;
    for(i = 1; i < length; i++){
        if(t->revtrans[(guchar)aaseq[i]])
            total *= t->revtrans[(guchar)aaseq[i]]->len;
        }
    prod[0] = 1;
    for(i = 1; i < length; i++)
        prod[i] = prod[i-1]
                * (t->revtrans[(guchar)aaseq[i-1]]
                ? t->revtrans[(guchar)aaseq[i-1]]->len
                : 1);
    for(i = 0; i < total; i++){
        for(j = 0; j < length; j++){
            if(t->revtrans[(guchar)aaseq[j]]){
                id = (i/prod[j])%t->revtrans[(guchar)aaseq[j]]->len;
                codon = GPOINTER_TO_INT(
                        t->revtrans[(guchar)aaseq[j]]->pdata[id]);
                nucleotide[j*3]     = t->nt[1<<((codon>>4)&3)];
                nucleotide[(j*3)+1] = t->nt[1<<((codon>>2)&3)];
                nucleotide[(j*3)+2] = t->nt[1<<(codon&3)];
            } else {
                nucleotide[j*3]     =
                nucleotide[(j*3)+1] =
                nucleotide[(j*3)+2] = 'N';
                }
            }
        trf(nucleotide, length*3, user_data);
        }
    g_free(nucleotide);
    g_free(prod);
    return;
    }

/* TODO: Add alternative genetic codes (when required) from:
 * http://www.ncbi.nlm.nih.gov/htbin-post/Taxonomy/wprintgc?mode=c
 */

