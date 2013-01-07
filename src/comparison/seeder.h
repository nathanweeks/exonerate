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

#ifndef INCLUDED_SEEDER_H
#define INCLUDED_SEEDER_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <glib.h>
#include "argument.h"
#include "recyclebin.h"
#include "hspset.h"
#include "fsm.h"
#include "vfsm.h"
#include "comparison.h"

typedef struct {
       gsize  fsm_memory_limit;
       gchar *force_fsm;
        gint  word_jump;
        gint  word_ambiguity;
} Seeder_ArgumentSet;

Seeder_ArgumentSet *Seeder_ArgumentSet_create(Argument *arg);

/**/

typedef struct {
  HSP_Param *hsp_param;
       gint  tpos_modifier;
    gdouble  saturate_expectation;
     size_t  hspset_offset;
} Seeder_Loader;

typedef struct {
      Sequence *query;
    Comparison *curr_comparison;
} Seeder_QueryInfo;

typedef struct {
          Seeder_Loader *loader;
       Seeder_QueryInfo *query_info;
} Seeder_Context;

/**/

#if 0

typedef struct {
              gint  query_pos;
    Seeder_Context *context;
} Seeder_Seed;

typedef struct Seeder_Seed_List {
    struct Seeder_Seed_List *next;
                 SeederSeed  seed;
} Seeder_Seed_List;

typedef struct {
    union {
           Seeder_Seed *seed;
      CompoundFile_Pos *offset;
    } data;
                  gint  total; /* -ve when extmem */
} Seeder_Seed_Array;

typedef struct Seeder_Neighbour_Array {
    union {
     struct Seeder_WordInfo *word_info;
           CompoundFile_Pos *offset;
    } data;
                       gint  total; /* -ve when extmem */
} Seeder_Neighbour_Array;

typedef struct Seeder_Neighbour_List {
    struct Seeder_Neighbour_List *next;
          struct Seeder_WordInfo *word_info;
} Seeder_Neighbour_List ;

typedef struct Seeder_WordInfo {
    union {
          Seeder_Seed_List *list;
         Seeder_Seed_Array *array;
    } seed;
    union {
          Seeder_Neighbour_List *list;
         Seeder_Neighbour_Array *array;
    } neighbour; /* Absent w/o wordhood */
    gint  match_count;     /* Absent w/o satn thold      */
    gint  match_mailbox;   /* Absent w/o satn thold      */
} Seeder_WordInfo;

#endif /* 0 */

/**/

typedef struct Seeder_Seed { /* For each occurence of the word in a query */
    struct Seeder_Seed *next;
                  gint  query_pos;
        Seeder_Context *context;
} Seeder_Seed;

typedef struct Seeder_Neighbour { /* For each neighbour of the FSM word */
    struct Seeder_Neighbour *next;
     struct Seeder_WordInfo *word_info;
} Seeder_Neighbour;

typedef struct Seeder_WordInfo_no_ST {
         Seeder_Seed *seed_list;
    Seeder_Neighbour *neighbour_list;
} Seeder_WordInfo_no_ST;

typedef struct Seeder_WordInfo { /* For each word in the FSM */
         Seeder_Seed *seed_list;
    Seeder_Neighbour *neighbour_list; /* Absent w/o wordhood   */
                gint  match_count;    /* Absent w/o satn thold */
                gint  match_mailbox;  /* Absent w/o satn thold */
} Seeder_WordInfo;

/* FIXME: remove neighbour_list and match_{count,mailbox} when unnecessary
 */

/*  add compact list formats for {seed,neighbour}_list
 *  Seeder_Seed_Compact {SeederSeed *data, gint length}
 *  Seeder_Neighbour_Compact {SeederSeed *data, gint length}
 */

/**/

typedef struct {
    FSM *fsm;
} Seeder_FSM;

typedef struct {
               VFSM  *vfsm;
    Seeder_WordInfo **leaf;
} Seeder_VFSM;

/**/

typedef void (*Seeder_ReportFunc)(Comparison *comparison,
                                  gpointer user_data);

typedef struct {
    Seeder_ArgumentSet *sas;
              gboolean  is_prepared;
                  gint  verbosity;
     Seeder_ReportFunc  report_func;
              gpointer  user_data;
                  gint  total_query_length;
                  gint  comparison_count;
           Match_Score  saturate_threshold;
          Match_Strand *target_strand;
      Comparison_Param *comparison_param;
             GPtrArray *query_info_list;
            RecycleBin *recycle_wordinfo;
            RecycleBin *recycle_seed;
            RecycleBin *recycle_neighbour;
            RecycleBin *recycle_context;
            Seeder_FSM *seeder_fsm;  /* NULL when not used */
           Seeder_VFSM *seeder_vfsm; /* NULL when not used */
         Seeder_Loader *dna_loader;
         Seeder_Loader *protein_loader;
         Seeder_Loader *codon_loader;
             GPtrArray *active_queryinfo_list;
             HSP_Param *any_hsp_param;
} Seeder;

Seeder *Seeder_create(gint verbosity,
                      Comparison_Param *comparison_param,
                      Match_Score saturate_threshold,
                      Seeder_ReportFunc report_func,
                      gpointer user_data);
    void  Seeder_destroy(Seeder *seeder);
   gsize  Seeder_memory_usage(Seeder *seeder);
    void  Seeder_memory_info(Seeder *seeder);
gboolean  Seeder_add_query(Seeder *seeder, Sequence *query);
    void  Seeder_add_target(Seeder *seeder, Sequence *target);
/* Must add all queries before adding any targets */

/* FIXME: add Seeder_create_external();
 *        to add a dataset to an empty seeder
 */

/**/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* INCLUDED_SEEDER_H */

