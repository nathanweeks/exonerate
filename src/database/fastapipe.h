/****************************************************************\
*                                                                *
*  Library for All-vs-All Fasta Database Comparisons             *
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

#ifndef INCLUDED_FASTAPIPE_H
#define INCLUDED_FASTAPIPE_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <glib.h>

#include "fastadb.h"

typedef gboolean (*FastaPipe_NextSeq_Func)(FastaDB_Seq *fdbs,
                                           gpointer user_data);
typedef void (*FastaPipe_Boundary_Func)(gpointer user_data);

typedef struct {
                   FastaDB *query_db;
                   FastaDB *target_db;
              FastaDB_Mask  mask;
   FastaPipe_Boundary_Func  init_func; /* pre_query_load */
   FastaPipe_Boundary_Func  prep_func; /* post_query_load */
   FastaPipe_Boundary_Func  term_func; /* post_target_pass */
    FastaPipe_NextSeq_Func  query_func;
    FastaPipe_NextSeq_Func  target_func;
                  gboolean  is_full;
                  gboolean  revcomp_query;
                  gboolean  revcomp_target;
               FastaDB_Seq *prev_query;
               FastaDB_Seq *prev_target;
} FastaPipe;

FastaPipe *FastaPipe_create(FastaDB *query_fdb, FastaDB *target_fdb,
                            FastaPipe_Boundary_Func init_func,
                            FastaPipe_Boundary_Func prep_func,
                            FastaPipe_Boundary_Func term_func,
                            FastaPipe_NextSeq_Func query_func,
                            FastaPipe_NextSeq_Func target_func,
                            FastaDB_Mask mask, gboolean translate_both,
                            gboolean use_revcomp);
/* query_func can return TRUE to finish loading current pipeline
 * target_func can return TRUE to finish processing of the pipeline
 */

/* init_func is called before query pipeline loading
 * prep_func is called at the end of query pipeline loading
 * term_func is called at the end of query pipeline analysis
 */

     void  FastaPipe_destroy(FastaPipe *fasta_pipe);
 gboolean  FastaPipe_process(FastaPipe *fasta_pipe, gpointer user_data);

/* FastaPipe_process()
 * will call query_func until end of db, or TRUE is returned,
 * then calls target_func until end of db, or TRUE is returned.
 * If not at end of target_db, will rewind query_func and continue.
 * Returns TRUE when true is returned from the target_func
 */

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* INCLUDED_FASTAPIPE_H */

