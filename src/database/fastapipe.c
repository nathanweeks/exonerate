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

#include "fastapipe.h"

FastaPipe *FastaPipe_create(FastaDB *query_fdb, FastaDB *target_fdb,
               FastaPipe_Boundary_Func init_func,
               FastaPipe_Boundary_Func prep_func,
               FastaPipe_Boundary_Func term_func,
               FastaPipe_NextSeq_Func query_func,
               FastaPipe_NextSeq_Func target_func,
               FastaDB_Mask mask, gboolean translate_both,
               gboolean use_revcomp){
    register FastaPipe *fasta_pipe = g_new(FastaPipe, 1);
    g_assert(init_func);
    g_assert(prep_func);
    g_assert(term_func);
    g_assert(query_func);
    g_assert(target_func);
    fasta_pipe->query_db = FastaDB_share(query_fdb);
    fasta_pipe->target_db = FastaDB_share(target_fdb);
    fasta_pipe->init_func = init_func;
    fasta_pipe->prep_func = prep_func;
    fasta_pipe->term_func = term_func;
    fasta_pipe->query_func = query_func;
    fasta_pipe->target_func = target_func;
    fasta_pipe->mask = mask;
    fasta_pipe->is_full = FALSE;
    if(use_revcomp){
        fasta_pipe->revcomp_query
            = (query_fdb->alphabet->type == Alphabet_Type_DNA);
        fasta_pipe->revcomp_target
            = (((query_fdb->alphabet->type == Alphabet_Type_PROTEIN)
                && (target_fdb->alphabet->type == Alphabet_Type_DNA))
              || translate_both);
    } else {
        fasta_pipe->revcomp_query = FALSE;
        fasta_pipe->revcomp_target = FALSE;
        }
    fasta_pipe->prev_query = NULL;
    fasta_pipe->prev_target = NULL;
    return fasta_pipe;
    }

void FastaPipe_destroy(FastaPipe *fasta_pipe){
    FastaDB_close(fasta_pipe->query_db);
    FastaDB_close(fasta_pipe->target_db);
    g_free(fasta_pipe);
    return;
    }

static FastaDB_Seq *FastaPipe_next_seq(FastaDB *fdb, FastaDB_Mask mask,
                    gboolean use_revcomp, FastaDB_Seq **prev_seq){
    register FastaDB_Seq *fdbs, *revcomp_fdbs;
    if(*prev_seq){
        fdbs = FastaDB_Seq_revcomp(*prev_seq);
        FastaDB_Seq_destroy(*prev_seq);
        (*prev_seq) = NULL;
    } else {
        fdbs = FastaDB_next(fdb, mask);
        if(use_revcomp && fdbs){
            if((fdbs->seq->annotation)
            && (fdbs->seq->annotation->strand == Sequence_Strand_REVCOMP)){
                revcomp_fdbs = FastaDB_Seq_revcomp(fdbs);
                FastaDB_Seq_destroy(fdbs);
                return revcomp_fdbs;
                }
            if((!fdbs->seq->annotation)
             || (fdbs->seq->annotation->strand == Sequence_Strand_REVCOMP))
                (*prev_seq) = FastaDB_Seq_share(fdbs);
            }
        }
    return fdbs;
    }

static gboolean FastaPipe_process_database(FastaDB *fdb,
                          FastaPipe_NextSeq_Func next_seq_func,
                          FastaDB_Mask mask, gboolean use_revcomp,
                          FastaDB_Seq **prev_seq, gpointer user_data){
    register FastaDB_Seq *fdbs;
    register gboolean stop_requested = FALSE;
    while((fdbs = FastaPipe_next_seq(fdb, mask,
                                     use_revcomp, prev_seq))){
        stop_requested = next_seq_func(fdbs, user_data);
        FastaDB_Seq_destroy(fdbs);
        if(stop_requested)
            break;
        }
    return stop_requested;
    }
/* Will return FALSE at end of database file
 */

gboolean FastaPipe_process(FastaPipe *fasta_pipe, gpointer user_data){
    register gboolean target_stop = FALSE;
    do {
        if(!fasta_pipe->is_full){
            if(FastaDB_is_finished(fasta_pipe->query_db)
            && (!fasta_pipe->prev_query))
                break;
            fasta_pipe->init_func(user_data);
            FastaPipe_process_database(fasta_pipe->query_db,
                                       fasta_pipe->query_func,
                                       fasta_pipe->mask,
                                       fasta_pipe->revcomp_query,
                                       &fasta_pipe->prev_query,
                                       user_data);
            fasta_pipe->is_full = TRUE;
            fasta_pipe->prep_func(user_data);
            }
        target_stop = FastaPipe_process_database(
                                   fasta_pipe->target_db,
                                   fasta_pipe->target_func,
                                   fasta_pipe->mask,
                                   fasta_pipe->revcomp_target,
                                   &fasta_pipe->prev_target,
                                   user_data);
        if(!target_stop){
            FastaDB_rewind(fasta_pipe->target_db);
            fasta_pipe->is_full = FALSE;
            fasta_pipe->term_func(user_data);
            }
    } while(!target_stop);
    return target_stop;
    }

