/****************************************************************\
*                                                                *
*  Library for FSM-based word matching.                          *
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

#ifndef INCLUDED_FSM_H
#define INCLUDED_FSM_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <limits.h>
#include <glib.h> /* For CHAR_BIT */
#include "recyclebin.h"

#ifndef ALPHABETSIZE
#define ALPHABETSIZE (1<<(CHAR_BIT))
#endif /* ALPHABETSIZE */

typedef struct FSM_Node {
    struct FSM_Node *next;
           gpointer  data;
} FSM_Node;

typedef gpointer (*FSM_Join_Func)(gpointer a, gpointer b,
                                  gpointer user_data);

typedef struct {
         guchar  index[ALPHABETSIZE]; /* Map chars to the array      */
         guchar  insertion_filter[ALPHABETSIZE];
         guchar  traversal_filter[ALPHABETSIZE];
         guchar  width;        /* Alphabet size for this FSM         */
       FSM_Node *root;         /* Base of the finite state machine   */
  FSM_Join_Func  merge_func;   /* Function for node data merging     */
  FSM_Join_Func  combine_func; /* Function for node data combination */
     RecycleBin *recycle;      /* RecycleBin for node chunks         */
         gulong  chunk_count;  /* Number of chunks used              */
       gpointer  user_data;
       gboolean  is_compiled;
} FSM;
/* merge_func is called when two identical words are submitted.
 * combine_func is called when one word is a subseq of another.
 */

  FSM *FSM_create(gchar *alphabet, FSM_Join_Func merge_func,
                 FSM_Join_Func combine_func, gpointer user_data);
gsize  FSM_memory_usage(FSM *f);
 void  FSM_info(FSM *f);
 void  FSM_destroy(FSM *f);

typedef void (*FSM_Destroy_Func)(gpointer data, gpointer user_data);
void FSM_destroy_with_data(FSM *f, FSM_Destroy_Func fdf,
                                   gpointer user_data);

gpointer FSM_add(FSM *f, gchar *seq, guint len, gpointer node_data);
void FSM_compile(FSM *f);

typedef void (*FSM_Traverse_Func)(guint seq_pos,
                                  gpointer node_data,
                                  gpointer user_data);

void FSM_traverse(FSM *f, gchar *seq, FSM_Traverse_Func ftf,
                  gpointer user_data);

void FSM_add_insertion_filter(FSM *f, guchar *filter);
void FSM_add_traversal_filter(FSM *f, guchar *filter);
/* Filters must be ALPHABETSIZE.
 * If filter is NULL, mapping is reset.
 */

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* INCLUDED_FSM_H */

