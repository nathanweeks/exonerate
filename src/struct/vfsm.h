/****************************************************************\
*                                                                *
*  VFSM Library : complete trie navigation                       *
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

#ifndef INCLUDED_VFSM_H
#define INCLUDED_VFSM_H

#include <glib.h>
#include <limits.h> /* For CHAR_BIT */

#ifdef G_HAVE_GINT64
typedef guint64 VFSM_Int;
#define VFSM_MAX G_GINT64_CONSTANT(18446744073709551615U)
#ifdef CUSTOM_GUINT64_FORMAT
#define VFSM_PRINT_FORMAT CUSTOM_GUINT64_FORMAT
#else /* CUSTOM_GUINT64_FORMAT */
#define VFSM_PRINT_FORMAT "llu"
#endif /* CUSTOM_GUINT64_FORMAT */
#else /* G_HAVE_GINT64 */
typedef guint32 VFSM_Int;
#define VFSM_MAX 4294967295U
#define VFSM_PRINT_FORMAT "u"
#endif /* G_HAVE_GINT64 */

#ifndef ALPHABETSIZE
#define ALPHABETSIZE (1<<(CHAR_BIT))
#endif /* ALPHABETSIZE */

typedef struct {
        gchar *alphabet;
        guint  alphabet_size;
        gchar  index[ALPHABETSIZE]; /* alphabet index             */
        guint  depth;               /* trie depth or wordsize     */
     gboolean  is_poweroftwo;
        guint  log_alphabet_size;   /* Used for &(base-1) tricks  */
     VFSM_Int  prs;                 /* Penultimate row start      */
     VFSM_Int  prw;                 /* Penultimate row width      */
     VFSM_Int  lrs;                 /* Last row start             */
     VFSM_Int  lrw;                 /* Last row width (leafcount) */
     VFSM_Int  total;               /* Total number of states     */
     VFSM_Int *branch_size;         /* For each depth position    */
} VFSM;

VFSM *VFSM_create(gchar *alphabet, guint depth);
/* Returns NULL if cannot build VFSM using basic types */

void VFSM_destroy(VFSM *vfsm);
void VFSM_info(VFSM *vfsm);

VFSM_Int VFSM_word2state(VFSM *vfsm, gchar *word);
/* Returns 0 if word is not valid */

gboolean VFSM_state2word(VFSM *vfsm, VFSM_Int state, gchar *word);
/* Returns FALSE if pos is not valid */
/* Requires vfsm->depth+1 char space for word */

VFSM_Int VFSM_change_state(VFSM *vfsm, VFSM_Int state, guchar ch);
VFSM_Int VFSM_change_state_POW2(VFSM *vfsm, VFSM_Int state, guchar ch);

#define VFSM_state_is_leaf(vfsm, state) ((state) >= (vfsm)->lrs)
#define VFSM_state2leaf(vfsm, state) ((state)-(vfsm)->lrs)
#define VFSM_leaf2state(vfsm, state) ((state)+(vfsm)->lrs)

#define VFSM_change_state_M(vfsm, state, ch)               \
    ((state) < (vfsm)->lrs)                                \
    ? ((state)*(vfsm)->alphabet_size)+(vfsm)->index[(ch)]  \
    : (((vfsm)->prs+(((state)-(vfsm)->lrs) % (vfsm)->prw)) \
    * (vfsm)->alphabet_size)+(vfsm)->index[(ch)];

#define VFSM_change_state_MPOW2(vfsm, state, ch)                \
    ((state) < (vfsm)->lrs)                                     \
    ? ((state)<<(vfsm)->log_alphabet_size)+(vfsm)->index[(ch)]  \
    : (((vfsm)->prs+(((state)-(vfsm)->lrs)&((vfsm)->prw-1)))    \
      << (vfsm)->log_alphabet_size)+(vfsm)->index[(ch)];

gchar VFSM_symbol_by_pos(VFSM *vfsm, VFSM_Int state, guint depth);
gchar VFSM_symbol_by_pos_POW2(VFSM *vfsm, VFSM_Int state, guint depth);

VFSM_Int VFSM_jump(VFSM *vfsm, VFSM_Int state,
                   guint depth, guchar next);

#define VFSM_jump_M(vfsm, state, depth, next, curr)  \
    ((state)+((vfsm)->branch_size[(depth)]           \
    *((vfsm)->index[(curr)]-(vfsm)->index[(next)])))

/* FIXME: some of these macros will have to change */

# endif /* INCLUDED_VFSM_H */

