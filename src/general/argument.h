/****************************************************************\
*                                                                *
*  Library for command line argument processing                  *
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

#ifndef INCLUDED_ARGUMENT_H
#define INCLUDED_ARGUMENT_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <stdio.h>
#include <limits.h> /* For CHAR_BIT */

#include <glib.h>

#ifndef ALPHABETSIZE
#define ALPHABETSIZE (1<<CHAR_BIT)
#endif /* ALPHABETSIZE */

typedef gchar *(*ArgumentHandler)(gchar *arg_string,
                                  gpointer data);

typedef void (*Argument_Cleanup_Func)(gpointer user_data);

typedef struct {
        gchar *desc;       /* Description */
    GPtrArray *arg_option; /* Sets of argument options */
} ArgumentSet;

typedef enum {
    ArgumentSource_NOT_SET = 0,
    ArgumentSource_FROM_LONG_ARGUMENT,
    ArgumentSource_FROM_SHORT_ARGUMENT,
    ArgumentSource_FROM_ENV_VAR
} ArgumentSource;

typedef struct {
        ArgumentSet *as;
              gchar  symbol;
              gchar *option;
              gchar *type; /* NULL for booleans */
              gchar *desc;
              gchar *default_string; /* NULL for mandatory arguments */
    ArgumentHandler  handler; /* NULL for lists */
           gpointer  handler_data;
              gchar *env_var;
     ArgumentSource  source;
              gchar *arg_value;
          GPtrArray *arg_list; /* Only used when a list argument */
} ArgumentOption;

typedef struct {
              gint   argc;
             gchar **argv;
             gchar  *name; /* Program name */
             gchar  *desc; /* Description  */
         GPtrArray  *arg_set; /* Set of argument sets */
         GPtrArray  *mandatory_set; /* Set of mandatory arguments */
    ArgumentOption  *symbol_registry[ALPHABETSIZE];
              void  *option_registry; /* root node of binary search tree */
          gboolean   show_short_help;
          gboolean   show_long_help;
         GPtrArray  *cleanup_list;
} Argument;

gint Argument_main(Argument *arg);
void Argument_process(Argument *arg, gchar *name, gchar *desc,
                      gchar *synopsis);
void Argument_info(Argument *arg);

ArgumentSet *ArgumentSet_create(gchar *desc);
       void  Argument_absorb_ArgumentSet(Argument *arg,
                                 ArgumentSet *as);
/* ArgumentSet does not need to be freed after absorption */

void ArgumentSet_add_option(ArgumentSet *as,
              gchar symbol,
              gchar *option,
              gchar *type,
              gchar *desc,
              gchar *default_string,
              ArgumentHandler handler,
              gpointer handler_data);
/* Set <symbol> to '\0' if a single-letter flag is not required.
 *
 * If no handler is supplied, the argument is assumed to be a list.
 */

/* Some common handlers */
gchar *Argument_parse_string(gchar *arg_string, gpointer data);
gchar *Argument_parse_char(gchar *arg_string, gpointer data);
gchar *Argument_parse_int(gchar *arg_string, gpointer data);
gchar *Argument_parse_float(gchar *arg_string, gpointer data);
gchar *Argument_parse_boolean(gchar *arg_string, gpointer data);

void Argument_add_cleanup(Argument *arg,
                          Argument_Cleanup_Func cleanup_func,
                          gpointer user_data);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* INCLUDED_ARGUMENT_H */

