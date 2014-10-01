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

#include "argument.h"
#include <search.h> /* For tdelete(), tfind(), tsearch(), twalk() */
#include <stdio.h>  /* For fprintf() */
#include <stdlib.h> /* For exit() */
#include <string.h> /* For strlen() */
#include <strings.h> /* For strcasecmp() */
#include <ctype.h>  /* For isalnum() */
#include <unistd.h> /* For gethostname() */

static void *_twalk_data; /* data to pass to twalk() */

static ArgumentOption *ArgumentOption_create(ArgumentSet *as,
                                             gchar symbol,
                                             gchar *option,
                                             gchar *type,
                                             gchar *desc,
                                             gchar *default_string,
                                             ArgumentHandler handler,
                                             gpointer handler_data){
    register ArgumentOption *ao= g_new(ArgumentOption, 1);
    ao->as = as;
    ao->symbol = symbol;
    ao->option = g_strdup(option);
    ao->arg_value = NULL;
    ao->desc = g_strdup(desc);
    ao->type = g_strdup(type);
    ao->default_string = g_strdup(default_string);
    ao->handler = handler;
    ao->handler_data = handler_data;
    ao->env_var = NULL;
    if(handler)
        ao->arg_list = NULL;
    else
        ao->arg_list = g_ptr_array_new();
    return ao;
    }

static void ArgumentOption_destroy(ArgumentOption *ao){
    g_free(ao->option);
    g_free(ao->desc);
    g_free(ao->type);
    g_free(ao->default_string);
    g_free(ao->env_var);
    if(ao->arg_list)
        g_ptr_array_free(ao->arg_list, TRUE);
    g_free(ao);
    return;
    }

static int Argument_strcmp_compare(const void *node1, const void *node2){
    return strcmp(((const ArgumentOption *)node1)->option, 
                  ((const ArgumentOption *)node2)->option);
    }

static void ArgumentSet_nodes_destroy(ArgumentSet *as, void **option_registry){
    int i;
    ArgumentOption *ao;
    for(i = 0; i < as->arg_option->len; i++){
        ao = as->arg_option->pdata[i];
        tdelete((void *)ao, option_registry, Argument_strcmp_compare);
        }
    }
static void ArgumentSet_destroy(ArgumentSet *as){
    register gint i;
    register ArgumentOption *ao;
    for(i = 0; i < as->arg_option->len; i++){
        ao = as->arg_option->pdata[i];
        ArgumentOption_destroy(ao);
        }
    g_ptr_array_free(as->arg_option, TRUE);
    g_free(as->desc);
    g_free(as);
    return;
    }

static void ArgumentOption_print_values(ArgumentOption *ao,
                                        gint use_long){
    register gint i;
    if(ao->arg_list){
        if(ao->arg_list->len){
            if((!ao->default_string)
            || ((ao->arg_list->len == 1)
                && strcmp(ao->default_string,
                          ao->arg_list->pdata[0]))){
                g_print("{%s", (gchar*)ao->arg_list->pdata[0]);
                for(i = 1; i < ao->arg_list->len; i++)
                    g_print(":%s", (gchar*)ao->arg_list->pdata[i]);
                g_print("}");
                }
        } else {
            g_print(" <*** empty list ***>");
            }
    } else {
        if(ao->arg_value){
            if((!ao->default_string)
            || strcmp(ao->default_string, ao->arg_value)){
                if(use_long)
                    g_print("Using: ");
                g_print("<%s>", ao->arg_value);
                }
        } else {
            g_print(" <*** not set ***>");
            }
        }
    return;
    }

static gboolean ArgumentOption_is_set(ArgumentOption *ao){
    if(ao->arg_list)
        return ao->arg_list->len?TRUE:FALSE;
    return ao->arg_value?TRUE:FALSE;
    }

static void ArgumentOption_add_value(ArgumentOption *ao,
                               GPtrArray *error_queue, gchar *value){
    if(ao->arg_list){
        g_ptr_array_add(ao->arg_list, value);
    } else {
        if(ao->arg_value)
            g_ptr_array_add(error_queue,
                g_strdup_printf("Already set --%s to \"%s\"",
                    ao->option, ao->arg_value));
        else
            ao->arg_value = value;
        }
    return;
    }

static gboolean ArgumentOption_is_mandatory(ArgumentOption *ao){
    return ao->default_string?FALSE:TRUE;
    }

/**/

static gboolean ArgumentParse_boolean(gchar *arg_string){
    register gint i;
    gchar *true_string[6]  = {"Y", "T", "TRUE",  "YES", "ON", "1"},
          *false_string[6] = {"N", "F", "FALSE", "NO",  "OFF", "0"};
    for(i = 0; i < 6; i++){
        if(!strcasecmp(arg_string, true_string[i]))
            return TRUE;
        if(!strcasecmp(arg_string, false_string[i]))
            return FALSE;
        }
    g_error("Cannot parse boolean \"%s\"", arg_string);
    return FALSE;
    }
/* FIXME: change to work with error listing
 */

static gchar *ArgumentHandler_short_help_func(gchar *arg_string,
                                              gpointer data){
    register Argument *arg = (Argument*)data;
    if(ArgumentParse_boolean(arg_string))
        arg->show_short_help = TRUE;
    return NULL;
    }

static gchar *ArgumentHandler_long_help_func(gchar *arg_string,
                                             gpointer data){
    register Argument *arg = (Argument*)data;
    if(ArgumentParse_boolean(arg_string))
        arg->show_long_help = TRUE;
    return NULL;
    }

static void Argument_show_version(Argument *arg){
    register gchar *branch = "$Name:  $";
    g_print("%s from %s version %s\n",
            arg->name, PACKAGE, VERSION);
    g_print("Using glib version %d.%d.%d\n",
            GLIB_MAJOR_VERSION,
            GLIB_MINOR_VERSION,
            GLIB_MICRO_VERSION);
    g_print("Built on %s\n", __DATE__);
    if(strlen(branch) >= 10)
        g_print("Branch: %.*s\n", (gint)(strlen(branch)-9), branch+7);
    return;
    }

static gchar *ArgumentHandler_version_func(gchar *arg_string,
                                           gpointer data){
    register Argument *arg = (Argument*)data;
    if(ArgumentParse_boolean(arg_string)){
        Argument_show_version(arg);
        exit(1);
        }
    return NULL;
    }

static void Argument_add_standard_options(Argument *arg){
    register ArgumentSet *as = ArgumentSet_create("General Options");
    ArgumentSet_add_option(as, 'h', "shorthelp", NULL,
        "Display compact help text",
        "FALSE", ArgumentHandler_short_help_func, arg);
    ArgumentSet_add_option(as, '\0', "help", NULL,
        "Displays verbose help text",
        "FALSE", ArgumentHandler_long_help_func, arg);
    ArgumentSet_add_option(as, 'v', "version", NULL,
                           "Show version number for this program",
                           "FALSE", ArgumentHandler_version_func, arg);
    Argument_absorb_ArgumentSet(arg, as);
    return;
    }

/**/

typedef struct {
    Argument_Cleanup_Func cleanup_func;
                 gpointer user_data;
} Argument_Cleanup;

void Argument_add_cleanup(Argument *arg,
                          Argument_Cleanup_Func cleanup_func,
                          gpointer user_data){
    register Argument_Cleanup *cleanup = g_new(Argument_Cleanup, 1);
    cleanup->cleanup_func = cleanup_func;
    cleanup->user_data = user_data;
    g_ptr_array_add(arg->cleanup_list, cleanup);
    return;
    }

static void Argument_cleanup(Argument *arg){
    register Argument_Cleanup *cleanup;
    register gint i;
    for(i = 0; i < arg->cleanup_list->len; i++){
        cleanup = arg->cleanup_list->pdata[i];
        cleanup->cleanup_func(cleanup->user_data);
        g_free(cleanup);
        }
    return;
    }

/**/

static Argument *Argument_create(gint argc, gchar **argv){
    register Argument *arg = g_new0(Argument, 1);
    arg->arg_set = g_ptr_array_new();
    arg->mandatory_set = g_ptr_array_new();
    arg->cleanup_list = g_ptr_array_new();
    arg->argc = argc;
    arg->argv = argv;
    Argument_add_standard_options(arg);
    return arg;
    }

static void Argument_destroy(Argument *arg){
    register ArgumentSet *as;
    register gint i;

    for(i = 0; i < arg->arg_set->len; i++) {
        as = arg->arg_set->pdata[i];
        ArgumentSet_nodes_destroy(as, &arg->option_registry);
        }

    for(i = 0; i < arg->arg_set->len; i++){
        as = arg->arg_set->pdata[i];
        ArgumentSet_destroy(as);
        }
    Argument_cleanup(arg);
    g_ptr_array_free(arg->cleanup_list, TRUE);
    g_ptr_array_free(arg->mandatory_set, TRUE);
    g_ptr_array_free(arg->arg_set, TRUE);
    g_free(arg->name);
    g_free(arg->desc);
    g_free(arg);
    return;
    }

static gboolean Argument_assertion_warning(void){
    g_warning("Compiled with assertion checking - will run slowly");
    return TRUE;
    }

static void Argument_error_handler(const gchar *log_domain,
                                   GLogLevelFlags log_level,
                                   const gchar *message,
                                   gpointer user_data){
    register Argument *arg = user_data;
    register gchar
        *stack_trace_str = (gchar*)g_getenv("EXONERATE_DEBUG_STACK_TRACE"),
        *debug_str = (gchar*)g_getenv("EXONERATE_DEBUG");
    fprintf(stderr, "** FATAL ERROR **: %s\n", message);
    if(stack_trace_str
    && ArgumentParse_boolean(stack_trace_str)){
        fprintf(stderr, "Generating stack trace ...\n");
        g_on_error_stack_trace(arg->name);
        }
    if(debug_str
    && ArgumentParse_boolean(debug_str)){
        fprintf(stderr, "Calling abort...\n");
        abort();
        }
    fprintf(stderr, "exiting ...\n");
    Argument_destroy(arg);
    exit(1);
    return;
    }
/* This is probably not what one is supposed to do with glib,
 * but it is to stop g_error() calling abort() and core dumping.
 * Maybe switch to g_critical() after glib-2 migration.
 */

int main(int argc, char **argv){
    register Argument *arg;
    register gint retval;
#ifdef USE_PTHREADS
    if(!g_thread_supported())
        g_thread_init(NULL);
#endif /* USE_PTHREADS */
    arg = Argument_create(argc, argv);
    g_log_set_handler(NULL, G_LOG_LEVEL_ERROR|G_LOG_FLAG_FATAL,
                      Argument_error_handler, arg);
    g_assert(Argument_assertion_warning());
    retval = Argument_main(arg);
    Argument_destroy(arg);
    return retval;
    }

static void Argument_usage(Argument *arg, gchar *synopsis){
    register ArgumentOption *ao;
    register gint i;
    Argument_show_version(arg);
    g_print("\n%s: %s\n", arg->name, arg->desc);
    if(synopsis){
        g_print("%s\n", synopsis);
    } else {
        g_print("Synopsis:\n--------\n%s", arg->name);
        for(i = 0; i < arg->mandatory_set->len; i++){
            ao = arg->mandatory_set->pdata[i];
            g_print(" <%s>", ao->type);
            }
        g_print("\n\n");
        }
    return;
    }

static void Argument_short_help(Argument *arg){
    register ArgumentSet *as;
    register ArgumentOption *ao;
    register gint i, j;
    for(i = 0; i < arg->arg_set->len; i++){
        as = arg->arg_set->pdata[i];
        if(as->arg_option->len){
            g_print("%s:\n", as->desc);
            for(j = strlen(as->desc); j > 0; j--)
                g_print("-");
            g_print("\n");
            for(j = 0; j < as->arg_option->len; j++){
                ao = as->arg_option->pdata[j];
                if(ao->symbol)
                    g_print("-%c ", ao->symbol);
                else
                    g_print("   ");
                g_print("--%s", ao->option);
                if(ao->default_string)
                    g_print(" [%s]", ao->default_string);
                else
                    g_print(" [mandatory]");
                g_print(" ");
                ArgumentOption_print_values(ao, FALSE);
                g_print("\n");
                }
            g_print("\n");
            }
        }
    g_print("--\n");
    return;
    }

static void Argument_long_help(Argument *arg){
    register ArgumentSet *as;
    register ArgumentOption *ao;
    register gint i, j;
    register gchar *env_value;
    for(i = 0; i < arg->arg_set->len; i++){
        as = arg->arg_set->pdata[i];
        if(as->arg_option->len){
            g_print("%s:\n", as->desc);
            for(j = strlen(as->desc); j > 0; j--)
                g_print("-");
            g_print("\n\n");
            for(j = 0; j < as->arg_option->len; j++){
                ao = as->arg_option->pdata[j];
                if(ao->symbol)
                    g_print("-%c ", ao->symbol);
                g_print("--%s", ao->option);
                if(ao->type)
                    g_print(" <%s>", ao->type);
                g_print("\n%s\n", ao->desc);
                g_print("Environment variable: $%s", ao->env_var);
                env_value = (gchar*)g_getenv(ao->env_var);
                if(env_value){
                    g_print(" (Set to \"%s\")\n", env_value);
                } else {
                    g_print(" (Not set)\n");
                    }
                if(ao->default_string)
                    g_print("Default: \"%s\"\n", ao->default_string);
                else
                    g_print("*** This argument is mandatory ***\n");
                ArgumentOption_print_values(ao, TRUE);
                g_print("\n");
                }
            }
        }
    g_print("--\n");
    return;
    }

/* set _twalk_data before calling twalk(tree, Argument_set_env_var_func) */
static void Argument_set_env_var_func(const void *ptr,
                                      VISIT order,
                                      int level){
    ArgumentOption *ao = *(ArgumentOption **)ptr;
    int i;
    ao->env_var = g_strdup_printf("%s_%s_%s",
                  PACKAGE, (char *)_twalk_data, ao->option);
    for(i = strlen(ao->env_var)-1; i >= 0; i--)
        if(isalnum(ao->env_var[i]))
            ao->env_var[i] = toupper(ao->env_var[i]);
        else
            ao->env_var[i] = '_';
    }

/* set _twalk_data before calling twalk(tree,Argument_traverse_registry_func) */
static void Argument_traverse_registry_func(const void *ptr,
                                            VISIT order,
                                            int level){
    ArgumentOption *ao = *(ArgumentOption**)ptr;
    GPtrArray *error_queue = (GPtrArray*)_twalk_data;
    char *err_msg, *env_var;
    if(!ArgumentOption_is_set(ao)){
        env_var = (gchar*)g_getenv(ao->env_var);
        if(env_var)
            ArgumentOption_add_value(ao, error_queue, env_var);
        }
    if(!ArgumentOption_is_set(ao)){
        if(ArgumentOption_is_mandatory(ao)){
            g_ptr_array_add(error_queue,
                g_strdup_printf(
                    "No value set for mandatory argument --%s <%s>",
                                            ao->option, ao->type));
        } else {
            ArgumentOption_add_value(ao, error_queue,
                                     ao->default_string);
            }
        }
    if(ao->handler){
        err_msg = ao->handler(ao->arg_value, ao->handler_data);
        if(err_msg)
            g_ptr_array_add(error_queue, err_msg);
    } else { /* List */
        (*((GPtrArray**)ao->handler_data)) = ao->arg_list;
        }
    }

static ArgumentOption *Argument_process_get_option(Argument *arg,
                           gchar *string, GPtrArray *error_queue){
    register ArgumentOption *ao = NULL;
    register gint i;
    if(string[0] == '-'){
        if(string[1] == '-'){
            void *tree_node = tfind((void*)&(ArgumentOption){.option=string+2},
                                    &arg->option_registry,
                                    Argument_strcmp_compare);
            ao = tree_node ? *(ArgumentOption **)tree_node : NULL;
            if(!ao)
                g_ptr_array_add(error_queue,
                     g_strdup_printf("Unrecognised option \"%s\"",
                                      string));
        } else {
            for(i = 1; string[i]; i++){
                ao = arg->symbol_registry[(guchar)string[i]];
                if(!ao)
                    g_error("Unknown flag [%c] in argument [%s]",
                            string[i], string);
                if(string[i+1]) /* If not last symbol */
                    ArgumentOption_add_value(ao, error_queue, "TRUE");
                }
            }
        }
    return ao;
    }

void Argument_info(Argument *arg){
    register gchar *cl = g_strjoinv(" ", arg->argv);
    gchar hostname[1024];
    g_print("Command line: [%s]\n", cl);
    g_free(cl);
    /**/
    gethostname(hostname, 1024);
    g_print("Hostname: [%s]\n", hostname);
    return;
    }

void Argument_process(Argument *arg, gchar *name, gchar *desc,
                      gchar *synopsis){
    register gint i;
    register GPtrArray *unflagged_arg = g_ptr_array_new();
    register ArgumentOption *ao;
    register GPtrArray *error_queue = g_ptr_array_new();
    arg->desc = desc?g_strdup(desc):g_strdup("");
    arg->name = g_strdup(name);
    _twalk_data = (void *)arg->name; /* used in Argument_set_env_var_func() */ 
    twalk(arg->option_registry, Argument_set_env_var_func);
    if(arg->mandatory_set->len && (arg->argc <= 1)){
        Argument_usage(arg, synopsis);
        exit(1);
        }
    for(i = 1; i < arg->argc; i++){
        ao = Argument_process_get_option(arg, arg->argv[i],
                                         error_queue);
        if(ao){ /* -? */
            if(ao->type){ /* Not boolean */
                if((i+1) == arg->argc){
                    g_ptr_array_add(error_queue,
                        g_strdup_printf("No argument supplied with %s",
                            arg->argv[i]));
                } else {
                    if(ao->handler){ /* Not list */
                        ArgumentOption_add_value(ao, error_queue,
                                                 arg->argv[++i]);
                    } else { /* Is list */
                        do {
                            ArgumentOption_add_value(ao, error_queue,
                                                     arg->argv[i+1]);
                            i++;
                        } while(((i+1) < arg->argc)
                             && (arg->argv[i+1][0] != '-'));
                        }
                    }
            } else { /* Is boolean */
                if((i+1) == arg->argc){
                    ArgumentOption_add_value(ao, error_queue, "TRUE");
                } else {
                    if(arg->argv[i+1][0] == '-'){
                        ArgumentOption_add_value(ao, error_queue,
                                                 "TRUE");
                    } else {
                        ArgumentOption_add_value(ao, error_queue,
                                                 arg->argv[i+1]);
                        i++;
                        }
                    }
                }
        } else { /* unflagged argument */
            g_ptr_array_add(unflagged_arg, arg->argv[i]);
            }
        }
    if(arg->mandatory_set->len < unflagged_arg->len){
        g_ptr_array_add(error_queue,
          g_strdup_printf("Too many unflagged arguments"));
    } else {
        for(i = 0; i < unflagged_arg->len; i++){
            ao = arg->mandatory_set->pdata[i];
            ArgumentOption_add_value(ao, error_queue,
                                     unflagged_arg->pdata[i]);
            }
        }
    g_ptr_array_free(unflagged_arg, TRUE);
    _twalk_data = (void *)error_queue;
    twalk(arg->option_registry, Argument_traverse_registry_func);
    if((!(arg->show_short_help | arg->show_long_help)
        && error_queue->len)){
        Argument_usage(arg, synopsis);
        g_print("--\n"
                "%d ERROR%s encountered in argument processing\n"
                "--\n",
                error_queue->len,
                (error_queue->len > 1)?"S were":" was");
        for(i = 0; i < error_queue->len; i++){
            g_print("[ %d ] : %s\n",
                    i+1, (gchar*)error_queue->pdata[i]);
            g_free(error_queue->pdata[i]);
            }
        g_print("--\n\n"
                "Use -h or --help for more information on usage\n"
                "\n");
        exit(1);
        }
    if(arg->show_short_help){
        Argument_usage(arg, synopsis);
        Argument_short_help(arg);
        exit(1);
        }
    if(arg->show_long_help){
        Argument_usage(arg, synopsis);
        Argument_long_help(arg);
        exit(1);
        }
    g_ptr_array_free(error_queue, TRUE);
    return;
    }
/* FIXME: tidy */

ArgumentSet *ArgumentSet_create(gchar *desc){
    register ArgumentSet *as = g_new(ArgumentSet, 1);
    g_assert(desc);
    as->desc = g_strdup(desc);
    as->arg_option = g_ptr_array_new();
    return as;
    }

void Argument_absorb_ArgumentSet(Argument *arg, ArgumentSet *as){
    register ArgumentOption *ao;
    register gint i;
    g_ptr_array_add(arg->arg_set, as);
    for(i = 0; i < as->arg_option->len; i++){
        ao = as->arg_option->pdata[i];
        if(!ao->default_string)
            g_ptr_array_add(arg->mandatory_set, ao);
        if(ao->symbol){ /* Check option not already used */
            g_assert(!arg->symbol_registry[(guchar)ao->symbol]);
            arg->symbol_registry[(guchar)ao->symbol] = ao;
            }
        g_assert(!tfind((void *)ao, &arg->option_registry, 
                      Argument_strcmp_compare));
        tsearch((void*)ao, &arg->option_registry, Argument_strcmp_compare);
        }
    return;
    }

void ArgumentSet_add_option(ArgumentSet *as,
                            gchar symbol,
                            gchar *option,
                            gchar *type,
                            gchar *desc,
                            gchar *default_string,
                            ArgumentHandler handler,
                            gpointer handler_data){
    register ArgumentOption *ao = ArgumentOption_create(
             as, symbol, option, type, desc,
             default_string, handler, handler_data);
    g_ptr_array_add(as->arg_option, ao);
    return;
    }

gchar *Argument_parse_string(gchar *arg_string, gpointer data){
    register gchar **dst_string = (gchar**)data;
    if(!strcasecmp(arg_string, "NULL"))
        (*dst_string) = NULL;
    else
        (*dst_string) = arg_string;
    return NULL;
    }

gchar *Argument_parse_char(gchar *arg_string, gpointer data){
    register gchar *dst_char = (gchar*)data;
    if((!arg_string[0]) || (arg_string[1]))
        return g_strdup_printf(
                "Expected single character argument not [%s]",
                arg_string);
    (*dst_char) = arg_string[0];
    return NULL;
    }

gchar *Argument_parse_int(gchar *arg_string, gpointer data){
    register gint *dst_int = (gint*)data;
    (*dst_int) = atoi(arg_string);
    return NULL;
    }

gchar *Argument_parse_float(gchar *arg_string, gpointer data){
    register gfloat *dst_float = (gfloat*)data;
    (*dst_float) = atof(arg_string);
    return NULL;
    }

gchar *Argument_parse_boolean(gchar *arg_string, gpointer data){
    register gboolean *dst_boolean = (gboolean*)data;
    (*dst_boolean) = ArgumentParse_boolean(arg_string);
    return NULL;
    }

/**/

