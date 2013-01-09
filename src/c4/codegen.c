/****************************************************************\
*                                                                *
*  Code generation module for C4                                 *
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

#include <stdlib.h> /* For system() */
#include <string.h> /* For strlen() */
#include <ctype.h>  /* For isalnum() */

#include <sys/stat.h>  /* For stat(), mkdir() */
#include <sys/types.h> /* For stat(), mkdir() */

#include "codegen.h"

Codegen_ArgumentSet *Codegen_ArgumentSet_create(Argument *arg){
    register ArgumentSet *as;
    static Codegen_ArgumentSet cas = {TRUE};
    if(arg){
        as = ArgumentSet_create("Code generation options");
        ArgumentSet_add_option(as, 'C', "compiled", NULL,
                "Use compiled viterbi implementations", "TRUE",
                Argument_parse_boolean, &cas.use_compiled);
        Argument_absorb_ArgumentSet(arg, as);
        }
    return &cas;
    }

/**/

gchar *Codegen_clean_path_component(gchar *name){
    register gint i;
    register GString *file_name_str = g_string_sized_new(strlen(name));
    register gchar *file_name, *tmp;
    g_assert(name);
    for(i = 0; name[i]; i++){
        if(isalnum(name[i]) || (name[i] == '_')){
            g_string_append_c(file_name_str, name[i]);
        } else {
            tmp = g_strdup_printf("_%d_", name[i]);
            g_string_append(file_name_str, tmp);
            g_free(tmp);
            }
        }
    file_name = file_name_str->str;
    g_string_free(file_name_str, FALSE);
    return file_name;
    }

gboolean Codegen_directory_exists(gchar *path){
    struct stat result;
    if(!stat(path, &result))
        if(S_ISDIR(result.st_mode))
            return TRUE;
    return FALSE;
    }

gboolean Codegen_file_exists(gchar *path){
    struct stat result;
    if(!stat(path, &result))
        if(S_ISREG(result.st_mode))
            return TRUE;
    return FALSE;
    }

static void Codegen_directory_create(gchar *path){
    register mode_t mode = S_IRUSR|S_IWUSR|S_IXUSR
                         | S_IRGRP|S_IXGRP
                         | S_IROTH|S_IXOTH;
    g_warning("Creating codegen directory [%s]", path);
    if(mkdir(path, mode))
        g_error("Could not create directory [%s]", path);
    else
        g_warning("Creating codegen directory [%s]", path);
    return;
    }

static gchar *Codegen_get_code_dir(gchar *directory){
    register gchar *platform_directory, *clean_hosttype
        = Codegen_clean_path_component(HOSTTYPE);
    char code_dir[_XOPEN_PATH_MAX];
    if(!directory){
        if (getenv("C4_CODEGEN_DIRECTORY"))
            strncpy(code_dir, getenv("C4_CODEGEN_DIRECTORY"),
                    sizeof(code_dir)-1);
        else
            snprintf(code_dir, sizeof(code_dir), "%s/codegen",
                     SOURCE_ROOT_DIR);
        }
    else
        strncpy(code_dir, directory, sizeof(code_dir)-1);
    if(!Codegen_directory_exists(code_dir))
        Codegen_directory_create(code_dir);
    platform_directory = g_strconcat(code_dir, G_DIR_SEPARATOR_S,
                                     clean_hosttype, NULL);
    if(!Codegen_directory_exists(platform_directory))
        Codegen_directory_create(platform_directory);
    g_free(clean_hosttype);
    return platform_directory;
    }

Codegen *Codegen_create(gchar *directory, gchar *name){
    register Codegen *c = g_new(Codegen, 1);
    register gchar *code_dir = Codegen_get_code_dir(directory);
    c->name = Codegen_clean_path_component(name);
    c->code_path = g_strconcat(code_dir, G_DIR_SEPARATOR_S,
                               c->name, ".c", NULL);
    c->object_path = g_strconcat(code_dir, G_DIR_SEPARATOR_S,
                               c->name, ".o", NULL);
    c->fp = fopen(c->code_path, "w");
    if(!c->fp){
        perror("Writing codegen file");
        g_error("Could not write codegen code to [%s]", c->code_path);
        }
    c->indent = 0;
    g_free(code_dir);
    return c;
    }

void Codegen_destroy(Codegen *c){
    g_assert(c);
    if(c->fp)
        fclose(c->fp);
    /* remove(c->code_path); */
    g_free(c->code_path);
    g_free(c->name);
    g_free(c);
    return;
    }

static gchar *Codegen_expand_format(Codegen *c, gchar *format){
    register GString *str = g_string_sized_new(strlen(format));
    register gint i, j;
    register gchar *expanded_format;
    for(i = 0; format[i]; i++){
        g_string_append_c(str, format[i]);
        if((format[i] == '\n')     /* If there is a newline */
        && format[i+1]             /* not at the end */
        && (format[i+1] != '\n')){ /* or preceeding another newline */
            for(j = 0; j < ((c->indent)<<2); j++)
                g_string_append_c(str, ' ');
            }
        }
    expanded_format = str->str;
    g_string_free(str, FALSE);
    return expanded_format;
    }
/* Allows tidy printing of multi-line statements */

void Codegen_printf(Codegen *c, gchar *format, ...){
    va_list args;
    register gchar *code;
    format = Codegen_expand_format(c, format);
    /* Reuse expanded format to stop compiler complaining */
    va_start(args, format);
    code = g_strdup_vprintf(format, args);
    va_end(args);
    fprintf(c->fp, "%*s%s", (c->indent<<2), "", code);
    g_free(format); /* Free expanded format */
    g_free(code);
    return;
    }

void Codegen_indent(Codegen *c, gint indent_change){
    g_assert((c->indent += indent_change) || TRUE);
    return;
    }
/* Code indentation is only active when assertion checking is on.
 */

/**/

void Codegen_compile(Codegen *c,
                    gchar *add_ccflags, gchar *add_ldflags){
    register gchar *cc_command = "gcc";
/* Optimisation is turned off when assertions
 * are on for faster compilation
 */
#ifdef G_DISABLE_ASSERT
    register gchar *cc_flags   = "-O2 -Wall";
#else /* G_DISABLE_ASSERT */
    register gchar *cc_flags   = "-g -Wall";
#endif /* G_DISABLE_ASSERT */
    /* FIXME: optimisation
     *        try some more aggressive compiler optimisations
     *        -fomit-frame-pointer -O3 etc
     *        or --arch ev6 / --arch ev67 with native cc on OSF1
     *        -O3 -tpp6 -xK
     */
    gchar *compile_command;
    register gchar *tmp;
    /* Allow customistation of compilation with environment variables */
    g_assert(c->indent == 0);
    fclose(c->fp);
    c->fp = NULL;
    /* If output already present, do not compile */
    if(Codegen_file_exists(c->object_path)){
        /* FIXME: unless object_path is older than model.o */
        g_warning("Reusing codegen object [%s]", c->object_path);
        return;
        }
    tmp = (gchar*)g_getenv("CC");
    if(tmp)
        cc_command = tmp;
    tmp = (gchar*)g_getenv("CFLAGS");
    if(tmp)
        cc_flags = tmp;
    /* FIXME: these enviroment variables should be merged
     *        with the command line options for this module.
     */
    /**/
    if(add_ccflags)
        cc_flags = g_strconcat(cc_flags, " ", add_ccflags, NULL);
    /**/
    compile_command = g_strconcat(cc_command, " ", cc_flags,
            " -o ", c->object_path, " -c ", c->code_path, NULL);
    g_message("Compiling codegen:\n%s", compile_command);
    if(system(compile_command))
        g_error("Problem compiling with\n%s", compile_command);
    g_free(compile_command);
    if(add_ccflags)
        g_free(cc_flags);
    return;
    }

