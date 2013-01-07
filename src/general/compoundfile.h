/****************************************************************\
*                                                                *
*  Library for reading large and/or split files.                 *
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

#ifndef INCLUDED_COMPOUNDFILE_H
#define INCLUDED_COMPOUNDFILE_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <glib.h>
#include <stdio.h>

#include <sys/types.h> /* For lseek() and off_t */
#include <unistd.h>    /* For lseek() and off_t */
#include <inttypes.h>  /* For PRIuPTR           */

#ifdef USE_PTHREADS
#include <pthread.h>
#endif /* USE_PTHREADS */

#define COMPOUND_FILE_BUFFER_SIZE BUFSIZ
/* Using BUFSIZ from stdio.h, although I'm not sure
 * why the stdio buffered IO with getc isn't just as
 * fast as this implementation.
 */

/* G_GNUC_EXTENSION is not defined on all platforms */
#ifdef G_GNUC_EXTENSION
#define CompoundFile_Pragma G_GNUC_EXTENSION
#else /* G_GNUC_EXTENSION */
#define CompoundFile_Pragma
#endif /* G_GNUC_EXTENSION */

typedef off_t CompoundFile_Pos;

#if __STDC_VERSION__ >= 199901L
#define CompoundFile_Pos_print(fp, pos) \
    CompoundFile_Pragma fprintf((fp), "%jd", (intmax_t)(pos))
typedef intmax_t CompoundFile_Pos_scan_type;
#define CompoundFile_Pos_scan(src, dst) \
    CompoundFile_Pragma sscanf((src), "%jd", (intmax_t*)(dst))
#else /* __STDC_VERSION__  */
#define CompoundFile_Pos_print(fp, pos) \
    CompoundFile_Pragma fprintf((fp), "%lld", (long long int)(pos))
typedef long long int CompoundFile_Pos_scan_type;
#define CompoundFile_Pos_scan(src, dst) \
    CompoundFile_Pragma sscanf((src), "%lld", (dst))
#endif /* __STDC_VERSION__  */

typedef struct {
              gchar  *path;
    CompoundFile_Pos  length;
} CompoundFile_Element;

typedef struct {
                           guint  ref_count;
                       GPtrArray *element_list;
                            FILE *fp;
                            gint  curr_element_id;
                           gchar  in_buffer[COMPOUND_FILE_BUFFER_SIZE];
                            gint  in_buffer_full;
                            gint  in_buffer_pos;
                CompoundFile_Pos  in_buffer_start;
    struct CompoundFile_Location *start_limit;
    struct CompoundFile_Location *stop_limit;
#ifdef USE_PTHREADS
                 pthread_mutex_t  compoundfile_mutex;
#endif /* USE_PTHREADS */
} CompoundFile;

CompoundFile *CompoundFile_create(GPtrArray *path_list,
                                  gboolean sort_on_file_size);
        void  CompoundFile_destroy(CompoundFile *cf);
CompoundFile *CompoundFile_share(CompoundFile *cf);
CompoundFile *CompoundFile_dup(CompoundFile *cf);
       gchar *CompoundFile_current_path(CompoundFile *cf);

#define CompoundFile_buffer_is_empty(cf)              \
        ((cf)->in_buffer_pos == (cf)->in_buffer_full)

#define CompoundFile_getc(cf)                      \
        (CompoundFile_buffer_is_empty(cf)          \
         ? CompoundFile_buffer_reload(cf)          \
         : (cf)->in_buffer[(cf)->in_buffer_pos++])

gint CompoundFile_buffer_reload(CompoundFile *cf);
/* This should only be called by the CompoundFile_getc() macro */

#define CompoundFile_ftell(cf) \
    ((cf)->in_buffer_start + (cf)->in_buffer_pos)
/* Alternatively, could use lseek(fileno(cf->fp), 0, SEEK_CUR) */

gboolean CompoundFile_is_finished(CompoundFile *cf);
    void CompoundFile_rewind(CompoundFile *cf);
    void CompoundFile_seek(CompoundFile *cf, CompoundFile_Pos pos);
CompoundFile_Pos CompoundFile_get_length(CompoundFile *cf);
CompoundFile_Pos CompoundFile_get_max_element_length(CompoundFile *cf);

 gint  CompoundFile_read(CompoundFile *cf, gchar *buf, gint length);

/**/

typedef struct CompoundFile_Location {
                gint  ref_count;
        CompoundFile *cf;
    CompoundFile_Pos  pos;
                gint  element_id;
} CompoundFile_Location;

CompoundFile_Location *CompoundFile_Location_create(CompoundFile *cf,
                             CompoundFile_Pos pos, gint element_id);
CompoundFile_Location *CompoundFile_Location_tell(CompoundFile *cf);
                 void  CompoundFile_Location_destroy(
                       CompoundFile_Location *cfl);
CompoundFile_Location *CompoundFile_Location_share(
                       CompoundFile_Location *cfl);
                 void  CompoundFile_Location_seek(
                       CompoundFile_Location *cfl);

/**/

CompoundFile_Location *CompoundFile_Location_from_pos(
                       CompoundFile *cf, CompoundFile_Pos pos);
/* Converts a position (in compound file coordinates) to a location.
 */

void CompoundFile_set_limits(CompoundFile *cf,
                             CompoundFile_Location *start,
                             CompoundFile_Location *stop);
/* Sets limits and rewinds truncated file */

/**/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* INCLUDED_COMPOUNDFILE_H */

