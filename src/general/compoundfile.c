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

#include "compoundfile.h"

#include <stdlib.h> /* For qsort() */

#include <errno.h>
#include <string.h> /* For strerror() */

/**/

static CompoundFile_Element *CompoundFile_Element_create(
                             gchar *path,
                             CompoundFile_Pos length){
    register CompoundFile_Element *cfe
     = g_new(CompoundFile_Element, 1);
    cfe->path = g_strdup(path);
    cfe->length = length;
    return cfe;
    }

static void CompoundFile_Element_destroy(CompoundFile_Element *cfe){
    g_free(cfe->path);
    g_free(cfe);
    return;
    }

/**/

static FILE *CompoundFile_open_file(gchar *path){
    register FILE *fp = fopen(path, "r");
    if(!fp)
        g_error("Could not open [%s] : %s", path, strerror(errno));
    return fp;
    }

static void CompoundFile_clear_buffers(CompoundFile *cf,
                                       CompoundFile_Pos pos){
    cf->in_buffer_full = 0;
    cf->in_buffer_pos = 0;
    cf->in_buffer_start = pos;
    return;
    }

static void CompoundFile_set_file(CompoundFile *cf, gint element_id){
    register CompoundFile_Element *cfe;
    if(cf->curr_element_id == element_id)
        return;
    if(cf->fp)
        fclose(cf->fp);
    cfe = cf->element_list->pdata[element_id];
    cf->fp = CompoundFile_open_file(cfe->path);
    cf->curr_element_id = element_id;
    CompoundFile_clear_buffers(cf, 0);
    return;
    }

static CompoundFile_Pos CompoundFile_size_from_path(gchar *path){
    register FILE *fp;
    register CompoundFile_Pos size;
    fp = CompoundFile_open_file(path);
    size = lseek(fileno(fp), 0, SEEK_END);
    if(size == -1)
        g_error("Could not seek in file [%s] (%s)", path, strerror(errno));
    fclose(fp);
    return size;
    }

static int CompoundFile_sort_path_list(const void *a, const void *b){
    register CompoundFile_Element
        **element_a = (CompoundFile_Element**)a,
        **element_b = (CompoundFile_Element**)b;
    register gint ret_val = (*element_a)->length - (*element_b)->length;
    if(!ret_val)
        ret_val = strcmp((*element_b)->path, (*element_a)->path);
    if(!ret_val)
        g_error("Duplicate file [%s] in input", (*element_a)->path);
    return ret_val;
    }
/* Sorted incase shell order is not preserved,
 * to check against duplicate files,
 * and to search smallest files first for best performance
 * with dynamic thresholds.
 */

static void CompoundFile_sort(CompoundFile *cf){
    qsort(cf->element_list->pdata, cf->element_list->len,
          sizeof(gpointer), CompoundFile_sort_path_list);
    return;
    }

CompoundFile *CompoundFile_create(GPtrArray *path_list,
                                  gboolean sort_on_file_size){
    register CompoundFile *cf = g_new(CompoundFile, 1);
    register CompoundFile_Element *cfe;
    register CompoundFile_Pos length;
    register gint i;
    register gchar *path;
    g_assert(path_list);
    g_assert(path_list->len);
    cf->ref_count = 1;
    cf->element_list = g_ptr_array_new();
    cf->fp = NULL;
    cf->curr_element_id = -1;
    /* Expand any path list to include any directories */
    for(i = 0; i < path_list->len; i++){
        path = path_list->pdata[i];
        /* Finding the size also checks that the file exists */
        length = CompoundFile_size_from_path(path);
        cfe = CompoundFile_Element_create(path, length);
        g_ptr_array_add(cf->element_list, cfe);
        }
    cf->start_limit = cf->stop_limit = NULL;
    if(sort_on_file_size)
        CompoundFile_sort(cf);
    CompoundFile_rewind(cf);
#ifdef USE_PTHREADS
    pthread_mutex_init(&cf->compoundfile_mutex, NULL);
#endif /* USE_PTHREADS */
    return cf;
    }

void CompoundFile_destroy(CompoundFile *cf){
    register gint i;
    register CompoundFile_Element *cfe;
    if(--cf->ref_count)
        return;
    for(i = 0; i < cf->element_list->len; i++){
        cfe = cf->element_list->pdata[i];
        CompoundFile_Element_destroy(cfe);
        }
    g_ptr_array_free(cf->element_list, TRUE);
    fclose(cf->fp);
    if(cf->start_limit)
        CompoundFile_Location_destroy(cf->start_limit);
    if(cf->stop_limit)
        CompoundFile_Location_destroy(cf->stop_limit);
#ifdef USE_PTHREADS
    pthread_mutex_destroy(&cf->compoundfile_mutex);
#endif /* USE_PTHREADS */
    g_free(cf);
    return;
    }

CompoundFile *CompoundFile_share(CompoundFile *cf){
    cf->ref_count++;
    return cf;
    }

static GPtrArray *CompoundFile_get_path_list(CompoundFile *cf){
    register gint i;
    register gchar *path;
    register GPtrArray *path_list = g_ptr_array_new();
    register CompoundFile_Element *element;
    for(i = 0; i < cf->element_list->len; i++){
        element = cf->element_list->pdata[i];
        path = g_strdup(element->path);
        g_ptr_array_add(path_list, path);
        }
    return path_list;
    }

CompoundFile *CompoundFile_dup(CompoundFile *cf){
    register CompoundFile *ncf;
    register CompoundFile_Location *nstart = NULL, *nstop = NULL;
    register gint i;
    register gchar *path;
    register GPtrArray *path_list = CompoundFile_get_path_list(cf);
    g_assert(path_list);
    ncf = CompoundFile_create(path_list, FALSE);
    if(cf->start_limit)
        nstart = CompoundFile_Location_create(ncf,
                                              cf->start_limit->pos,
                                              cf->start_limit->element_id);
    if(cf->stop_limit)
        nstop  = CompoundFile_Location_create(ncf,
                                              cf->stop_limit->pos,
                                              cf->stop_limit->element_id);
    CompoundFile_set_limits(ncf, nstart, nstop);
    if(nstart)
        CompoundFile_Location_destroy(nstart);
    if(nstop)
        CompoundFile_Location_destroy(nstop);
    for(i = 0; i < path_list->len; i++){
        path = path_list->pdata[i];
        g_free(path);
        }
    g_ptr_array_free(path_list, TRUE);
    return ncf;
    }

gchar *CompoundFile_current_path(CompoundFile *cf){
    register CompoundFile_Element *cfe
        = cf->element_list->pdata[cf->curr_element_id];
    return cfe->path;
    }

gint CompoundFile_buffer_reload(CompoundFile *cf){
    register gsize read_size = sizeof(gchar)*COMPOUND_FILE_BUFFER_SIZE;
    g_assert(cf);
    cf->in_buffer_start += cf->in_buffer_full;
    if(cf->stop_limit)
        if(cf->curr_element_id == cf->stop_limit->element_id)
            if((cf->in_buffer_start + read_size) > cf->stop_limit->pos)
                read_size = cf->stop_limit->pos - cf->in_buffer_start;
    cf->in_buffer_full = read(fileno(cf->fp), cf->in_buffer, read_size);
    cf->in_buffer_pos = 0;
    if(!cf->in_buffer_full){
        if(cf->stop_limit){
            if(cf->curr_element_id == cf->stop_limit->element_id)
                return EOF;
        } else {
            if((cf->curr_element_id + 1) == cf->element_list->len)
                return EOF;
            }
        CompoundFile_set_file(cf, cf->curr_element_id+1);
        return CompoundFile_getc(cf);
        }
    return cf->in_buffer[cf->in_buffer_pos++];
    }

gboolean CompoundFile_is_finished(CompoundFile *cf){
    g_assert(cf);
    if(!CompoundFile_buffer_is_empty(cf))
        return FALSE;
    if(cf->in_buffer_full
       == (sizeof(gchar)*COMPOUND_FILE_BUFFER_SIZE))
        return FALSE; /* At end of file */
    if((cf->curr_element_id+1) < cf->element_list->len)
        return FALSE;
    return TRUE;
    }

void CompoundFile_rewind(CompoundFile *cf){
    if(cf->start_limit){
        CompoundFile_Location_seek(cf->start_limit);
    } else {
        CompoundFile_set_file(cf, 0);
        CompoundFile_seek(cf, 0);
        }
    return;
    }

void CompoundFile_seek(CompoundFile *cf, CompoundFile_Pos pos){
    register CompoundFile_Element *cfe;
    if(lseek(fileno(cf->fp), pos, SEEK_SET) != pos){
        cfe = cf->element_list->pdata[cf->curr_element_id];
        g_error("Could not seek in file [%s] (%s)",
                cfe->path, strerror(errno));
        }
    CompoundFile_clear_buffers(cf, pos);
    return;
    }
/* FIXME: optimisation: should not move within current buffer
 *        when possible
 */

CompoundFile_Pos CompoundFile_get_length(CompoundFile *cf){
    register gint i;
    register CompoundFile_Pos length = 0;
    register CompoundFile_Element *cfe;
    register gint start_element, stop_element;
    start_element = cf->start_limit?cf->start_limit->element_id
                                   :0;
    stop_element = cf->stop_limit?(cf->stop_limit->element_id + 1)
                                 :cf->element_list->len;
    for(i = start_element; i < stop_element; i++){
        cfe = cf->element_list->pdata[i];
        length += cfe->length;
        }
    if(cf->start_limit)
        length -= cf->start_limit->pos;
    if(cf->stop_limit){
        cfe = cf->element_list->pdata[cf->stop_limit->element_id];
        length -= (cfe->length - cf->stop_limit->pos);
        }
    return length;
    }

CompoundFile_Pos CompoundFile_get_max_element_length(CompoundFile *cf){
    register CompoundFile_Pos max = 0;
    register CompoundFile_Element *cfe;
    register gint i;
    for(i = 0; i < cf->element_list->len; i++){
        cfe = cf->element_list->pdata[i];
        if(max < cfe->length)
            max = cfe->length;
        }
    return max;
    }

/**/

gint CompoundFile_read(CompoundFile *cf, gchar *buf, gint length){
    g_error("ni");
    /* FIXME: do unbuffered read, line reload() */
    return -1;
    }

/**/

CompoundFile_Location *CompoundFile_Location_create(CompoundFile *cf,
                                                    CompoundFile_Pos pos,
                                                    gint element_id){
    register CompoundFile_Location *cfl
     = g_new(CompoundFile_Location, 1);
    g_assert(cf);
    cfl->ref_count = 1;
    cfl->cf = CompoundFile_share(cf);
    cfl->pos = pos,
    cfl->element_id = element_id;
    return cfl;
    }

CompoundFile_Location *CompoundFile_Location_tell(CompoundFile *cf){
    return CompoundFile_Location_create(cf, CompoundFile_ftell(cf),
                                        cf->curr_element_id);
    }

void CompoundFile_Location_destroy(CompoundFile_Location *cfl){
    g_assert(cfl);
    if(--cfl->ref_count)
        return;
    CompoundFile_destroy(cfl->cf);
    g_free(cfl);
    return;
    }

CompoundFile_Location *CompoundFile_Location_share(
                       CompoundFile_Location *cfl){
    g_assert(cfl);
    cfl->ref_count++;
    return cfl;
    }

void CompoundFile_Location_seek(CompoundFile_Location *cfl){
    g_assert(cfl);
    CompoundFile_set_file(cfl->cf, cfl->element_id);
    CompoundFile_seek(cfl->cf, cfl->pos);
    return;
    }

/**/

CompoundFile_Location *CompoundFile_Location_from_pos(
                       CompoundFile *cf, CompoundFile_Pos pos){
    register CompoundFile_Location *cfl
     = g_new(CompoundFile_Location, 1);
    register gint i;
    register CompoundFile_Element *cfe;
    cfl->ref_count = 1;
    cfl->cf = CompoundFile_share(cf);
    cfl->pos = pos;
    cfl->element_id = -1;
    /* Find which file this position should be in */
    for(i = 0; i < cf->element_list->len; i++){
        cfe = cf->element_list->pdata[i];
        if((cfl->pos - cfe->length) < 0){
            cfl->element_id = i;
            break;
        } else {
            cfl->pos -= cfe->length;
            }
        }
    g_assert(cfl->element_id != -1);
    return cfl;
    }

void CompoundFile_set_limits(CompoundFile *cf,
                             CompoundFile_Location *start,
                             CompoundFile_Location *stop){
    if(cf->start_limit)
        CompoundFile_Location_destroy(cf->start_limit);
    if(cf->stop_limit)
        CompoundFile_Location_destroy(cf->stop_limit);
    if(start)
        cf->start_limit = CompoundFile_Location_share(start);
    if(stop)
        cf->stop_limit = CompoundFile_Location_share(stop);
    CompoundFile_rewind(cf);
    return;
    }

/**/

