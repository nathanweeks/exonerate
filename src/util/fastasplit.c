/****************************************************************\
*                                                                *
*  fastasplit: split a fasta format file into chunks             *
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

#include <libgen.h> /* For basename() */
#include <sys/types.h> /* For fstat() */
#include <sys/stat.h>  /* For fstat() */
#include <unistd.h>    /* For fstat() */

#include "argument.h"
#include "fastadb.h"

static void fasta_chunk_write(FastaDB *fdb, CompoundFile_Pos start,
                              CompoundFile_Pos stop, gchar *path){
    register FILE *fp;
    register CompoundFile_Pos pos;
    if(start == stop) /* Don't write an empty file */
        return;
    fp = fopen(path, "w");
    if(!fp)
        g_error("Could not open [%s] for writing fasta chunk", path);
    CompoundFile_seek(fdb->cf, start);
    for(pos = start; pos < stop; pos++)
        putc(CompoundFile_getc(fdb->cf), fp);
    fclose(fp);
    return;
    }

static CompoundFile_Pos fasta_split_get_file_size(FILE *fp){
    struct stat buf;
    fstat(fileno(fp), &buf);
    return buf.st_size;
    }

static void fasta_split(FastaDB *fdb, gchar *output_stem,
                        gint num_chunks){
    register CompoundFile_Pos total
        = fasta_split_get_file_size(fdb->cf->fp),
        chunk_size = total/num_chunks;
    register CompoundFile_Pos *position = g_new(CompoundFile_Pos,
                                                num_chunks+1);
    register guint i;
    register gchar *chunk_path;
    for(i = 1; i < num_chunks; i++)
        position[i] = FastaDB_find_next_start(fdb, i*chunk_size);
    position[0] = 0;
    position[num_chunks] = total;
    for(i = 0; i < num_chunks; i++){
        chunk_path = g_strdup_printf("%s_chunk_%07d", output_stem, i);
        fasta_chunk_write(fdb, position[i],
                               position[i+1], chunk_path);
        g_free(chunk_path);
        }
    g_free(position);
    return;
    }

int Argument_main(Argument *arg){
    register FastaDB *fdb;
    register ArgumentSet *as
           = ArgumentSet_create("Sequence Input Options");
    gchar *query_path, *output_dir;
    gint num_chunks;
    register gchar *output_stem;
    ArgumentSet_add_option(as, 'f', "fasta", "path",
        "Fasta input file", NULL,
        Argument_parse_string, &query_path);
    ArgumentSet_add_option(as, 'o', "output", "dirpath",
        "Output directory", NULL,
        Argument_parse_string, &output_dir);
    ArgumentSet_add_option(as, 'c', "chunk", NULL,
        "Number of chunks to generate", "2",
        Argument_parse_int, &num_chunks);
    Argument_absorb_ArgumentSet(arg, as);
    Argument_process(arg, "fastasplit",
        "A utility to split fasta format sequences\n"
        "Guy St.C. Slater. guy@ebi.ac.uk. 2000-2003.\n", NULL);
    if(num_chunks < 1)
        g_error("Must request at least 1 chunk");
    fdb = FastaDB_open(query_path, NULL);
    output_stem = g_strdup_printf("%s%c%s",
            output_dir, G_DIR_SEPARATOR, basename(query_path));
    fasta_split(fdb, output_stem, num_chunks);
    FastaDB_close(fdb);
    return 0;
    }

