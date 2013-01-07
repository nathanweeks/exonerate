/****************************************************************\
*                                                                *
*  Library for manipulation of FASTA format databases            *
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

#include "fastadb.h"

gint Argument_main(Argument *arg){
    register FastaDB *fdb;
    register FastaDB_Seq *fdbs;
    register gchar *path;
    register Alphabet *alphabet;
    if(arg->argc == 2){
        path = arg->argv[1];
        alphabet = Alphabet_create(Alphabet_Type_UNKNOWN, FALSE);
        fdb = FastaDB_open(path, alphabet);
        Alphabet_destroy(alphabet);
        g_message("opened file [%s]", path);
        while((fdbs = FastaDB_next(fdb, FastaDB_Mask_ALL))){
            g_message("read [%s]", fdbs->seq->id);
            FastaDB_Seq_print(fdbs, stdout, FastaDB_Mask_ALL);
            FastaDB_Seq_destroy(fdbs);
            }
        FastaDB_close(fdb);
        g_message("closed file [%s]", path);
        }
    return 0;
    }

