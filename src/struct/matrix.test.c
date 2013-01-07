/****************************************************************\
*                                                                *
*  Simple matrix creation routines                               *
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

#include "matrix.h"

int main(void){
    register gint **m2d, ***m3d;
    register gshort ***ms3d, ****ms4d;
    register gint i, j, k, l;
    g_message("testing 2d matrix routines\n");
    m2d = (gint**)Matrix2d_create(5, 10, sizeof(gint));
    for(i = 0; i < 5; i++){
        for(j = 0; j < 10; j++){
            g_assert(!m2d[i][j]);
            m2d[i][j] = i+j;
            }
        }
    for(i = 0; i < 5; i++){
        for(j = 0; j < 10; j++){
            g_print("[%2d]", m2d[i][j]);
            }
        g_print("\n");
        }
    g_print("\n");
    g_free(m2d);
    /**/
    g_print("\n");
    g_message("testing 3d matrix routines\n");
    m3d = (gint***)Matrix3d_create(4, 3, 2, sizeof(gint));
    for(i = 0; i < 4; i++){
        for(j = 0; j < 3; j++){
            for(k = 0; k < 2; k++){
                g_assert(!m3d[i][j][k]);
                m3d[i][j][k] = i+j+k;
                }
            }
        }
    for(i = 0; i < 4; i++){
        for(j = 0; j < 3; j++){
            g_print(" {");
            for(k = 0; k < 2; k++){
                g_print("[%2d]", m3d[i][j][k]);
                }
            g_print("}");
            }
        g_print("\n");
        }
    g_print("\n");
    g_free(m3d);
    /**/
    g_message("testing 3d matrix routines on gshorts\n");
    ms3d = (gshort***)Matrix3d_create(4, 3, 2, sizeof(gshort));
    for(i = 0; i < 4; i++){
        for(j = 0; j < 3; j++){
            for(k = 0; k < 2; k++){
                g_assert(!ms3d[i][j][k]);
                ms3d[i][j][k] = i+j+k;
                }
            }
        }
    for(i = 0; i < 4; i++){
        for(j = 0; j < 3; j++){
            g_print(" {");
            for(k = 0; k < 2; k++){
                g_print("[%2d]", ms3d[i][j][k]);
                }
            g_print("}");
            }
        g_print("\n");
        }
    g_print("\n");
    g_free(ms3d);
    /**/
    g_message("testing 4d matrix routines on gshorts\n");
    ms4d = (gshort****)Matrix4d_create(5, 4, 3, 2, sizeof(gshort));
    for(i = 0; i < 5; i++){
        for(j = 0; j < 4; j++){
            for(k = 0; k < 3; k++){
                for(l = 0; l < 2; l++){
                    g_assert(!ms4d[i][j][k][l]);
                    ms4d[i][j][k][l] = i+j+k+l;
                    }
                }
            }
        }
    for(i = 0; i < 5; i++){
        for(j = 0; j < 4; j++){
            g_print(" {");
            for(k = 0; k < 3; k++){
                for(l = 0; l < 2; l++){
                    g_print("[%2d]", ms4d[i][j][k][l]);
                    }
                g_print(",");
                }
            g_print("}\n");
            }
        g_print("\n");
        }
    g_print("\n");
    g_free(ms4d);
    /**/
    return 0;
    }

