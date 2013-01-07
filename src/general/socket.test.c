/****************************************************************\
*                                                                *
*  Simple client-server code library                             *
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <glib.h>

#include "socket.h"

static gchar *read_next_line(gchar *prompt){
    gchar buffer[1024];
    register gchar *msg = NULL;
    g_print("%s", prompt);
    if(!fgets(buffer, 1024, stdin))
          g_error("Didn't read anything");
    if(buffer[0] != '\n'){
        msg = g_strdup(buffer);
        g_message("ret [%s]", msg);
        }
    return msg;
    }

static void run_client(gchar *host, gint port){
    register SocketClient *sc = SocketClient_create(host, port);
    register gchar *msg, *reply;
    while(TRUE){
        msg = read_next_line("> ");
        g_message("client going to send [%s]", msg);
        if(!msg)
            continue;
        if(!strncmp(msg, "quit client", 11))
            break;
        reply = SocketClient_send(sc, msg);
        if(reply){
            g_print(" reply [%s]", reply);
            g_free(reply);
            }
        g_free(msg);
        }
    SocketClient_destroy(sc);
    return;
    }

static gboolean test_server_func(gchar *msg, gchar **reply,
                                 gpointer connection_data,
                                 gpointer user_data){
    g_message("server received [%s]", msg);
    (*reply) = g_strdup_printf("msg received[%s]", msg);
    return strncmp(msg, "close", 11)?FALSE:TRUE;
    }

static void run_server(gint port){
    register SocketServer *ss = SocketServer_create(port, 2,
                        test_server_func, NULL, NULL, NULL);
    while(SocketServer_listen(ss));
    SocketServer_destroy(ss);
    return;
    }

int main(int argc, char **argv){
    register gboolean be_client;
    register gint port;
    register gchar *host;
    if(argc != 4){
        g_warning("Usage: socket.test <c|s> host port");
        return 0; /* exit quietly to please make check */
        }
    be_client = (argv[1][0] == 'c');
    host = argv[2];
    port = atoi(argv[3]);
    g_message("client [%s] host [%s] port [%d]",
              be_client?"yes":"no", host, port);
    if(be_client)
        run_client(host, port);
    else
        run_server(port);
    return 0;
    }

