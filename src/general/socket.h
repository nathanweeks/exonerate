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

#ifndef INCLUDED_SOCKET_H
#define INCLUDED_SOCKET_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <glib.h>

#ifdef USE_PTHREADS
#include <pthread.h>
#endif /* USE_PTHREADS */

typedef gboolean SocketProcessFunc(gchar *msg, gchar **reply,
                                   gpointer connection_data,
                                   gpointer user_data);
/* Return FALSE to close the connection
 * If reply is set, it will be freed by the server.
 */

typedef gpointer SocketConnectionOpenFunc(gpointer user_data);
typedef void SocketConnectionCloseFunc(gpointer connection_data,
                                       gpointer user_data);

typedef struct {
      int  sock;
    gchar *host;
     gint  port;
} SocketConnection;

typedef struct {
    SocketConnection *connection;
#ifdef USE_PTHREADS
     pthread_mutex_t  connection_mutex;
#endif /* USE_PTHREADS */
} SocketClient;

#ifdef USE_PTHREADS
typedef struct {
               gboolean  in_use;
    struct SocketServer *server;
                    int  msgsock;
              pthread_t  thread;
} SocketServer_pthread_Data;
#endif /* USE_PTHREADS */

typedef struct SocketServer {
            SocketConnection *connection;
           SocketProcessFunc *server_process_func;
    SocketConnectionOpenFunc *connection_open_func;
   SocketConnectionCloseFunc *connection_close_func;
                    gpointer  user_data;
                        gint  max_connections;
#ifdef USE_PTHREADS
   SocketServer_pthread_Data *sspd;
             pthread_mutex_t  connection_mutex;
#endif /* USE_PTHREADS */
} SocketServer;

SocketClient *SocketClient_create(gchar *host, gint port);
       gchar *SocketClient_send(SocketClient *client, gchar *msg);
        void  SocketClient_destroy(SocketClient *client);

SocketServer *SocketServer_create(gint port, gint max_connections,
                           SocketProcessFunc server_process_func,
                           SocketConnectionOpenFunc connection_open_func,
                           SocketConnectionCloseFunc connection_close_func,
                           gpointer user_data);
    gboolean  SocketServer_listen(SocketServer *server);
        void  SocketServer_destroy(SocketServer *server);


/**/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* INCLUDED_SOCKET_H */

