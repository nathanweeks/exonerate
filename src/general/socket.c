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
#include <string.h>    /* For strlen() */
#include <unistd.h>    /* For close() */
#include <sys/types.h>
#include <sys/socket.h>
#include <sys/time.h>
#include <netinet/tcp.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <netdb.h>
#include <fcntl.h>
#include <signal.h>   /* For sigaction() */
#include <sys/wait.h> /* For wait()   */
#include <errno.h>

#include "socket.h"

#define Socket_BUFSIZE BUFSIZ

static SocketConnection *SocketConnection_create(gchar *host,
                                                 gint  port){
    register SocketConnection *connection
     = g_new(SocketConnection, 1);
    register gint result;
    int flag = 1;
    if((connection->sock = socket(AF_INET, SOCK_STREAM, 0)) < 0){
        perror("opening stream socket");
        exit(1);
        }
    connection->host = g_strdup(host);
    connection->port = port;
    result = setsockopt(connection->sock, IPPROTO_TCP, TCP_NODELAY,
                                 (char*)&flag, sizeof(int));
    if(result < 0)
        perror("setting socket options");
    return connection;
    }

SocketClient *SocketClient_create(gchar *host, gint port){
    register SocketClient *client = g_new(SocketClient, 1);
    struct sockaddr_in server;
    struct hostent *hp = gethostbyname(host);
    register gchar *reply;
#ifdef USE_PTHREADS
    pthread_mutex_init(&client->connection_mutex, NULL);
#endif /* USE_PTHREADS */
    if(!hp){
        perror("looking up hostname");
        exit(1);
        }
    client->connection = SocketConnection_create(host, port);
#if 0
    /* Make non-blocking */
    if(fcntl(client->connection->sock, F_SETFL, O_NDELAY) == -1){
        perror("Tcp: Could not set O_NDELAY on socket with fcntl");
        exit(1);
        }
#endif /* 0 */
    server.sin_family = AF_INET;
    server.sin_port = htons(port);
    memmove(&server.sin_addr, hp->h_addr_list[0], hp->h_length);
    if(connect(client->connection->sock, (struct sockaddr*)&server,
               sizeof(server)) < 0){
        perror("connecting client socket");
        exit(1);
        }
    /* Send a packet to the server to test the connection.
     * If we get any reply, the connection is OK.
     *
     * FIXME: There should be a better way of doing this,
     *        but I can't find a method that works portably.
     */
    reply = SocketClient_send(client, "ping");
    if(reply){
        g_free(reply);
    } else {
        SocketClient_destroy(client);
        return NULL;
        }
    return client;
    }

static void SocketConnection_destroy(SocketConnection *connection){
    if(close(connection->sock) == -1){
        perror("closing socket");
        exit(1);
        }
    g_free(connection->host);
    g_free(connection);
    return;
    }

static gchar *SocketConnection_read(gint sock){
    register gint i, len = 0, line_complete = 0, line_expect = 1;
    register gchar *reply;
    register GString *string = g_string_sized_new(Socket_BUFSIZE);
    register gboolean line_count_given = FALSE;
    gchar buffer[Socket_BUFSIZE+1];
    do {
        if((len = recv(sock, buffer, Socket_BUFSIZE, 0)) <= 0){
            len = 0;
            break;
            }
        buffer[len] = '\0';
        for(i = 0; i < len; i++)
            if(buffer[i] == '\n')
                line_complete++;
        if((!string->len) && (len > 10))
            if(!strncmp(buffer, "linecount:", 10)){
                line_expect = atoi(&buffer[11]);
                line_count_given = TRUE;
                if(line_expect < 2)
                    g_error("linecount: must be > 1");
                }
        if(line_complete > line_expect){
            if(line_count_given){
                g_error("Received [%d] socket message lines, but expected [%d]",
                   line_complete, line_count_given);
            } else {
                g_error("Multiline socket messages must use linecount:");
                }
            }
        g_string_append(string, buffer);
    } while(line_complete < line_expect);
    if(string->len){
        reply = string->str;
        g_string_free(string, FALSE);
        return reply;
        }
    g_string_free(string, TRUE);
    return NULL;
    }

static void Socket_send_msg(gint sock, gchar *msg, gchar *err_msg){
    register gint start = 0, len, msglen = strlen(msg);
    do {
        if((len = send(sock, msg+start, msglen-start, 0)) < 0){
            perror(err_msg);
            exit(1);
            }
        start += len;
    } while(start < msglen);
    return;
    }

static void Socket_send(gint sock, gchar *msg, gchar *err_msg){
    register gint i, line_count = 0;
    register GString *full_msg = g_string_sized_new(strlen(msg)+32);
    for(i = 0; msg[i]; i++)
        if(msg[i] == '\n')
            line_count++;
    if(line_count > 1)
        g_string_sprintfa(full_msg, "linecount: %d\n", line_count+1);
    g_string_sprintfa(full_msg, "%s", msg);
    if(!line_count)
        g_string_append_c(full_msg, '\n');
    Socket_send_msg(sock, full_msg->str, err_msg);
    g_string_free(full_msg, TRUE);
    return;
    }

gchar *SocketClient_send(SocketClient *client, gchar *msg){
    register gchar *reply;
#ifdef USE_PTHREADS
    pthread_mutex_lock(&client->connection_mutex);
#endif /* USE_PTHREADS */
    Socket_send(client->connection->sock, msg, "writing client message");
    reply = SocketConnection_read(client->connection->sock);
    g_assert(reply);
#ifdef USE_PTHREADS
    pthread_mutex_unlock(&client->connection_mutex);
#endif /* USE_PTHREADS */
    return reply;
    }

void SocketClient_destroy(SocketClient *client){
    SocketConnection_destroy(client->connection);
#ifdef USE_PTHREADS
    pthread_mutex_destroy(&client->connection_mutex);
#endif /* USE_PTHREADS */
    g_free(client);
    return;
    }

static gint global_connection_count = 0;

static void SocketServer_reap_dead_children(int signum){
    while(waitpid(-1, NULL, WNOHANG) > 0)
        global_connection_count--;
    return;
    }

static void SocketServer_shutdown(int signum){
    g_message("Server shutting down");
    exit(0);
    return;
    }

#ifdef USE_PTHREADS
static void SocketServer_broken_pipe(int signum){
    const char error_message[] = "Server detected broken pipe - closing thread\n";
    write(STDERR_FILENO, error_message, sizeof(error_message));
    pthread_exit(NULL);
    return;
    }
#endif

SocketServer *SocketServer_create(gint port, gint max_connections,
              SocketProcessFunc server_process_func,
              SocketConnectionOpenFunc connection_open_func,
              SocketConnectionCloseFunc connection_close_func,
              gpointer user_data){
    register SocketServer *server = g_new(SocketServer, 1);
    struct sockaddr_in sock_server;
    struct sigaction sa;
    socklen_t len = sizeof(sock_server);
    server->connection = SocketConnection_create("localhost", port);
    g_assert(server_process_func);
    server->server_process_func = server_process_func;
    server->connection_open_func = connection_open_func;
    server->connection_close_func = connection_close_func;
    server->user_data = user_data;
    server->max_connections = max_connections;
    sock_server.sin_family = AF_INET;
    sock_server.sin_addr.s_addr = INADDR_ANY;
    sock_server.sin_port = htons(port);
    if(bind(server->connection->sock,
            (struct sockaddr*)&sock_server, len)){
        perror("binding stream socket");
        exit(1);
        }
    if(getsockname(server->connection->sock,
                   (struct sockaddr*)&sock_server, &len)){
        perror("getting socket name");
        exit(1);
        }
    server->connection->port = ntohs(sock_server.sin_port);
    sa.sa_handler = SocketServer_reap_dead_children;
    sigemptyset(&sa.sa_mask);
    sa.sa_flags = SA_RESTART;
    if (sigaction(SIGCHLD, &sa, NULL) == -1) {
        perror("sigaction");
        exit(1);
    }
    sa.sa_handler = SocketServer_shutdown;
    sa.sa_flags = 0;
    if (sigaction(SIGTERM, &sa, NULL) == -1) {
        perror("sigaction");
        exit(1);
    }
#ifdef USE_PTHREADS
    server->sspd = g_new0(SocketServer_pthread_Data, max_connections);
    pthread_mutex_init(&server->connection_mutex, NULL);
#endif /* USE_PTHREADS */
    return server;
    }

static void SocketServer_process_connection(SocketServer *server, int msgsock){
    register gpointer connection_data = server->connection_open_func
                    ? server->connection_open_func(server->user_data)
                    : NULL;
    register gchar *msg;
    register gboolean ok = TRUE;
    gchar *reply;
    do {
#ifdef USE_PTHREADS
        pthread_mutex_lock(&server->connection_mutex);
#endif /* USE_PTHREADS */
        msg = SocketConnection_read(msgsock);
#ifdef USE_PTHREADS
        pthread_mutex_unlock(&server->connection_mutex);
#endif /* USE_PTHREADS */
        reply = NULL;
        if(!msg)
            break;
        ok = server->server_process_func(msg, &reply, connection_data,
                                         server->user_data);
        g_free(msg);
        if(reply){
#ifdef USE_PTHREADS
            pthread_mutex_lock(&server->connection_mutex);
#endif /* USE_PTHREADS */
            Socket_send(msgsock, reply, "writing reply");
#ifdef USE_PTHREADS
            pthread_mutex_unlock(&server->connection_mutex);
#endif /* USE_PTHREADS */
            g_free(reply);
        } else {
            g_error("no reply from server");
            }
    } while(ok);
    if(server->connection_close_func)
        server->connection_close_func(connection_data, server->user_data);
    return;
    }

#ifdef USE_PTHREADS
static void *SocketServer_pthread_func(void* data){
    register SocketServer_pthread_Data *sspd = (SocketServer_pthread_Data*)data;
    signal(SIGPIPE, SocketServer_broken_pipe);
    SocketServer_process_connection(sspd->server, sspd->msgsock);
    pthread_mutex_lock(&sspd->server->connection_mutex);
    g_message("cleaning up connection [%d]", global_connection_count);
    global_connection_count--;
    close(sspd->msgsock);
    sspd->in_use = FALSE;
    pthread_mutex_unlock(&sspd->server->connection_mutex);
    pthread_exit(NULL);
    return NULL;
    }
#endif /* USE_PTHREADS */

gboolean SocketServer_listen(SocketServer *server){
    register int msgsock;
    struct sockaddr_in client_addr;
    socklen_t client_len = sizeof(struct sockaddr_in);
#ifdef USE_PTHREADS
    register gint i;
    pthread_attr_t pt_attr;
    pthread_attr_init(&pt_attr);
    pthread_attr_setdetachstate(&pt_attr, PTHREAD_CREATE_DETACHED);
#endif /* USE_PTHREADS */
    /**/
    listen(server->connection->sock, server->max_connections);
    while ((msgsock = accept(server->connection->sock,
                     (struct sockaddr*)&client_addr, &client_len)) == -1) {
        if(errno == EINTR)
            continue;
        perror("server accept");
        close(msgsock);
        exit(1);
        }
    /* FIXME: send an busy warning message back to the client ? */
#ifdef USE_PTHREADS
    pthread_mutex_lock(&server->connection_mutex);
    if(global_connection_count >= server->max_connections){
        pthread_mutex_unlock(&server->connection_mutex);
        g_message("Max connections reached");
    } else {
        for(i = 0; i < server->max_connections; i++)
            if(!server->sspd[i].in_use){
                server->sspd[i].in_use = TRUE;
                break;
                }
        server->sspd[i].server = server;
        server->sspd[i].msgsock = dup(msgsock);
        if(server->sspd[i].msgsock < 0)
            perror("duplicating socket for pthread");
        global_connection_count++;
        g_message("opened connection [%d/%d] from [%s]",
                   global_connection_count,
                   server->max_connections,
                   inet_ntoa(client_addr.sin_addr));
        pthread_mutex_unlock(&server->connection_mutex);
        pthread_create(&server->sspd[i].thread, &pt_attr,
                       SocketServer_pthread_func, (void*)&server->sspd[i]);
        }
#else /* USE_PTHREADS */
    if(global_connection_count >= server->max_connections){
        g_message("Max connections reached");
    } else {
        if(fork() == 0){
            SocketServer_process_connection(server, msgsock);
            exit(0);
        } else {
            global_connection_count++;
            g_message("opened connection [%d/%d] from [%s]",
                    global_connection_count,
                    server->max_connections,
                    inet_ntoa(client_addr.sin_addr));
            }
        }
#endif /* USE_PTHREADS */
    close(msgsock);
    return TRUE; /* Keep server running */
    }

void SocketServer_destroy(SocketServer *server){
    SocketConnection_destroy(server->connection);
#ifdef USE_PTHREADS
    g_free(server->sspd);
    pthread_mutex_destroy(&server->connection_mutex);
#endif /* USE_PTHREADS */
    g_free(server);
    return;
    }


