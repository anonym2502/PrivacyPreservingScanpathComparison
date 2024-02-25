#ifndef SOCKET_COMM_HPP
#define SOCKET_COMM_HPP

#include <string>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <unistd.h>

class SocketComm {
private:
    int sockfd;
    struct sockaddr_in server_addr, client_addr;
    bool is_server;

public:
    SocketComm();
    ~SocketComm();

    bool init_as_server(int port);
    bool init_as_client(const std::string& ip, int port);
    bool listen_for_client();
    bool send_message(const std::string& msg);
    std::string receive_message();
    void close_socket();
};

#endif // SOCKET_COMM_HPP
