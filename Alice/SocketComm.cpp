#include "SocketComm.hpp"
#include <netinet/tcp.h>  // Required for TCP_NODELAY
#include <cstring>


SocketComm::SocketComm() : sockfd(-1), is_server(false) {
    memset(&server_addr, 0, sizeof(server_addr));
    memset(&client_addr, 0, sizeof(client_addr));
}

SocketComm::~SocketComm() {
    if (sockfd != -1) {
        close(sockfd);
    }
}

bool SocketComm::init_as_server(int port) {
    sockfd = socket(AF_INET, SOCK_STREAM, 0);
    
    if (sockfd < 0) return false;

    int optval = 1;
    setsockopt(sockfd, IPPROTO_TCP, TCP_NODELAY, &optval, sizeof(optval));

    server_addr.sin_family = AF_INET;
    server_addr.sin_addr.s_addr = INADDR_ANY;
    
    server_addr.sin_port = htons(port);

    if (bind(sockfd, (struct sockaddr*)&server_addr, sizeof(server_addr)) < 0) {
        return false;
    }

    is_server = true;
    return true;
}

bool SocketComm::init_as_client(const std::string& ip, int port) {
    sockfd = socket(AF_INET, SOCK_STREAM, 0);
    if (sockfd < 0) return false;

    int optval = 1;
    setsockopt(sockfd, IPPROTO_TCP, TCP_NODELAY, &optval, sizeof(optval));

    server_addr.sin_family = AF_INET;
    server_addr.sin_port = htons(port);
    server_addr.sin_addr.s_addr = inet_addr(ip.c_str());

    if (connect(sockfd, (struct sockaddr*)&server_addr, sizeof(server_addr)) < 0) {
        return false;
    }

    is_server = false;
    return true;
}

bool SocketComm::listen_for_client() {
    if (!is_server) return false;

    listen(sockfd, 5);
    socklen_t clilen = sizeof(client_addr);
    int new_sockfd = accept(sockfd, (struct sockaddr*)&client_addr, &clilen);

    if (new_sockfd < 0) {
        return false;
    }

    close(sockfd);
    sockfd = new_sockfd;

    return true;
}

bool SocketComm::send_message(const std::string& msg) {
    // Send the message length first
    uint32_t msg_length = htonl(msg.length());  // Convert to network byte order
    send(sockfd, &msg_length, sizeof(msg_length), 0);

    // Send the actual message in chunks
    size_t bytesSent = 0;
    while(bytesSent < msg.length()) {
        int n = send(sockfd, msg.c_str() + bytesSent, msg.length() - bytesSent, 0);
        if(n <= 0) return false;
        bytesSent += n;
    }
    return true;
}

std::string SocketComm::receive_message() {
    // Receive the length of the message first
    uint32_t msg_length = 0;
    int n = recv(sockfd, &msg_length, sizeof(msg_length), 0);
    if (n <= 0) return "";
    msg_length = ntohl(msg_length);  // Convert from network byte order

    // Receive the actual message in chunks
    std::string msg;
    char buffer[256];
    size_t bytesReceived = 0;
    while(bytesReceived < msg_length) {
        int n = recv(sockfd, buffer, std::min(sizeof(buffer), msg_length - bytesReceived), 0);
        if(n <= 0) return "";
        msg.append(buffer, n);
        bytesReceived += n;
    }
    return msg;
}
void SocketComm::close_socket() {
    close(sockfd);
}
