#include <iostream>
#include <string>
#include <sys/socket.h>
#include <netinet/in.h>
#include <unistd.h>
#include <arpa/inet.h>
#include "SocketComm.hpp"
#include "AliceProcessor.hpp"

// Main entry point of the program.
int main(int argc, char* argv[]) {
    // Check if the correct number of arguments is provided.
    if(argc != 7) {
        // If not, print the correct usage and exit with an error code.
        std::cerr << "Usage: " << argv[0] << " <root folder> <filename> <log filename> <paillier key size> <port number> <file index>" << std::endl;
        return 1;
    }

    // Parse command line arguments.
    std::string rootFolder = argv[1]; // Root folder name.
    std::string filenamePart = argv[2]; // The main part of dataset filename such as "scanpaths_Salient360_Alice".
    std::string logFilenamePart = argv[3]; // Output log filename.
    int paillierKeySize = std::stoi(argv[4]); // Number of bits for the Paillier key for encryption.
    int portNumber = std::stoi(argv[5]); // Port number for socket communication.
    int fileIndex = std::stoi(argv[6]); // Index of the file to be processed - person id.

    // Display the GMP (GNU Multiple Precision Arithmetic Library) version.
    printf("GMP Version: %s\n", gmp_version);

    // Initialize socket communication and the main processor.
    SocketComm aliceSocket;
    AliceProcessor aliceProcessor(aliceSocket, rootFolder, filenamePart, logFilenamePart, paillierKeySize, portNumber, fileIndex);

    // Start the communication with Bob.
    aliceProcessor.startProcess();

    // Close the socket after processing is complete.
    aliceSocket.close_socket();

    // Return 0 to indicate successful execution.
    return 0;
}
