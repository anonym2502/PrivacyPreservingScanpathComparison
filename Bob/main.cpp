#include "SocketComm.hpp"
#include <iostream>
#include "BobProcessor.hpp"
#include <iostream>
#include <vector>
#include <set>
#include <random>
#include <iterator>
#include <iostream>
#include <vector>
#include <set>
#include <ctime>
#include <cstdlib>
#include <map>
#include <iostream>
#include <vector>
#include <map>
#include <fstream>
#include <ctime>
#include <set>
#include <cstdlib>
#include <cmath>
#include <numeric>




int main(int argc, char* argv[]) {

    if(argc != 6) {  // Adjust this check based on the total number of arguments your program requires
        std::cerr << "Usage: " << argv[0] << " <root folder> <filename> <log filename> <port number> <fileinxex>" << std::endl;
        return 1;
    }
    std::string rootFolder = argv[1];
    std::string filename = argv[2];
    std::string logFilename = argv[3];
    int portNumber = std::stoi(argv[4]);  // Assuming port number is the 6th argument
    int fileIndex = std::stoi(argv[5]);  // Convert string to integer

    SocketComm bobSocket;

    BobProcessor bobProcessor(bobSocket, rootFolder, filename, logFilename, portNumber,fileIndex);
    bobProcessor.startProcess();


    bobSocket.close_socket();

    return 0;
}


