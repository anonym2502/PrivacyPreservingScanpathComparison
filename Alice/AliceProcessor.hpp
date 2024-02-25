#ifndef ALICE_PROCESSOR_HPP
#define ALICE_PROCESSOR_HPP

#include "SocketComm.hpp"
#include "Paillier.hpp"
#include "utils.hpp"
#include <string>
#include <vector>
#include <unordered_map>
#include <gmp.h>

class AliceProcessor {
private:
    SocketComm socketComm;  // Handles socket communication
    Paillier paillier;      // Paillier cryptosystem for encryption/decryptio
    std::string bobScanpathNumberMessage;
    std::string rootFolder; // Root folder path for file operations
    std::string filename;   // Path to the file containing scanpaths
    std::string logFilePath; // File path for logging
    int fileIndex;          // Index of the file being processed
    int paillierKeySize;    // Size of the Paillier key
    int portNumber;         // Network port number for communication

public:
    mpz_t public_n; // Public component 'n' of the Paillier key
    mpz_t public_g; // Public component 'g' of the Paillier key

    // Statistical metrics for various operations
    OnlineStats statsMinCalculation;
    OnlineStats statsDistanceMatrixCalculation;
    OnlineStats statsInitialization;
    OnlineStats statsIndividualScanpathComparison;
    OnlineStats statsUpdateKey;

    // Map to store encrypted matrices associated with Alice's scanpaths
    std::unordered_map<std::string, std::vector<MPZ_Wrapper>> encryptedMatrixMap;

    int iterationCounter; // Counter for the number of iterations in the process

    // Constructor
    AliceProcessor(const SocketComm& socketCommObj, const std::string& rootFolder,
                   const std::string& filenamePart, const std::string& logFilenamePart,
                   int paillierKeySize, int portNumber, int index);

    // Processes the input string based on operation flags and generates the result string
    bool process(std::string& resultString, const std::string& inputString);

    // Handles the operation when the flag is zero (initial substitution cost matrix generation)
    bool handleFlagZero(const std::string& inputString, std::string& resultString);

    // Handles the operation when the flag is one (minimum computation of encrypted values)
    bool handleFlagOne(const std::string& inputString, std::string& resultString);

    // Handles the operation when the flag is two (decryption of final similarity cost)
    bool handleFlagTwo(const std::string& inputString, std::string& resultString);

    // Initiates the scanpath comparison process and handles all processes for scanpath comparison
    void startProcess();

    // Logs the accumulated statistics of the scanpath comparison process
    void logStatistics();

    // Creates a vector of public key components for the Paillier cryptosystem
    std::vector<MPZ_Wrapper> createPublicKeyVector();

    // Updates the Paillier key used for encryption/decryption
    void updatePaillierKey();

    // Resets statistical metrics for a new iteration in the scanpath comparison process
    void resetStatistics();

    // Logs the completion details of an iteration, including time taken
    void logIterationCompletion(Timer timerStartLoop);

    // Logs the start of an iteration, including details about the current Alice scanpath
    void logIterationStart(std::string aliceScanpath);
};

#endif // ALICE_PROCESSOR_HPP
