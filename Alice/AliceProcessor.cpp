#include "AliceProcessor.hpp"
#include <iostream>
#include "Paillier.hpp"
#include "utils.hpp"
#include <string>
#include <chrono>
#include <iomanip>
#include <ctime>
#include <fstream>


// Constructor for AliceProcessor class
AliceProcessor::AliceProcessor(const SocketComm& socketCommObj, const std::string& rootFolder,
                               const std::string& filenamePart, const std::string& logFilenamePart,
                               int paillierKeySize, int portNumber, int index)
    : socketComm(socketCommObj), paillier(paillierKeySize), rootFolder(rootFolder), portNumber(portNumber), fileIndex(index)
{
    // Initialize the file paths for datasets and logs
    filename = rootFolder + "Datasets/" + filenamePart + "_" + std::to_string(fileIndex) + ".txt";
    logFilePath = rootFolder + "Output/" + logFilenamePart;

    // Initialize public key components
    this->paillierKeySize = paillierKeySize;
    mpz_init(public_n);
    mpz_init(public_g);
}


// This function initiates the scanpath comparison process and handles all processes of the scanpath comparison.
void AliceProcessor::startProcess(){
    // Variable declarations
    std::string outputAlice, inputAlice;

    // Start Alice server and wait for connection
    startAliceServerAndListen(socketComm, portNumber);

    // Load scanpaths and send the size message to Bob
    std::vector<std::string> allScanpaths = loadScanpathsFromFile(filename);
    sendSizeMessage(allScanpaths, socketComm, 3);

    // Logging the start of communication
    log("Communication", "Receiving number of Bob scanpaths.", logFilePath, "Alice", iterationCounter, "Start");

    // Receive a message from Bob
    bobScanpathNumberMessage = socketComm.receive_message();
        
    // Process each scanpath
    iterationCounter = 0;
    for (const std::string& aliceScanpath : allScanpaths) {
        
        /*
         // Uncomment the following block if you need to update the public key for different Alice scanpaths.      updatePaillerKey();
         */
        
        // Iterate over the number of scanpaths received from Bob
        for (int k = 0; k < getNumberOfScanpaths(bobScanpathNumberMessage); ++k) {
            // Reset statistical counters for the new iteration
            resetStatistics();
            
            // Log the start of this iteration, including details about the current Alice scanpath
            logIterationStart(aliceScanpath);
            
            // Start the timer to measure the duration of processing the current scanpath
            Timer timerStartLoop;
            
            // Set the input for processing as the current Alice scanpath and add the operation flag for initialization.
            inputAlice = aliceScanpath;
            addOperationFlag(0, inputAlice);
            
            // Process the current scanpath. Initialization process with flag 0.
            process(outputAlice, inputAlice);
            
            // Update the initialization statistics with the duration from the timer
            statsInitialization.update(timerStartLoop.getDuration());
            
            // Loop to handle continuous processing based on incoming messages
            bool controlWhile = true;
            while(controlWhile){
                // Receive messages from Bob and process them for minimum computation.
                inputAlice = socketComm.receive_message();
                controlWhile = process(outputAlice, inputAlice);
            }
            
            // Log the completion of the iteration
            logIterationCompletion(timerStartLoop);
            
            // Increment the counter for the next iteration
            iterationCounter = iterationCounter + 1;
        }
        }
    // Log final statistics after all processing is complete.
    logStatistics();
}



// Process input string based on operation flag
bool AliceProcessor::process(std::string& resultString, const std::string& inputString) {
    int operation_flag = getOperationFlag(inputString);

    // Handle different operations based on the flag
    switch (operation_flag) {
        // Initial substitution cost matrix generation and send the encryted matrix to Bob.
        case 0:
            return handleFlagZero(inputString, resultString);
        //It is to find minimum for encrypted edit sequence costs. It is called in each iteration.
        case 1:
            return handleFlagOne(inputString, resultString);
        // It is to decrypt final similarity cost and share it with Bob.
        case 2:
            return handleFlagTwo(inputString, resultString);
        default:
            // Handle unknown operation flag
            return false;
    }
}



// This function is compute substitution cost matrix for Alice scanpath and send it to Bob.
bool AliceProcessor::handleFlagZero(const std::string& inputString, std::string& resultString) {
    // Extract Alice's scanpath from the input string
    std::string aliceScanpath = inputString.substr(2);


    // Start the timer for measuring the encryption duration
    Timer timerEncryptedDistanceMatrix;

    // Vector to store the encrypted distance matrix
    std::vector<MPZ_Wrapper> encryptedDistanceMatrix;

    // Uncomment the following block if you would like to store substitution cost matrix to make it more efficient, you can use following block. However, you should not update the Paillier key for same scanpath.
    /*
    // Check if the encrypted distance matrix is already calculated for this scanpath
    bool scanpathFoundControl = getEncryptedDistanceMatrix(aliceScanpath, encryptedMatrixMap, encryptedDistanceMatrix);

    // If not already calculated, compute and store it
    if (!scanpathFoundControl) {
        encryptedDistanceMatrix = encryptAndFlattenDifferenceMatrix(this->paillier, aliceScanpath);
        encryptedMatrixMap[aliceScanpath] = encryptedDistanceMatrix;
    }
    */
    
    // Compute substitution cost matrix for Alice scanpath
    encryptedDistanceMatrix = encryptAndFlattenDifferenceMatrix(this->paillier, aliceScanpath);


    // Create a vector containing the public key components
    std::vector<MPZ_Wrapper> publicKeyVector = createPublicKeyVector();
    // Combine the public key vector with the encrypted distance matrix
    std::vector<MPZ_Wrapper> response0 = appendVectors(publicKeyVector, encryptedDistanceMatrix);

    // Update statistics with the encryption duration
    statsDistanceMatrixCalculation.update(timerEncryptedDistanceMatrix.getDuration());

    // Send the combined vector to Bob
    sendEncryptedResult(socketComm, response0, 0);
    
    encryptedDistanceMatrix.clear();

    return true;
}


// This function handles minimum computation of the values sent by Bob.
bool AliceProcessor::handleFlagOne(const std::string& inputString, std::string& resultString) {
    // Start the timer for min value calculation
    Timer timerMinCalculation;

    // Convert the input string to a vector of MPZ_Wrapper objects
    std::vector<MPZ_Wrapper> inputValues = convertStringToMPZ_Wrapper(inputString);

    // Initialize variables for calculation
    mpz_t minValue, tempPlaintext, encryptedMinValue;
    mpz_inits(minValue, tempPlaintext, encryptedMinValue, NULL);

    // Set minValue to a large number (public_n)
    mpz_set(minValue, public_n);

    // Find the minimum value among decrypted input values
    for (size_t i = 1; i < inputValues.size(); ++i) {
        this->paillier.decryptMessage(tempPlaintext, inputValues[i].value);
        if (mpz_cmp(tempPlaintext, minValue) < 0) {
            mpz_set(minValue, tempPlaintext);
        }
    }

    // Encrypt the min value
    paillier.encryptMessage(encryptedMinValue, minValue);

    // Push the encrypted min value into a vector
    std::vector<MPZ_Wrapper> minValueVector;
    minValueVector.push_back(MPZ_Wrapper(encryptedMinValue));

    // Clear the used variables
    mpz_clears(encryptedMinValue, minValue, tempPlaintext, NULL);

    // Send the encrypted min value vector as a result
    sendEncryptedResult(socketComm, minValueVector, 1);

    return true;
}

// This function is used at the end. It gets the encrypted final similarity value from Bob and then send the decrypted one.
bool AliceProcessor::handleFlagTwo(const std::string& inputString, std::string& resultString) {
    // Initialize a variable for plaintext
    mpz_t tempPlaintext;
    mpz_init(tempPlaintext);

    // Convert the input string to a vector of MPZ_Wrapper objects
    std::vector<MPZ_Wrapper> inputValues = convertStringToMPZ_Wrapper(inputString);

    // Decrypt the second value in the input values
    this->paillier.decryptMessage(tempPlaintext, inputValues[1].value);

    // Print the decrypted final results
    gmp_printf("Decrypted Final Results: %Zd\n", tempPlaintext);

    // Push the decrypted value into a result vector
    std::vector<MPZ_Wrapper> resultVector;
    resultVector.push_back(MPZ_Wrapper(tempPlaintext));

    // Send the result vector
    sendEncryptedResult(socketComm, resultVector, 2);

    // Clear the plaintext variable
    mpz_clear(tempPlaintext);

    return false;
}

void AliceProcessor::logStatistics() {
    // Log overall statistics for  scanpath comparison
    log("Scanpath comparison is done.",
        buildLogMessage("Current statistics for individual scanpath comparison: ", statsIndividualScanpathComparison),
        logFilePath, "Alice", iterationCounter, "Value");

    // Log statistics specific to initialization processes
    log("Scanpath comparison is done.",
        buildLogMessage("Statistics for initialization: ", statsInitialization),
        logFilePath, "Alice", iterationCounter, "Value");

    // Log statistics for distance matrix(encrypted substitution cost matrix for Alice's scanpath) calculations
    log("Scanpath comparison is done.",
        buildLogMessage("Current statistics for distance matrix calculation: ", statsDistanceMatrixCalculation),
        logFilePath, "Alice", iterationCounter, "Value");

    // Log statistics for each minimum value calculation in the process
    log("Scanpath comparison is done.",
        buildLogMessage("Current statistics for each min calculation(each small step): ", statsMinCalculation),
        logFilePath, "Alice", iterationCounter, "Value");

    // Log statistics for the updates on the Paillier key
    log("Scanpath comparison is done.",
        buildLogMessage("Current statistics for Stats Update Key (For each Alice's Scanpath): ", statsUpdateKey),
        logFilePath, "Alice", iterationCounter, "Value");

    // Flush all remaining logs to the log file
    flush_remaining_logs(logFilePath);

    // Calculate and output the average time taken for each iteration
    double averageTime = static_cast<double>(statsIndividualScanpathComparison.getMean());
    std::cout << "Average time for each iteration: " << averageTime << " microseconds." << std::endl;
}




std::vector<MPZ_Wrapper> AliceProcessor::createPublicKeyVector() {
    // Get the public key components from the Paillier cryptosystem
    this->paillier.getPublicKey(this->public_n, this->public_g);

    // Create a vector to hold the public key components
    std::vector<MPZ_Wrapper> publicKeyVector;
    publicKeyVector.push_back(MPZ_Wrapper(this->public_n));
    publicKeyVector.push_back(MPZ_Wrapper(this->public_g));

    // Return the vector containing the public key components
    return publicKeyVector;
}

void AliceProcessor::updatePaillierKey() {
    // Start a timer to measure the key update duration
    Timer timerUpdateKey;

    // Reset the statistics for the key update
    statsUpdateKey.reset();

    // Update the Paillier key
    paillier.updateEncryptionKey(paillierKeySize);

    // Retrieve the updated public key components
    paillier.getPublicKey(public_n, public_g);

    // Update the statistics with the duration of the key update
    statsUpdateKey.update(timerUpdateKey.getDuration());
}


void AliceProcessor::resetStatistics() {
    // Reset statistical metrics for various calculations and comparisons in each scanpath comparison.
    statsMinCalculation.reset();
    statsDistanceMatrixCalculation.reset();
    statsInitialization.reset();
    statsIndividualScanpathComparison.reset();
}
void AliceProcessor::logIterationStart(std::string aliceScanpath) {
    // Log the start of a new loop iteration with Alice's scanpath details
    log("Loop", "Loop.", logFilePath, "Alice", iterationCounter, "Start");
    log("Scanpath", aliceScanpath, logFilePath, "Alice", iterationCounter, "Value");
    log("Scanpath Length", std::to_string(aliceScanpath.length()), logFilePath, "Alice", iterationCounter, "Value");
}



void AliceProcessor::logIterationCompletion(Timer timerStartLoop) {
    // Logging at the end of a loop
    // Calculate the duration of the iteration
    std::chrono::duration<double> iterationTime = timerStartLoop.getDuration();

    // Update statistics with the iteration duration
    statsIndividualScanpathComparison.update(iterationTime);

    // Log the completion of the loop iteration and the time it took
    log("Loop", "Loop.", logFilePath, "Alice", iterationCounter, "End");
    log("Loop", "Loop ended. It took " + durationToString(iterationTime, 5) + " microseconds.", logFilePath, "Bob", iterationCounter, "End");

    // Log the overall statistics at the end of the iteration
    logStatistics();
}


