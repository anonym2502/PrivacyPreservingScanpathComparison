#include "BobProcessor.hpp"
#include "PublicPaillier.hpp"
#include "utils.hpp"
#include <string>
#include <chrono>
#include <iomanip>
#include <ctime>
#include <fstream>

// This function initiates the processing on Bob's side.
void BobProcessor::startProcess() {
    iterationCounter = 0;  // Initialize iteration counter

    alphabet = createAlphabetMapping();  // Create a mapping for an alphabet used in processing
    
    initializeConnection();  // Set up a network connection
    std::string aliceScanpathNumberMessage = socketComm.receive_message();  // Receive a message from Alice
    allScanpaths = loadScanpathsFromFile(filename);  // Load scanpaths from a file
    std::string sizeMessage = std::to_string(allScanpaths.size());  // Convert the size of allScanpaths to a string
    addOperationFlag(3, sizeMessage);  // Add an operation flag to the size message
    // Sending the size message to Alice
    if (!socketComm.send_message(sizeMessage)) {
        std::cerr << "Error sending message to Alice." << std::endl;  // Handle sending error
    }
    
    std::string outputBob;
    std::string inputBob;
    
    // Iterate over the number of scanpaths Alice sent
    for (int k = 0; k < getNumberOfScanpaths(aliceScanpathNumberMessage); ++k) {
        for (const std::string& tempScanpath : allScanpaths) {
            processScanpath(tempScanpath, outputBob, inputBob);  // Process each scanpath
        }
    }
    
    logStatistics();  // Log statistical information about the process
}

// Function definition for process
bool BobProcessor::process(
    std::string& resultString, // Reference to a string for storing the result
    const std::string& inputString // Input string received from Alice
){
    std::string responseToAlice; // Variable to store the response to be sent to Alice

    // Initialize various timers to measure the durations of different parts of the process
    Timer timerInitialization;
    Timer timerPublicPaillier;

    // Process the initial message received from Alice
    std::vector<std::vector<MPZ_Wrapper>> encryptedAliceDistanceMatrix = processInitialMessage(inputString);
    // Update statistics for the creation of the public Paillier cryptosystem
    statsCreationOfPublicPaillier.update(timerPublicPaillier.getDuration());
    // Determine the length of Alice's scanpath from the received matrix
    int aliceScanpathLength = encryptedAliceDistanceMatrix.size();

    // Start timer for encrypting Bob's scanpath
    Timer timerScanpathEncryption;
    // Encrypt Bob's scanpath
    std::vector<MPZ_Wrapper> encryptedBobScanpath = encryptScanpath(public_paillier, bobScanpath);
    // Update encryption statistics
    statsEncryptionOfScanpath.update(timerScanpathEncryption.getDuration());

    // Initialize cost variables for the computation
    initializeCosts();
    // Create a vector of indices corresponding to letters in Bob's scanpath
    std::vector<int> bobIndicesVector = createBobScanpathLetterIndicesVector();

    // Determine the length of Bob's scanpath
    int bobScanpathLength = encryptedBobScanpath.size();

    // Initialize the alignment matrix with the lengths of Alice's and Bob's scanpaths
    Matrix encryptedDistanceMatrix = initializeAlignmentMatrix(aliceScanpathLength, bobScanpathLength, encryptedGapValue);

    // Initialize a set of candidates for processing
    std::set<std::pair<int, int>, std::less<>> candidates;
    candidates.insert({1, 1});

    // Update initialization statistics
    statsInitialization.update(timerInitialization.getDuration());

    // Vector to store elapsed times for each step
    std::vector<std::chrono::duration<double>> elapsedTimesEachStep;

    // Counter for the number of candidates processed
    int candidateWhileCounter = 0;
    // Process each candidate until there are no more candidates
    while (!candidates.empty()) {
        oneStepProcess(candidates, encryptedDistanceMatrix, candidateWhileCounter, encryptedDeletionCost, encryptedInsertionCost, bobIndicesVector, encryptedAliceDistanceMatrix, aliceScanpathLength, bobScanpathLength);
    }
    // Finalize the candidate while counter
    candidateWhileCounterFinal = candidateWhileCounter;

    // Get the final result from the computation
    getFinalResult(encryptedDistanceMatrix, aliceScanpathLength, bobScanpathLength);
    // Output the final comparison result
    std::cout << finalComparisonResult << std::endl;

    // Return true to indicate successful completion of the process
    return true;
}



// This function sets up a connection to Alice.
void BobProcessor::initializeConnection() {
    std::cout << "Trying to connect to Alice on port: " << portNumber << std::endl;
    if (!socketComm.init_as_client("127.0.0.1", portNumber)) {
        std::cerr << "Error connecting to server." << std::endl;
        return;  // Exit or handle error appropriately
    }
}


// Function definition for processScanpath
void BobProcessor::processScanpath(
    const std::string& tempScanpath, // The scanpath for Bob to process
    std::string outputBob, // A string to store the output for Bob
    std::string inputBob // A string to receive input, presumably from Alice
){
    // Set the current scanpath to the one provided in the argument
    bobScanpath = tempScanpath;

    // Log the start of the comparison process
    logsComparisonStart();

    // Start a timer to measure the duration of the algorithm
    Timer timerStartAlgorithm;

    // Receive a message from Alice and store it in inputBob
    inputBob = socketComm.receive_message();

    // Process the received message along with the current scanpath
    process(outputBob, inputBob);

    // Log the end of the comparison process and the duration it took
    logsComparisonEnd(timerStartAlgorithm);

    // Increment the iteration counter
    iterationCounter = iterationCounter + 1;
}



// This function logs the accumulated statistics of the processing.
void BobProcessor::logStatistics() {

    log("Scanpath comparison is done.",
        buildLogMessage(" Final statistics for individual scanpath comparison: ", statsInidividualScanpathComparison),
        logFilePath, "Bob", iterationCounter, "Value");
    
    log("Scanpath comparison is done.",
        buildLogMessage("Final statistics for initialization: ", statsInitialization),
        logFilePath, "Bob", iterationCounter, "Value");
    
    log("Scanpath comparison is done.",
        buildLogMessage("Final statistics for encryption of scanpath: ", statsEncryptionOfScanpath),
        logFilePath, "Bob", iterationCounter, "Value");
    
    log("Scanpath comparison is done.",
        buildLogMessage("Final statistics for each cell calculcation(each small step): ", statsEachStep),
        logFilePath, "Bob", iterationCounter, "Value");
    
    log("Scanpath comparison is done.",
        buildLogMessage("Final statistics for creation of public paillier object: ", statsCreationOfPublicPaillier),
        logFilePath, "Bob", iterationCounter, "Value");
    
    log("Scanpath comparison is done.",
        buildLogMessage("Final statistics for minimum communication with Alice: ", statsMinCommunication),
        logFilePath, "Bob", iterationCounter, "Value");
    
    // Optional: If you have any final logs or clean up, you can perform it here.
    flush_remaining_logs(logFilePath); // Assuming you have a mechanism to flush pending logs.
    
    double averageTime = static_cast<double>(statsInidividualScanpathComparison.getMean());
    std::cout << "Average time for each iteration: " << averageTime << " microseconds." << std::endl;
}


void BobProcessor::initializeCosts(){
    mpz_t deletionCost, insertionCost, gapValue; // Declare mpz_t variables for costs
    mpz_inits(encryptedDeletionCost,encryptedInsertionCost, encryptedGapValue,NULL);

    mpz_inits(deletionCost, insertionCost, gapValue, NULL); // Initialize the mpz_t variables

    // Set the deletion, insertion, and gap costs to 1
    // "10" here refers to the base (decimal) for the number representation
    mpz_set_str(deletionCost, "1", 10); // Set deletion cost to 1
    mpz_set_str(insertionCost, "1", 10); // Set insertion cost to 1
    mpz_set_str(gapValue, "1", 10); // Set gap value to 1

    // Encrypt the gap value using the Paillier cryptosystem
    public_paillier.encrypt(encryptedGapValue, gapValue);

    // Encrypt the deletion cost using the Paillier cryptosystem
    public_paillier.encrypt(encryptedDeletionCost, deletionCost);

    // Encrypt the insertion cost using the Paillier cryptosystem
    public_paillier.encrypt(encryptedInsertionCost, insertionCost);

    // Clear the mpz_t variables to free up resources
    mpz_clears(deletionCost, insertionCost, gapValue, NULL);
}



std::vector<std::vector<MPZ_Wrapper>> BobProcessor::processInitialMessage(
    const std::string& inputString // Input string that contains the initial message
){
    // Convert the input string to a vector of MPZ_Wrapper objects
    // MPZ_Wrapper is a wrapper around the mpz_t type from GMP for handling big numbers
    std::vector<MPZ_Wrapper> inputValues = convertStringToMPZ_Wrapper(inputString);


    // Check if the public parameters n and g need to be updated
    // These parameters are part of the Paillier cryptosystem's public key
    if (mpz_cmp(public_n, inputValues[1].value) != 0 || mpz_cmp(public_g, inputValues[2].value) != 0) {
        // Update public parameters n and g
        mpz_inits(public_n,public_g,NULL);
        mpz_set(public_n, inputValues[1].value);
        mpz_set(public_g, inputValues[2].value);
        public_paillier.updatePublicParams(public_n, public_g);
    }

    // Extract the encrypted Alice distance vector from the input values
    std::vector<MPZ_Wrapper> encryptedAliceDistanceVector;
    for (size_t i = 3; i < inputValues.size(); ++i) {
        MPZ_Wrapper value;
        mpz_init_set(value.value, inputValues[i].value);
        encryptedAliceDistanceVector.push_back(value);
    }

    // Convert the flat vector of encrypted distances into a matrix
    // The size of the matrix is determined by the size of the alphabet used in the process
    std::vector<std::vector<MPZ_Wrapper>> encryptedAliceDistanceMatrix = unflattenVector(encryptedAliceDistanceVector, alphabet.size());

    // Clear and shrink the input vectors to free up memory
    inputValues.clear();
    inputValues.shrink_to_fit();
    encryptedAliceDistanceVector.clear();
    encryptedAliceDistanceVector.shrink_to_fit();

    // Return the matrix of encrypted Alice distances
    return encryptedAliceDistanceMatrix;
}



std::vector<int> BobProcessor::createBobScanpathLetterIndicesVector(){
    std::vector<int> bobIndicesVector;
    for (char c : bobScanpath) {
        auto it = alphabet.find(c);
        if (it != alphabet.end()) {
            bobIndicesVector.push_back(it->second);
        } else {
            // Handle case where the character is not found in the alphabet, if needed
            bobIndicesVector.push_back(-1); // For example, you can push -1 to indicate it's not found
        }
    }
    return bobIndicesVector;
}

Matrix BobProcessor::initializeAlignmentMatrix(
    int aliceScanpathLength, // Length of Alice's scanpath
    int bobScanpathLength, // Length of Bob's scanpath
    mpz_t encryptedGapValue // Encrypted value of the gap, used in the matrix initialization
){
    // Create a matrix with dimensions (aliceScanpathLength + 1) x (bobScanpathLength + 1)
    Matrix encryptedDistanceMatrix(aliceScanpathLength + 1, bobScanpathLength + 1);

    mpz_t tempValue; // Temporary variable for initialization values
    mpz_init(tempValue); // Initialize the temporary variable

    // Initialize all elements of the matrix to UNCOMPUTED_VALUE
    for (int i = 0; i < encryptedDistanceMatrix.getRows(); ++i) {
        for (int j = 0; j < encryptedDistanceMatrix.getCols(); ++j) {
            mpz_set_si(tempValue, UNCOMPUTED_VALUE); // Set tempValue to UNCOMPUTED_VALUE
            encryptedDistanceMatrix.set(i, j, tempValue); // Set the matrix cell to tempValue
        }
    }

    // Set tempValue to 0 and encrypt it
    mpz_set_str(tempValue, "0", 10); // Set tempValue to "0"
    public_paillier.encrypt(tempValue, tempValue); // Encrypt tempValue

    mpz_t tempResult; // Temporary variable for calculations
    mpz_init(tempResult); // Initialize the temporary variable

    // Initialize the first column of the matrix
    int j = 0;
    for (int i = 0; i < aliceScanpathLength + 1; ++i) {
        encryptedDistanceMatrix.set(i, j, tempValue); // Set the value in the matrix
        public_paillier.add(tempResult, tempValue, encryptedGapValue); // Add encryptedGapValue to tempValue
        mpz_set(tempValue, tempResult); // Update tempValue for the next iteration
    }

    // Reset tempValue to 0 and encrypt it again
    mpz_set_str(tempValue, "0", 10);
    public_paillier.encrypt(tempValue, tempValue);

    // Initialize the first row of the matrix
    int i = 0;
    for (int j = 0; j < bobScanpathLength + 1; ++j) {
        encryptedDistanceMatrix.set(i, j, tempValue); // Set the value in the matrix
        public_paillier.add(tempResult, tempValue, encryptedGapValue); // Add encryptedGapValue to tempValue
        mpz_set(tempValue, tempResult); // Update tempValue for the next iteration
    }

    mpz_clears(tempValue, tempResult, NULL); // Clear the temporary mpz_t variables

    return encryptedDistanceMatrix; // Return the initialized matrix
}




void BobProcessor::updateCandidates(
    int i, int j, // Current indices in the matrix
    int aliceScanpathLength, // The length of Alice's scanpath
    int bobScanpathLength, // The length of Bob's scanpath
    Matrix encryptedDistanceMatrix, // The matrix containing encrypted distances
    std::set<std::pair<int, int>, std::less<>>& candidates // Set of candidate index pairs
){
    mpz_t tempValue; // Temporary variable to store matrix values
    mpz_init(tempValue); // Initialize the temporary variable

    // Define potential next positions to move in the matrix
    std::vector<std::pair<int, int>> potential_next = {{i+1, j+1}, {i+1, j}, {i, j+1}};

    // Iterate through each potential next position
    for (auto& p : potential_next) {
        int pi = p.first, pj = p.second; // Extract the indices of the potential next position
        
        // Check if the potential position is within the bounds of the matrix
        if (pi <= aliceScanpathLength && pj <= bobScanpathLength) {
            // Check if the left position is computed
            encryptedDistanceMatrix.get(pi-1, pj, tempValue);
            bool leftComputed = mpz_cmp_si(tempValue, UNCOMPUTED_VALUE) != 0;
            
            // Check if the top position is computed
            encryptedDistanceMatrix.get(pi, pj-1, tempValue);
            bool topComputed = mpz_cmp_si(tempValue, UNCOMPUTED_VALUE) != 0;
            
            // Check if the diagonal position is computed
            encryptedDistanceMatrix.get(pi-1, pj-1, tempValue);
            bool diagComputed = mpz_cmp_si(tempValue, UNCOMPUTED_VALUE) != 0;
            
            // If all three positions (left, top, diagonal) are computed, add this position to the candidates
            if (leftComputed && topComputed && diagComputed) {
                candidates.insert(p);
            }
        }
    }
    mpz_clear(tempValue); // Clear the temporary mpz_t variable
}




void BobProcessor::inverseAffineTransformation(
    mpz_t& result, // Variable to store the result of the inverse transformation
    mpz_t& value, // The value to be transformed back
    mpz_t& delta1Inverse, // Inverse of delta1 used in the original affine transformation
    mpz_t& delta2Inverse, // Inverse of delta2 used in the original affine transformation
    mpz_t& rho2_inverse // Inverse of rho2 used in the original affine transformation
){
    mpz_t temp; // Temporary variable for intermediate calculations
    mpz_init(temp); // Initialize the temporary variable

    // First part of reversing the affine transformation:
    // temp = value - delta2 (addition with additive inverse in the context of the Paillier cryptosystem)
    public_paillier.scalarAdd(temp, value, delta2Inverse);

    // Second part of the transformation:
    // temp = temp / rho2 (multiplication with multiplicative inverse in the context of the Paillier cryptosystem)
    public_paillier.scalarMultiply(temp, temp, rho2_inverse);

    // Final part of the transformation:
    // result = value - delta1 (addition with additive inverse in the context of the Paillier cryptosystem)
    public_paillier.scalarAdd(result, temp, delta1Inverse);

    // Clear the temporary mpz_t variable to free up resources
    mpz_clear(temp);
}


void BobProcessor::inverseOrderPreservingMasking(
    mpz_t& xMinValue, // The value to which the inverse masking will be applied
    mpz_t& encryptedSum, // The sum of encrypted values used in the masking process
    mpz_t& rhoInverse, // The inverse of the rho value used in the masking process
    bool& randomDefinitionControl // A boolean flag indicating which masking method was used
){
    // Check if the random definition control flag is true
    if (randomDefinitionControl) {
        // If true, it indicates that the order preserving masking method was used

        // Add the encrypted sum to the xMinValue
        public_paillier.add(xMinValue, xMinValue, encryptedSum);

        // Multiply the result by the inverse of rho (unmasking the value)
        public_paillier.scalarMultiply(xMinValue, xMinValue, rhoInverse);
    } else {
        // If the flag is false, a simpler masking method was used
        // Directly multiply xMinValue by the inverse of rho to unmask the value
        public_paillier.scalarMultiply(xMinValue, xMinValue, rhoInverse);
    }
}



// Function definition for BobProcessor::getMinimumFromAlice
void BobProcessor::getMinimumFromAlice(
    mpz_t& xMinValue, // Reference to an mpz_t variable to store the minimum value
    const std::vector<MPZ_Wrapper> encryptedMinCalculationVector // Vector of encrypted values for finding the minimum
){
    Timer timerMinCommunication; // Start a timer to measure communication duration

    // Send the vector of encrypted values to Alice
    sendEncryptedResult(socketComm, encryptedMinCalculationVector, 1);

    // Receive the response from Alice contains the minimum value
    std::string outputAlice = socketComm.receive_message();

    // Update the communication statistics with the duration of the communication
    statsMinCommunication.update(timerMinCommunication.getDuration());

    // Convert the received message (string) back to a vector of MPZ_Wrapper objects
    std::vector<MPZ_Wrapper> tempOutput = convertStringToMPZ_Wrapper(outputAlice);

    // Set the second element in the received vector to the minimum value
    mpz_set(xMinValue, tempOutput[1].value);

    // Clear the temporary vector and release its memory
    tempOutput.clear();
    tempOutput.shrink_to_fit();
}




void BobProcessor::applyAfineTransform( mpz_t& delta1Inverse, mpz_t& delta2Inverse,mpz_t& rho_2_inverse, mpz_t& x1Prime,  mpz_t& x2Prime,  mpz_t& x3Prime, std::vector<MPZ_Wrapper>& encryptedMinCalculationVector){
    
    mpz_t rho_2,delta_1,delta_2;
    mpz_inits(rho_2,delta_1,delta_2,NULL);
    
    

    // Generate random values for rho_2, delta_1, and delta_2
    public_paillier.randomWithInverse(rho_2,rho_2_inverse);
    public_paillier.randomWithAdditiveInverse(delta_1,delta1Inverse);
    public_paillier.randomWithAdditiveInverse(delta_2,delta2Inverse);
    
    // Now, apply the affine transformation
    
    // x1'' = rho_2 * (x1' + delta_1) + delta_2
    public_paillier.scalarAdd(x1Prime, x1Prime, delta_1);
    public_paillier.scalarMultiply(x1Prime, x1Prime, rho_2);
    public_paillier.scalarAdd(x1Prime, x1Prime, delta_2);
    
    // x2'' = rho_2 * (x2' + delta_1) + delta_2
    public_paillier.scalarAdd(x2Prime, x2Prime, delta_1);
    public_paillier.scalarMultiply(x2Prime, x2Prime, rho_2);
    public_paillier.scalarAdd(x2Prime, x2Prime, delta_2);
    
    // x3'' = rho_2 * (x3' + delta_1) + delta_2
    public_paillier.scalarAdd(x3Prime, x3Prime, delta_1);
    public_paillier.scalarMultiply(x3Prime, x3Prime, rho_2);
    public_paillier.scalarAdd(x3Prime, x3Prime, delta_2);
    
    // Store the results in a vector
    encryptedMinCalculationVector.push_back(MPZ_Wrapper(x1Prime));
    encryptedMinCalculationVector.push_back(MPZ_Wrapper(x2Prime));
    encryptedMinCalculationVector.push_back(MPZ_Wrapper(x3Prime));
    mpz_clears(rho_2,delta_1,delta_2,NULL);

}



void BobProcessor::applyOrderPreservingMasking(
    const mpz_t x1, const mpz_t x2, const mpz_t x3, // Input values (mpz_t types) to be masked
    mpz_t& encryptedSum,  // Reference to a variable to store the sum of the original values
    mpz_t& rhoInverse, // Reference to a variable to store the inverse of rho
    mpz_t& x1Prime, mpz_t& x2Prime, mpz_t& x3Prime, // Output variables for the masked values
    bool& randomDefinitionControl // Boolean to control random decision in masking process
){
    // Temporary mpz_t variables for storing the inverses and a random value rho
    mpz_t x1_inverse, x2_inverse, x3_inverse, rho;
    mpz_inits(x1_inverse, x2_inverse, x3_inverse, rho, NULL);

    // Generate a random number to decide the masking technique
    double randomDefinition = randomDoubleFromUrandom();
    
    // Calculate the multiplicative inverses of x1, x2, and x3 under modulo public_paillier.n_square
    mpz_invert(x1_inverse, x1, public_paillier.n_square);
    mpz_invert(x2_inverse, x2, public_paillier.n_square);
    mpz_invert(x3_inverse, x3, public_paillier.n_square);
    
    if (randomDefinition >= 0.5) {
        // Apply a specific masking technique
        randomDefinitionControl = true;
        // Generate rho and its shifted inverse
        public_paillier.randomWithShiftedInverse(rho, rhoInverse);
        // Calculate the encrypted sum of x1, x2, and x3
        public_paillier.add(encryptedSum, x1, x2);
        public_paillier.add(encryptedSum, encryptedSum, x3);
        
        // Apply the masking operation to x1, x2, and x3
        // x1Prime = (rho * x1) - x2 - x3
        public_paillier.scalarMultiply(x1Prime, x1, rho);
        public_paillier.add(x1Prime, x1Prime, x2_inverse);
        public_paillier.add(x1Prime, x1Prime, x3_inverse);
        
        // x2Prime = (rho * x2) - x1 - x3
        public_paillier.scalarMultiply(x2Prime, x2, rho);
        public_paillier.add(x2Prime, x2Prime, x1_inverse);
        public_paillier.add(x2Prime, x2Prime, x3_inverse);
        
        // x3Prime = (rho * x3) - x1 - x2
        public_paillier.scalarMultiply(x3Prime, x3, rho);
        public_paillier.add(x3Prime, x3Prime, x1_inverse);
        public_paillier.add(x3Prime, x3Prime, x2_inverse);
        
    } else {
        // Apply a basic masking technique
        public_paillier.randomWithInverseMultiplication(rho, rhoInverse);
        
        // Multiply each x1, x2, and x3 by rho
        public_paillier.scalarMultiply(x1Prime, x1, rho);
        public_paillier.scalarMultiply(x2Prime, x2, rho);
        public_paillier.scalarMultiply(x3Prime, x3, rho);
    }
    
    // Clear the temporary mpz_t variables
    mpz_clears(x1_inverse, x2_inverse, x3_inverse, rho, NULL);
}


void BobProcessor::calculateSequenceEditingCosts(
    int i, int j,
    mpz_t& x1, mpz_t& x2, mpz_t& x3,
    const mpz_t encryptedDeletionCost,
    const mpz_t encryptedInsertionCost,
    const Matrix encryptedDistanceMatrix,
    const std::vector<int> bobIndicesVector,
    const std::vector<std::vector<MPZ_Wrapper>> encryptedAliceDistanceMatrix
){
    // Declare variables to hold temporary values for matrix positions and updates
    mpz_t tempLeftVal,tempUpperVal,tempDiagonalVal,tempUpdatedLeftVal,tempUpdatedUpperVal,tempUpdatedDiagonalVal;
    
    // Initialize these temporary variables to be used for arithmetic operations
    mpz_inits(tempLeftVal,tempUpperVal,tempDiagonalVal,tempUpdatedLeftVal,tempUpdatedUpperVal,tempUpdatedDiagonalVal,NULL);

    // Retrieve values from the encrypted distance matrix for the current (i,j) values
    // These values represent costs calculated in previous steps of the algorithm
    encryptedDistanceMatrix.get(i-1,j,tempLeftVal);
    encryptedDistanceMatrix.get(i,j-1,tempUpperVal);
    encryptedDistanceMatrix.get(i-1,j-1,tempDiagonalVal);

    // Add the encrypted deletion cost to the left value
    public_paillier.add(tempUpdatedLeftVal,tempLeftVal,encryptedDeletionCost);
    // Add the encrypted insertion cost to the upper value
    public_paillier.add(tempUpdatedUpperVal,tempUpperVal,encryptedInsertionCost);
    // Add the value from Alice's encrypted distance matrix to the diagonal value
    // This is based on the current indices i and j, adjusted for zero-based indexing
    // This is substitution cost.
    public_paillier.add(tempUpdatedDiagonalVal,tempDiagonalVal,encryptedAliceDistanceMatrix[i-1][bobIndicesVector[j-1]].value);

    // Set the values of x1, x2, and x3 to the respective updated values
    mpz_set(x1, tempUpdatedDiagonalVal);
    mpz_set(x2, tempUpdatedLeftVal);
    mpz_set(x3, tempUpdatedUpperVal);

    // Clear the temporary variables to free up resources
    mpz_clears(tempLeftVal,tempUpperVal,tempDiagonalVal,tempUpdatedLeftVal,tempUpdatedUpperVal,tempUpdatedDiagonalVal,NULL);
}



void BobProcessor::calculateCurrentCell(int i, int j, const mpz_t encryptedDeletionCost, const mpz_t encryptedInsertionCost, Matrix& encryptedDistanceMatrix, const std::vector<int> bobIndicesVector,const std::vector<std::vector<MPZ_Wrapper>> encryptedAliceDistanceMatrix){
    mpz_t x1,x2,x3,encryptedSum, x1Prime,x2Prime,x3Prime,xMinValue,rhoInverse, rho_2_inverse,delta1Inverse,delta2Inverse;
    mpz_inits(x1, x2, x3, encryptedSum, x1Prime, x2Prime, x3Prime,   xMinValue, rhoInverse,rho_2_inverse,delta1Inverse,delta2Inverse,NULL);
    calculateSequenceEditingCosts(i, j,  x1, x2, x3,  encryptedDeletionCost, encryptedInsertionCost,  encryptedDistanceMatrix,  bobIndicesVector, encryptedAliceDistanceMatrix);
    bool randomDefinitionControl = false;
    std::vector<MPZ_Wrapper> encryptedMinCalculationVector;
    encryptedMinCalculationVector.reserve(3);
    applyOrderPreservingMasking(x1,  x2, x3, encryptedSum, rhoInverse,  x1Prime,  x2Prime, x3Prime, randomDefinitionControl);
    applyAfineTransform( delta1Inverse, delta2Inverse, rho_2_inverse, x1Prime,  x2Prime,  x3Prime, encryptedMinCalculationVector);
    randomlyPermutateVector(encryptedMinCalculationVector);
    getMinimumFromAlice(xMinValue,encryptedMinCalculationVector);
    // Reverse the affine transformation for x1
    inverseAffineTransformation(xMinValue, xMinValue,   delta1Inverse,  delta2Inverse,  rho_2_inverse);
    inverseOrderPreservingMasking(xMinValue,  encryptedSum,  rhoInverse,  randomDefinitionControl);
    encryptedMinCalculationVector.clear();
    encryptedMinCalculationVector.shrink_to_fit();
    encryptedDistanceMatrix.set(i, j,xMinValue);
    mpz_clears(x1, x2, x3, encryptedSum, x1Prime, x2Prime, x3Prime,   xMinValue, rhoInverse,rho_2_inverse,delta1Inverse,delta2Inverse,NULL);

}


void BobProcessor::selectRandomCandidate(
    int& i, int& j, // References to integers to store the selected candidate's indices
    std::set<std::pair<int, int>, std::less<>>& candidates // Reference to a set of candidate pairs
){
    // Generate a random index based on the size of the candidates set
    size_t offset = random_index_from_urandom(candidates.size());

    // Get an iterator to the randomly selected candidate
    // std::next advances the iterator from the beginning of the set by 'offset' positions
    auto it = std::next(candidates.begin(), offset);
    
    // Assign the selected candidate's indices to i and j
    i = it->first;
    j = it->second;

    // Remove the selected candidate from the set
    candidates.erase(it);
}


// Function definition for BobProcessor::oneStepProcess
void BobProcessor::oneStepProcess(
    std::set<std::pair<int, int>, std::less<>>& candidates, // Set of candidate pairs for processing
    Matrix& encryptedDistanceMatrix, // Reference to the matrix storing encrypted distances
    int candidateWhileCounter, // Counter for the number of iterations in the while loop
    const mpz_t encryptedDeletionCost, // Encrypted cost of a deletion operation
    const mpz_t encryptedInsertionCost, // Encrypted cost of an insertion operation
    const std::vector<int> bobIndicesVector, // Vector of indices corresponding to Bob's scanpath
    const std::vector<std::vector<MPZ_Wrapper>> encryptedAliceDistanceMatrix, // Matrix of Alice's encrypted distances
    const int aliceScanpathLength, // Length of Alice's scanpath
    const int bobScanpathLength // Length of Bob's scanpath
){
    int i,j; // Declare variables to store the indices of the selected candidate

    candidateWhileCounter=candidateWhileCounter+1; // Increment the counter for each step

    Timer timerStartEachStep; // Start a timer to measure the duration of this step

    // Select a random candidate from the set of candidates
    // i and j will be set to the indices of this candidate
    selectRandomCandidate(i,j , candidates);

    // Calculate the current cell based on the selected candidate
    // This involves computing the costs for various operations (insertion, deletion, etc.)
    calculateCurrentCell( i,  j,   encryptedDeletionCost,   encryptedInsertionCost,  encryptedDistanceMatrix,   bobIndicesVector, encryptedAliceDistanceMatrix);

    // Update the set of candidates based on the current cell's computation
    // This involves adding new candidates
    updateCandidates( i, j,aliceScanpathLength,bobScanpathLength, encryptedDistanceMatrix, candidates);

    // Update the statistics for each step with the duration of the current step
    statsEachStep.update(timerStartEachStep.getDuration());
}


void BobProcessor::getFinalResult(const Matrix encryptedDistanceMatrix, const int aliceScanpathLength,const int bobScanpathLength){

    mpz_t tempResult;  // Declare a variable of type mpz_t to store temporary results
    mpz_init(tempResult);  // Initialize the mpz_t variable to ensure it's ready for use

    // Retrieve the value from the encrypted distance matrix at the last position(m,n)
    // by the lengths of the Alice and Bob scanpaths
    encryptedDistanceMatrix.get(aliceScanpathLength, bobScanpathLength,tempResult);

    std::vector<MPZ_Wrapper> finalResultValue;  // Create a vector to store the final result value
    // Wrap the temporary result in an MPZ_Wrapper and add it to the finalResultValue vector
    finalResultValue.push_back(MPZ_Wrapper(tempResult));

    // Send the final result to Alice using the encrypted communication channel
    // '2' is a flag for the type of message being sent
    sendEncryptedResult(socketComm, finalResultValue, 2);

    // Receive a message from Alice and store it in outputAlice
    std::string outputAlice = socketComm.receive_message();

    // Store Alice's response in the class member variable finalComparisonResult for later use
    finalComparisonResult = outputAlice;
    
    mpz_clear(tempResult);  // Clear the temporary mpz_t variable to free up resources
}


// Class method for starting the comparison of logs in BobProcessor
void BobProcessor::logsComparisonStart(){
    // Resetting statistics for individual scanpath comparison
    statsInidividualScanpathComparison.reset();
    // Resetting statistics for initialization process
    statsInitialization.reset();
    // Resetting statistics for encryption of the scanpath
    statsEncryptionOfScanpath.reset();
    // Resetting statistics for each step in the process
    statsEachStep.reset();
    // Resetting statistics for the creation of Public Paillier key
    statsCreationOfPublicPaillier.reset();
    // Resetting statistics for minimum communication required
    statsMinCommunication.reset();

    // Logging various events and data at the start of the loop
    log("Loop", "Loop started.", logFilePath, "Bob", iterationCounter,"Start");
    log("Scanpath", bobScanpath, logFilePath, "Bob", iterationCounter,"Value");
    log("Scanpath Length",  std::to_string(bobScanpath.length()), logFilePath, "Bob", iterationCounter,"Value");
    log("Algorithm", "Algorithm started.", logFilePath, "Bob", iterationCounter,"Start");
    
    // Starting a timer to measure the duration of the algorithm execution
    Timer timerStartAlgorithm ;
}

// Class method for ending the comparison of logs in BobProcessor
void BobProcessor::logsComparisonEnd(Timer timerStartAlgorithm){
    // Calculating the duration of the algorithm execution
    std::chrono::duration<double> iterationTime = timerStartAlgorithm.getDuration();
    
    // Updating the statistics with the duration of this iteration
    statsInidividualScanpathComparison.update(iterationTime);
    
    // Logging the end of the algorithm with detailed information
    log("Algorithm", "Algorithm ended. Final Decrypted Result: " + finalComparisonResult + ". Total while iteration: "+ std::to_string(candidateWhileCounterFinal) + " It took " + durationToString(iterationTime,5) + " microseconds.", logFilePath, "Bob", iterationCounter,"End");
    
    // Logging all the statistics collected during the process
    logStatistics();
}
