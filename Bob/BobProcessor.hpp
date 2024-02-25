#ifndef BOB_PROCESSOR_HPP
#define BOB_PROCESSOR_HPP

#include <string>
#include <vector>
#include <set>
#include "Paillier.hpp"
#include "utils.hpp"
#include "SocketComm.hpp"

// Include other necessary headers

class BobProcessor {
private:
    SocketComm socketComm;  // Private member for SocketComm
    std::string bobScanpath;
    mpz_t public_n;
    mpz_t public_g;
    mpz_t init_public_n;
    mpz_t init_public_g;
    mpz_t encryptedDeletionCost;
    mpz_t encryptedInsertionCost;
    mpz_t encryptedGapValue;
    PublicPaillier public_paillier;
    std::string rootFolder;
    std::string filename;
    std::string logFilePath;
    int portNumber;
    const int UNCOMPUTED_VALUE = -9999;
    std::vector<std::string> allScanpaths;
    std::map<char, int> alphabet;
    std::vector<int>createBobScanpathLetterIndicesVector();
    Matrix initializeAlignmentMatrix(int aliceScanpathLength,int bobScanpathLength,mpz_t encryptedGapValue );

    void updateCandidates(int i, int j,
                            int aliceScanpathLength,
                            int bobScanpathLength,
                            Matrix encryptedDistanceMatrix,
                          std::set<std::pair<int, int>, std::less<>>& candidates);
    
    void inverseAffineTransformation(mpz_t& result, mpz_t& value, mpz_t& delta1_inverse, mpz_t& delta2_inverse, mpz_t& rho2_inverse);

    
    void inverseOrderPreservingMasking(mpz_t& x_min_value, mpz_t& encryptedSum, mpz_t& rho_inverse, bool& randomDefinitionControl);


    PublicPaillier initializePublicPaillier() {
        // Call mpz initialization functions here
        mpz_inits(init_public_g,init_public_n,NULL);
        mpz_init_set_ui(init_public_n, 10);
        mpz_init_set_ui(init_public_g, 20);
        return PublicPaillier(init_public_n, init_public_g);
    }
    void iterateOverScanpaths(std::string outputBob, std::string inputBob);
    void initializeConnection();
    void processScanpath(const std::string& tempScanpath, std::string outputBob, std::string inputBob);

    void initializeCosts();
    void getFinalResult(const Matrix encryptedDistanceMatrix, const int aliceScanpathLength,const int bobScanpathLength);

    void getMinimumFromAlice(mpz_t& x_min_value,const std::vector<MPZ_Wrapper> encryptedMinCalculationVector);
    void applyAfineTransform( mpz_t& delta_1_inverse, mpz_t& delta_2_inverse,mpz_t& rho_2_inverse, mpz_t& x1_prime,  mpz_t& x2_prime,  mpz_t& x3_prime, std::vector<MPZ_Wrapper>& encryptedMinCalculationVector);

    void applyOrderPreservingMasking(const mpz_t x1, const mpz_t x2,const mpz_t x3, mpz_t& encryptedSum,  mpz_t& rho_inverse, mpz_t& x1_prime, mpz_t& x2_prime, mpz_t& x3_prime, bool& randomDefinitionControl);
    void calculateSequenceEditingCosts(int i, int j, mpz_t& x1, mpz_t& x2, mpz_t& x3, const mpz_t encryptedDeletionCost, const mpz_t encryptedInsertionCost, const Matrix encryptedDistanceMatrix, const std::vector<int> bobIndicesVector,const std::vector<std::vector<MPZ_Wrapper>> encryptedAliceDistanceMatrix);
    void calculateCurrentCell(int i, int j, const mpz_t encryptedDeletionCost, const mpz_t encryptedInsertionCost, Matrix& encryptedDistanceMatrix, const std::vector<int> bobIndicesVector,const std::vector<std::vector<MPZ_Wrapper>> encryptedAliceDistanceMatrix);
    void selectRandomCandidate(int& i, int& j ,std::set<std::pair<int, int>, std::less<>>& candidates);
    void oneStepProcess(std::set<std::pair<int, int>, std::less<>>& candidates,Matrix& encryptedDistanceMatrix, int candidateWhileCounter, const mpz_t encryptedDeletionCost, const mpz_t encryptedInsertionCost, const std::vector<int> bobIndicesVector, const std::vector<std::vector<MPZ_Wrapper>> encryptedAliceDistanceMatrix,const int aliceScanpathLength, const int bobScanpathLength );

    std::vector<std::vector<MPZ_Wrapper>> processInitialMessage(const std::string& inputString);
    void logsComparisonStart();
    void logsComparisonEnd(Timer timerStartAlgorithm);


public:
    
    // Adjusted constructor to take a SocketComm object
    BobProcessor(const SocketComm& socketCommObj, const std::string& rootFolder,
                 const std::string& filenamePart, const std::string& logFilenamePart,
                 int inputPortNumber, int fileIndex)
        : socketComm(socketCommObj), rootFolder(rootFolder),
    public_paillier(initializePublicPaillier())
    {
        mpz_inits(encryptedDeletionCost,encryptedInsertionCost, encryptedGapValue,NULL);

        portNumber = inputPortNumber;
        filename = rootFolder + "Datasets/" + filenamePart + "_" + std::to_string(fileIndex) + ".txt";
        logFilePath = rootFolder + "Output/" + logFilenamePart;

    }

    int iterationCounter;
    OnlineStats statsInidividualScanpathComparison;
    OnlineStats statsInitialization;
    OnlineStats statsEncryptionOfScanpath;
    OnlineStats statsEachStep;
    OnlineStats statsCreationOfPublicPaillier;
    OnlineStats statsMinCommunication;
    int candidateWhileCounterFinal;
    std::string finalComparisonResult;

    void startProcess();
    bool process(std::string& resultString, const std::string& inputString);
    void logStatistics();
    
};

#endif // BOB_PROCESSOR_HPP
