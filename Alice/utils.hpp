#ifndef UTILS_HPP
#define UTILS_HPP

#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <gmp.h>
#include <stdexcept>
#include <chrono>
#include <random>
#include <algorithm>
#include <sstream>
#include "Paillier.hpp"
#include "SocketComm.hpp"


// Declaration of MPZ_Wrapper class, which wraps GMP's mpz_t
class MPZ_Wrapper {
public:
    mpz_t value; // Wrapped GMP mpz_t variable

    // Default constructor
    MPZ_Wrapper();

    // Copy constructor from mpz_t
    MPZ_Wrapper(const mpz_t& other);

    // Copy constructor from another MPZ_Wrapper object
    MPZ_Wrapper(const MPZ_Wrapper& other);

    // Move constructor
    MPZ_Wrapper(MPZ_Wrapper&& other) noexcept;

    // Destructor
    ~MPZ_Wrapper();

    // Copy assignment operator from mpz_t
    MPZ_Wrapper& operator=(const mpz_t& other);

    // Copy assignment operator from another MPZ_Wrapper object
    MPZ_Wrapper& operator=(const MPZ_Wrapper& other);

    // Move assignment operator
    MPZ_Wrapper& operator=(MPZ_Wrapper&& other) noexcept;
};

// Definition of the Matrix class to store data in a matrix format
class Matrix {
private:
    unsigned int rows; // Number of rows in the matrix
    unsigned int cols; // Number of columns in the matrix
    std::vector<std::vector<MPZ_Wrapper>> data; // 2D vector to store MPZ_Wrapper values

public:
    // Constructor that initializes the Matrix with the specified number of rows and columns
    Matrix(unsigned int r, unsigned int c);

    // Copy constructor for the Matrix class
    Matrix(const Matrix& other);

    // Copy assignment operator for the Matrix class
    Matrix& operator=(const Matrix& other);

    // Function to set the value of an element in the matrix
    void set(int row, int col, const mpz_t& value);

    // Function to get the value of an element in the matrix
    void get(int row, int col, mpz_t& out_value) const;

    // Function to get the number of rows in the matrix
    int getRows() const;

    // Function to get the number of columns in the matrix
    int getCols() const;
};


// Timer class for measuring time durations
class Timer {
public:
    // Default constructor
    Timer();

    // Destructor
    ~Timer();

    // Function to get the duration since the timer started
    std::chrono::duration<double> getDuration() const;

private:
    std::chrono::time_point<std::chrono::high_resolution_clock> start; // Start time point
    mutable std::chrono::duration<double> duration; // Cached duration
};


// OnlineStats class for tracking statistics online
class OnlineStats {
public:
    // Default constructor
    OnlineStats();

    // Function to update statistics with a new duration value
    void update(std::chrono::duration<double> x);

    // Function to get the count of data points
    int getCount() const;

    // Function to get the mean value of the data points
    double getMean() const;

    // Function to get the variance of the data points
    double getVariance() const;

    // Function to get the standard deviation of the data points
    double getStandardDeviation() const;

    // Function to get the minimum value among the data points
    double getMin() const;

    // Function to get the maximum value among the data points
    double getMax() const;

    // Function to get the total sum of the data points
    double getTotalSum() const;

    // Function to reset the statistics
    void reset();

private:
    int n; // Count of data points
    double mean; // Mean of data points
    double M2; // Intermediate value for variance calculation
    double min; // Minimum value
    double max; // Maximum value
    double totalSum; // Sum of data points
};


// Function to create an alphabet mapping
std::map<char, int> createAlphabetMapping();
// Function to map a string to a vector of integers using an alphabet mapping
std::vector<int> mapStringUsingAlphabet(const std::map<char, int>& alphabet, const std::string& input);


// Function to randomly permute a vector of MPZ_Wrapper
void randomlyPermutateVector(std::vector<MPZ_Wrapper>& vec);
// Function to append two vectors of MPZ_Wrapper
std::vector<MPZ_Wrapper> appendVectors(const std::vector<MPZ_Wrapper>& vec1, const std::vector<MPZ_Wrapper>& vec2);

// Function to convert a string to a vector of MPZ_Wrapper
std::vector<MPZ_Wrapper> convertStringToMPZ_Wrapper(const std::string& str);
// Function to convert a vector of encrypted MPZ_Wrapper to a string representation
std::string encryptedResultToString(const std::vector<MPZ_Wrapper>& encryptedMappedResult);
// Function to add an operation flag to a string
void addOperationFlag(int number, std::string& str);
// Function to get an operation flag from a string
int getOperationFlag(const std::string input);

// Function to compute the difference matrix between a mapped scanpath and an alphabet mapping
std::vector<std::vector<int>> computeAlphabetDifferenceMatrix(const std::vector<int>& mappedScanpath, const std::map<char, int>& alphabet);
// Function to flatten a matrix into a vector
std::vector<int> flattenMatrix(const std::vector<std::vector<int>>& matrix);
// Function to load scanpaths from a file and return them as strings
std::vector<std::string> loadScanpathsFromFile(const std::string& filename);
// Function to get the number of scanpaths in an input string
int getNumberOfScanpaths(const std::string& input);

// Function to compute the encrypted distance matrix for a scanpath
bool getEncryptedDistanceMatrix(const std::string& scanpath, const std::unordered_map<std::string, std::vector<MPZ_Wrapper>>& encryptedMatrixMap, std::vector<MPZ_Wrapper>& encryptedDistanceMatrix);
// Function to store the encrypted distance matrix for a scanpath
void storeEncryptedDistanceMatrix(const std::vector<MPZ_Wrapper>& encryptedDistanceMatrix, const std::string& scanpath, std::unordered_map<std::string, std::vector<MPZ_Wrapper>>& encryptedMatrixMap);
// Function for logging messages
void log(const std::string& message_subject, const std::string& message,
         const std::string& log_file_path, const std::string& code_name,
         int iteration_id, const std::string& phase);


// Function to flush remaining log messages to a file
void flush_remaining_logs(const std::string& log_file_path);
// Function to revert the insertion of additional values into a vector
std::vector<MPZ_Wrapper> revertInsertionOfAdditionalValues(
    const std::vector<MPZ_Wrapper>& combinedVector,
    const std::vector<size_t>& insertionIndices);

// Function to convert a duration to a string with a specified precision
std::string durationToString(const std::chrono::duration<double>& duration, int precision);

// Function to generate a random mpz_t value from /dev/urandom
void random_from_urandom(mpz_t result, const mpz_t max_value);
// Function to generate a random mpz_t value with a specified number of bits from /dev/urandom
void random_from_urandom_bits(mpz_t result, unsigned int bits);
// Function to generate a random index from /dev/urandom in the range [0, max_value)
size_t random_index_from_urandom(size_t max_value);

// Function to send an encrypted result over a socket
void sendEncryptedResult(SocketComm& socketComm, const std::vector<MPZ_Wrapper>& resultVector, int flag);
// Function to start an Alice server and listen for connections on a specified port
void startAliceServerAndListen(SocketComm& socketComm, int portNumber);
// Function to send a size message over a socket
void sendSizeMessage(std::vector<std::string> allScanpaths, SocketComm& socketComm, int flag);

// Function to encrypt a vector of integers using the Paillier encryption scheme
std::vector<MPZ_Wrapper> encryptIntVector(const Paillier& paillier, const std::vector<int>& intVector);
// Function to encrypt and flatten the difference matrix for Alice's scanpath
std::vector<MPZ_Wrapper> encryptAndFlattenDifferenceMatrix(const Paillier& paillier, const std::string& aliceScanpath);

// Function to build a log message with a title and statistics
std::string buildLogMessage(const std::string& title, const OnlineStats& stats);




#endif // UTILS_HPP
