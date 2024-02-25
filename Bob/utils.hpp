// utils.hpp
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
#include "PublicPaillier.hpp"
#include "SocketComm.hpp"


// Declaration of MPZ_Wrapper class, a wrapper for GMP's mpz_t type
class MPZ_Wrapper {
public:
    // Default constructor
    MPZ_Wrapper();

    // Constructor from a GMP mpz_t type
    MPZ_Wrapper(const mpz_t& other);

    // Copy constructor
    MPZ_Wrapper(const MPZ_Wrapper& other);

    // Move constructor
    MPZ_Wrapper(MPZ_Wrapper&& other) noexcept;

    // Destructor
    ~MPZ_Wrapper();

    // Copy assignment operator from mpz_t type
    MPZ_Wrapper& operator=(const mpz_t& other);

    // Copy assignment operator from another MPZ_Wrapper
    MPZ_Wrapper& operator=(const MPZ_Wrapper& other);

    // Move assignment operator
    MPZ_Wrapper& operator=(MPZ_Wrapper&& other) noexcept;

    // The wrapped mpz_t value
    mpz_t value;
};





// Matrix class definition
class Matrix {
private:
    // Number of rows and columns in the matrix
    unsigned int rows, cols;
    
    // 2D vector storing the matrix data using MPZ_Wrapper for each element
    std::vector<std::vector<MPZ_Wrapper>> data;

public:
    // Constructor for creating a matrix of specific dimensions
    Matrix(unsigned int r, unsigned int c);

    // Copy constructor for deep copying another matrix
    Matrix(const Matrix& other);

    // Copy assignment operator for deep copying another matrix
    Matrix& operator=(const Matrix& other);

    // Destructor for the matrix
    ~Matrix();

    // Set a specific element in the matrix
    void set(int row, int col, const mpz_t& value);

    // Get a specific element from the matrix
    void get(int row, int col, mpz_t& out_value) const;

    // Get the number of rows in the matrix
    int getRows() const;

    // Get the number of columns in the matrix
    int getCols() const;
};



// OnlineStats class for maintaining and updating statistical measures
class OnlineStats {
public:
    // Constructor
    OnlineStats();

    // Updates statistics with a new value
    void update(std::chrono::duration<double> x);

    // Returns the number of values processed
    int getCount() const;

    // Returns the mean (average) of the values
    double getMean() const;

    // Returns the variance of the values
    double getVariance() const;

    // Returns the standard deviation of the values
    double getStandardDeviation() const;

    // Returns the minimum value encountered
    double getMin() const;

    // Returns the maximum value encountered
    double getMax() const;

    // Returns the total sum of all values
    double getTotalSum() const;

    // Resets all statistics to initial state
    void reset();

private:
    int n;           // Number of values processed
    double mean;     // Mean of the values
    double M2;       // Sum of squares of differences from the current mean
    double min;      // Minimum value
    double max;      // Maximum value
    double totalSum; // Total sum of all values
};



// Timer class for measuring time intervals
class Timer {
public:
    // Constructor that starts the timer
    Timer();

    // Destructor
    ~Timer();

    // Returns the duration from the start of the timer to the current moment
    std::chrono::duration<double> getDuration() const;

private:
    // Time point representing the start time
    std::chrono::time_point<std::chrono::high_resolution_clock> start;

    // Duration representing the measured time interval
    mutable std::chrono::duration<double> duration;
};



// Function declarations


// Declaration of the appendVectors function
// Appends two vectors of MPZ_Wrapper elements into a single vector

std::vector<MPZ_Wrapper> appendVectors(const std::vector<MPZ_Wrapper>& vec1, const std::vector<MPZ_Wrapper>& vec2);

// Declaration of the convertStringToMPZ_Wrapper function
// Converts a comma-separated string into a vector of MPZ_Wrapper objects

std::vector<MPZ_Wrapper> convertStringToMPZ_Wrapper(const std::string& str);



// Declaration of the randomlyPermutateVector function
// Randomly permutes the elements of a vector of MPZ_Wrapper objects.
// This function is used for cryptographically securely shuffling the elements of the vector.
void randomlyPermutateVector(std::vector<MPZ_Wrapper>& vec);

// Declaration of the encryptScanpath function
// Encrypts a scanpath string using the PublicPaillier encryption scheme.
// Each character in the string is mapped to an integer, and each integer is encrypted.
// The function returns a vector of encrypted MPZ_Wrapper objects.
std::vector<MPZ_Wrapper> encryptScanpath(const PublicPaillier& paillier, const std::string& bobScanpath);


// Declaration of the unflattenVector function
// Converts a 1D vector (flattened vector) into a 2D vector (matrix) of specified width.
// The height of the matrix is determined based on the size of the flattened vector and the given width.
std::vector<std::vector<MPZ_Wrapper>> unflattenVector(const std::vector<MPZ_Wrapper>& flattenedVector, int width);

// Declaration of the log function
// Logs a message with additional context information to a specified log file.
// The message includes the subject, the actual message, the path of the log file, a code name (identifier),
// an iteration ID, and the phase during which the log is made.
void log(const std::string& message_subject, const std::string& message,
         const std::string& log_file_path, const std::string& code_name,
         int iteration_id, const std::string& phase);



// Declaration of the encryptedResultToString function
// Converts a vector of MPZ_Wrapper objects (each representing an encrypted value) into a comma-separated string.
// This function is useful for serializing encrypted data for storage or transmission.
std::string encryptedResultToString(const std::vector<MPZ_Wrapper>& encryptedMappedResult);



// Declaration of the sendEncryptedResult function
// Sends an encrypted result vector over a socket connection.
// This function converts the encrypted data into a string, adds an operation flag, and sends it using the SocketComm object.
void sendEncryptedResult(SocketComm& socketComm, const std::vector<MPZ_Wrapper>& resultVector, int flag);



// Declaration of the createAlphabetMapping function
// Creates a map from characters to integers for the alphabet
std::map<char, int> createAlphabetMapping();



// Declaration of the mapString function
// Maps each character of the input string to its corresponding integer value based on the provided alphabet mapping
std::vector<int> mapString(const std::map<char, int>& alphabet, const std::string& input);



    
// Declaration of the loadScanpathsFromFile function
// Loads scanpaths (each line representing one scanpath) from a file and returns them as a vector of strings.
std::vector<std::string> loadScanpathsFromFile(const std::string& filename);

// Declaration of the getNumberOfScanpaths function
// Extracts and returns the number of scanpaths from a given input string, typically from a network message.
int getNumberOfScanpaths(const std::string& input);

// Declaration of the addOperationFlag function
// Prepends an operation flag (as an integer) to a string. Useful for encoding messages with a specific operation code.
void addOperationFlag(int number, std::string& str);

// Declaration of the getOperationFlag function
// Extracts and returns the operation flag (as an integer) from a string. Typically used to decode the type of operation from a network message.
int getOperationFlag(const std::string input);



// Declaration of the loadScanpathsFromFile function
// Loads scanpaths (each line representing one scanpath) from a file and returns them as a vector of strings.
std::vector<std::string> loadScanpathsFromFile(const std::string& filename);

// Declaration of the getNumberOfScanpaths function
// Extracts and returns the number of scanpaths from a given input string, typically from a network message.
int getNumberOfScanpaths(const std::string& input);

// Declaration of the addOperationFlag function
// Prepends an operation flag (as an integer) to a string. Useful for encoding messages with a specific operation code.
void addOperationFlag(int number, std::string& str);

// Declaration of the getOperationFlag function
// Extracts and returns the operation flag (as an integer) from a string. Typically used to decode the type of operation from a network message.
int getOperationFlag(const std::string input);




// Declaration of flush_remaining_logs
// Flushes any remaining logs in the buffer to the specified log file.
void flush_remaining_logs(const std::string& log_file_path);

// Declaration of durationToString
// Converts a std::chrono::duration to a string with specified precision.
std::string durationToString(const std::chrono::duration<double>& duration, int precisionVal);

// Declaration of buildLogMessage
// Constructs a log message string from a title and the statistics of an OnlineStats object.
std::string buildLogMessage(const std::string& title, const OnlineStats& stats);


// Declaration of random_from_urandom
// Generates a random number using the /dev/urandom device and ensures it is less than max_value.
void random_from_urandom(mpz_t result, const mpz_t max_value);

// Declaration of random_from_urandom_bits
// Generates a random number of a specified bit length using the /dev/urandom device.
void random_from_urandom_bits(mpz_t result, unsigned int bits);

// Declaration of random_double_from_urandom
// Generates a random double precision floating-point number using /dev/urandom.
double randomDoubleFromUrandom();

// Declaration of random_index_from_urandom
// Generates a random index (as size_t) less than a specified max_value using /dev/urandom.
size_t random_index_from_urandom(size_t max_value);
void print_mpz_alloc(const char* name, const mpz_t& mpzVar);




#endif // UTILS_HPP
