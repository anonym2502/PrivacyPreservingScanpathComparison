//
//  utils.cpp
//  PaillierCrypto
//
//  Created by Süleyman Özdel on 7.09.2023.
//

#include "utils.hpp"
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <set>
#include <fstream>
#include <ctime>
#include <iomanip>  // For std::put_time


// Default constructor for MPZ_Wrapper. Initializes an mpz_t value.
MPZ_Wrapper::MPZ_Wrapper() {
    mpz_init(value);
}

// Constructor that initializes the MPZ_Wrapper with an existing mpz_t value.
MPZ_Wrapper::MPZ_Wrapper(const mpz_t& other) {
    mpz_init_set(value, other);
}

// Copy constructor for MPZ_Wrapper. Initializes with another MPZ_Wrapper's value.
MPZ_Wrapper::MPZ_Wrapper(const MPZ_Wrapper& other) {
    mpz_init_set(value, other.value);
}

// Move constructor for MPZ_Wrapper. Transfers ownership of the mpz_t value.
MPZ_Wrapper::MPZ_Wrapper(MPZ_Wrapper&& other) noexcept {
    mpz_init_set(value, other.value);
    mpz_set_ui(other.value, 0); // Reset other's value.
}

// Destructor for MPZ_Wrapper. Clears the allocated mpz_t value.
MPZ_Wrapper::~MPZ_Wrapper() {
    mpz_clear(value);
}



// Assignment operator for assigning from an mpz_t value.
MPZ_Wrapper& MPZ_Wrapper::operator=(const mpz_t& other) {
    if (mpz_cmp(value, other) != 0) {
        mpz_set(value, other);
    }
    return *this;
}


// Assignment operator for assigning from another MPZ_Wrapper.
MPZ_Wrapper& MPZ_Wrapper::operator=(const MPZ_Wrapper& other) {
    if (this != &other) {
        mpz_set(value, other.value);
    }
    return *this;
}

// Move assignment operator for MPZ_Wrapper.
MPZ_Wrapper& MPZ_Wrapper::operator=(MPZ_Wrapper&& other) noexcept {
    if (this != &other) {
        mpz_set(value, other.value);
        mpz_set_ui(other.value, 0);
    }
    return *this;
}



// Matrix constructor: Initializes a matrix with specified number of rows and columns.
Matrix::Matrix(unsigned int r, unsigned int c) : rows(r), cols(c), data(r, std::vector<MPZ_Wrapper>(c)) {}

// Matrix copy constructor: Performs a deep copy of another matrix.
Matrix::Matrix(const Matrix& other) : rows(other.rows), cols(other.cols), data(other.data) {}

// Matrix assignment operator: Assigns this matrix with the values from another matrix.
Matrix& Matrix::operator=(const Matrix& other) {
    if (this == &other) return *this;
    rows = other.rows;
    cols = other.cols;
    data = other.data;
    return *this;
}

// Sets the value at a specific row and column in the matrix.
void Matrix::set(int row, int col, const mpz_t& value) {
    if (row < 0 || row >= rows || col < 0 || col >= cols) {
        throw std::out_of_range("Index out of range");
    }
    data[row][col] = value;
}

// Retrieves the value from a specific row and column in the matrix.
void Matrix::get(int row, int col, mpz_t& out_value) const {
    if (row < 0 || row >= rows || col < 0 || col >= cols) {
        throw std::out_of_range("Index out of range");
    }
    mpz_set(out_value, data[row][col].value);
}

// Returns the number of rows in the matrix.
int Matrix::getRows() const {
    return rows;
}

// Returns the number of columns in the matrix.
int Matrix::getCols() const {
    return cols;
}





// Constructor for OnlineStats. Initializes statistical measures to zero.
OnlineStats::OnlineStats() : n(0), mean(0), M2(0), min(0), max(0), totalSum(0) {}

// Updates statistical measures with a new duration value.
void OnlineStats::update(std::chrono::duration<double> x_duration) {
    double x = x_duration.count() * 1e6;  // Convert duration to microseconds.

    n += 1;  // Increment the count of measurements.
    double delta = x - mean;
    mean += delta / n;  // Update the running mean.
    M2 += delta * (x - mean);  // Update the sum of squares of differences from the current mean.
    totalSum += x;  // Update the total sum of measurements.
    
    // Update min and max values.
    if (n == 1) {
        min = x;
        max = x;
    } else {
        if (x < min) min = x;
        if (x > max) max = x;
    }
}

// Returns the count of measurements.
int OnlineStats::getCount() const {
    return n;
}

// Returns the mean of the measurements.
double OnlineStats::getMean() const {
    return mean;
}

// Returns the variance of the measurements.
double OnlineStats::getVariance() const {
    if (n < 2) return 0;  // Variance is undefined for n < 2.
    return M2 / n;  // Return the variance.
}


// Returns the standard deviation of the measurements.
double OnlineStats::getStandardDeviation() const {
    return std::sqrt(getVariance());
}

// Returns the minimum of the measurements.
double OnlineStats::getMin() const {
    return min;
}

// Returns the maximum of the measurements.
double OnlineStats::getMax() const {
    return max;
}

// Returns the total sum of all measurements.
double OnlineStats::getTotalSum() const {
    return totalSum;
}

// Resets all statistical measures to zero.
void OnlineStats::reset() {
    n = 0;
    mean = 0;
    M2 = 0;
    min = 0;
    max = 0;
    totalSum = 0;
}





// Constructor for Timer. Initializes the timer with the current time.
Timer::Timer() : start(std::chrono::high_resolution_clock::now()) {}

// Destructor for Timer. Does not need to perform any specific action.
Timer::~Timer() {
}

// Gets the duration between start and end times.
std::chrono::duration<double> Timer::getDuration() const {
    auto end = std::chrono::high_resolution_clock::now();
    duration = end - start;
    return duration;
}




// Stores the encrypted distance matrix in a map, using the scanpath as the key.
void storeEncryptedDistanceMatrix(const std::vector<MPZ_Wrapper>& encryptedDistanceMatrix, const std::string& scanpath, std::unordered_map<std::string, std::vector<MPZ_Wrapper>>& encryptedMatrixMap) {
    encryptedMatrixMap[scanpath] = encryptedDistanceMatrix;
}

// Retrieves the encrypted distance matrix from the map if it exists for the given scanpath.
bool getEncryptedDistanceMatrix(const std::string& scanpath, const std::unordered_map<std::string, std::vector<MPZ_Wrapper>>& encryptedMatrixMap, std::vector<MPZ_Wrapper>& encryptedDistanceMatrix) {
    auto it = encryptedMatrixMap.find(scanpath);
    if (it != encryptedMatrixMap.end()) {
        encryptedDistanceMatrix = it->second;
        return true;
    } else {
        encryptedDistanceMatrix.clear();
        return false;
    }
}

// Parses the number of scanpaths from a given input string.
int getNumberOfScanpaths(const std::string& input) {
    std::size_t commaPos = input.find(',');
    if (commaPos == std::string::npos) {
        std::cerr << "No comma found in input string." << std::endl;
        return -1;
    }
    std::string numberString = input.substr(commaPos + 1);
    return std::stoi(numberString);
}

// Creates a mapping from each character of the alphabet to an integer.
std::map<char, int> createAlphabetMapping() {
    std::map<char, int> alphabet;
    std::string allLetters = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
    for (int i = 0; i < allLetters.length(); ++i) {
        alphabet.insert(std::make_pair(allLetters[i], i));
    }
    return alphabet;
}

// Maps a string to a sequence of integers based on the provided alphabet mapping.
std::vector<int> mapStringUsingAlphabet(const std::map<char, int>& alphabet, const std::string& input) {
    std::vector<int> mappedString;
    for (char c : input) {
        auto it = alphabet.find(c);
        if (it != alphabet.end()) {
            mappedString.push_back(it->second);
        }
    }
    return mappedString;
}



// Randomly permutes the elements of a vector.
void randomlyPermutateVector(std::vector<MPZ_Wrapper>& vec) {
    size_t n = vec.size();
    for (size_t i = n - 1; i > 0; --i) {
        size_t j = random_index_from_urandom(i + 1);
        std::swap(vec[i], vec[j]);
    }
}

// Appends two vectors and returns the resulting vector.
std::vector<MPZ_Wrapper> appendVectors(const std::vector<MPZ_Wrapper>& vec1, const std::vector<MPZ_Wrapper>& vec2) {
    std::vector<MPZ_Wrapper> result;
    result.reserve(vec1.size() + vec2.size());
    result.insert(result.end(), vec1.begin(), vec1.end());
    result.insert(result.end(), std::make_move_iterator(vec2.begin()), std::make_move_iterator(vec2.end()));
    return result;
}

// Converts a string representation of MPZ_Wrappers into a vector of MPZ_Wrappers.
std::vector<MPZ_Wrapper> convertStringToMPZ_Wrapper(const std::string& str) {
    std::vector<MPZ_Wrapper> result;
    std::stringstream ss(str);
    std::string item;

    while (std::getline(ss, item, ',')) {
        MPZ_Wrapper wrapper;
        mpz_set_str(wrapper.value, item.c_str(), 10);
        result.push_back(wrapper);
    }

    return result;
}

// Loads scanpaths from a file into a vector.
std::vector<std::string> loadScanpathsFromFile(const std::string& filename) {
    std::vector<std::string> scanpaths;
    std::ifstream file(filename);
    std::string line;

    while (std::getline(file, line)) {
        scanpaths.push_back(line);
    }

    return scanpaths;
}

// Converts an encrypted result (vector of MPZ_Wrappers) to a string representation.
std::string encryptedResultToString(const std::vector<MPZ_Wrapper>& encryptedMappedResult) {
    std::ostringstream os;
    for (size_t i = 0; i < encryptedMappedResult.size(); i++) {
        char buffer[8192];
        gmp_sprintf(buffer, "%Zd", encryptedMappedResult[i].value);
        os << buffer;
        if (i != encryptedMappedResult.size() - 1) {
            os << ", ";
        }
    }
    return os.str();
}

// Adds an operation flag to the beginning of a string.
void addOperationFlag(int number, std::string& str) {
    str = std::to_string(number) + "," + str;
}

// Extracts the operation flag from a string.
int getOperationFlag(const std::string input) {
    std::stringstream ss(input);
    int firstValue;
    ss >> firstValue;
    return firstValue;
}




// Converts a duration to a string with specified precision.
std::string durationToString(const std::chrono::duration<double>& duration, int precisionVal = 8) {
    std::ostringstream oss;
    auto microsecs = std::chrono::duration_cast<std::chrono::microseconds>(duration);
    oss << std::fixed << std::setprecision(precisionVal) << microsecs.count();
    return oss.str();
}

// Builds a log message with various statistics.
std::string buildLogMessage(const std::string& title, const OnlineStats& stats) {
    std::ostringstream logMessage;
    logMessage << title << " Time: "
               << "Mean: " <<std::fixed << std::setprecision(2) << stats.getMean() << ", "
               << "Min: " <<std::fixed << std::setprecision(2) << stats.getMin() << ", "
               << "Max: " <<std::fixed << std::setprecision(2) << stats.getMax() << ", "
               << "Standard Deviation: " << std::fixed << std::setprecision(2) <<stats.getStandardDeviation() << ", "
               << "Total Time: " << std::fixed << std::setprecision(2) <<stats.getTotalSum() << ", "
                << "Count: " << std::fixed << std::setprecision(2) << stats.getCount() ;

    return logMessage.str();
}

// Computes the difference matrix between a mapped scanpath and an alphabet mapping.
std::vector<std::vector<int>> computeAlphabetDifferenceMatrix(const std::vector<int>& mappedScanpath, const std::map<char, int>& alphabet) {
    int n = mappedScanpath.size();
    int m = alphabet.size();
    std::vector<std::vector<int>> matrix(n, std::vector<int>(m));

    for (int i = 0; i < n; ++i) {
        int rowValue = mappedScanpath[i];
        int j = 0;
        for (const auto& [charKey, charValue] : alphabet) {
            if(rowValue - charValue != 0){
                matrix[i][j] = 1;
            }
            else {
                matrix[i][j] = 0;
            }
            ++j;
        }
    }

    return matrix;
}

// Flattens a 2D matrix into a 1D vector.
std::vector<int> flattenMatrix(const std::vector<std::vector<int>>& matrix) {
    std::vector<int> flattenedVector;
    
    // Reserve space in the vector for efficiency
    int totalSize = 0;
    for (const auto& row : matrix) {
        totalSize += row.size();
    }
    flattenedVector.reserve(totalSize);
    
    // Flatten the matrix into the vector
    for (const auto& row : matrix) {
        flattenedVector.insert(flattenedVector.end(), row.begin(), row.end());
    }

    return flattenedVector;
}

// Logs a message with additional information such as timestamp, code name, iteration, and phase.
void log(const std::string& message_subject, const std::string& message,
         const std::string& log_file_path, const std::string& code_name,
         int iteration_id, const std::string& phase) {
    
    static std::vector<std::string> log_buffer;
    
    // Get current time point
    auto now = std::chrono::system_clock::now();
    auto epoch_time = now.time_since_epoch();
    
    // Microseconds since epoch
    auto microsecs_since_epoch = std::chrono::duration_cast<std::chrono::microseconds>(epoch_time).count();
    
    // Convert to time_t for UTC time
    auto now_as_time_t = std::chrono::system_clock::to_time_t(now);
    auto utc_time = gmtime(&now_as_time_t);
    
    std::ostringstream oss;
    oss << "UTC Time: " << std::put_time(utc_time, "%Y-%m-%d %H:%M:%S")
        << ", " << microsecs_since_epoch << " ms since epoch"
        << ", Code: " << code_name
        << ", Iteration: " << iteration_id
        << ", Phase: " << phase
        << ", Subject: " << message_subject
        << " - " << message << "\n";
    
    log_buffer.push_back(oss.str());
    
    // Conditionally flush - example: every 10 logs
    if(log_buffer.size() >= 1) {
        std::ofstream log_file(log_file_path, std::ios_base::out | std::ios_base::app);
        if(!log_file.is_open()) {
            std::cerr << "Failed to open log file." << std::endl;
            return;
        }
        
        for(const auto& log_message : log_buffer) {
            log_file << log_message;
        }
        
        log_buffer.clear();
    }
}

// Flushes any remaining log messages to a log file.
void flush_remaining_logs(const std::string& log_file_path) {
    static std::vector<std::string> log_buffer;  // Re-declared for demonstration. Remove in actual code.
    // Write any remaining logs
    std::ofstream log_file(log_file_path, std::ios_base::out | std::ios_base::app);
    if(!log_file.is_open()) {
        std::cerr << "Failed to open log file." << std::endl;
        return;
    }
    
    for(const auto& log_message : log_buffer) {
        log_file << log_message;
    }
    
    log_buffer.clear();
}

// Generates a random number within a specified range using /dev/urandom.
void random_from_urandom(mpz_t result, const mpz_t max_value) {
    size_t bytes = (mpz_sizeinbase(max_value, 2) + 7) / 8; // Determine bytes based on max_value
    unsigned char *buf = new unsigned char[bytes];
    
    std::ifstream urandom("/dev/urandom", std::ios::in|std::ios::binary);
    if (!urandom) {
        throw std::runtime_error("Failed to open /dev/urandom");
    }
    urandom.read(reinterpret_cast<char*>(buf), bytes);
    urandom.close();

    mpz_import(result, bytes, 1, 1, 0, 0, buf);
    delete[] buf;

    // Ensure the result is within the range [0, max_value)
    mpz_mod(result, result, max_value);
}

// Generates a random number with a specified number of bits using /dev/urandom.
void random_from_urandom_bits(mpz_t result, unsigned int bits) {
    // Convert bits to bytes, then read from /dev/urandom
    unsigned int size_in_bytes = (bits + 7) / 8; // +7 to round up to the nearest byte

    std::vector<unsigned char> buffer(size_in_bytes);
    std::ifstream urandom("/dev/urandom", std::ios::in | std::ios::binary);

    if (!urandom) {
        throw std::runtime_error("Failed to open /dev/urandom");
    }

    urandom.read(reinterpret_cast<char*>(buffer.data()), size_in_bytes);
    urandom.close();

    mpz_import(result, size_in_bytes, 1, sizeof(buffer[0]), 0, 0, buffer.data());

    // Now, since we got a random number of size_in_bytes length, we might have a number
    // larger than our desired bit-length. So, we'll take it modulo 2^bits to ensure it's within range.
    mpz_t two_power_bits;
    mpz_init(two_power_bits);
    mpz_ui_pow_ui(two_power_bits, 2, bits);

    mpz_mod(result, result, two_power_bits);

    mpz_clear(two_power_bits);
}

// Generates a random index within a specified range using /dev/urandom.
size_t random_index_from_urandom(size_t max_value) {
    mpz_t result, mpz_max;
    mpz_init(result);
    mpz_init_set_ui(mpz_max, max_value);

    random_from_urandom(result, mpz_max);

    size_t index = mpz_get_ui(result);

    mpz_clear(result);
    mpz_clear(mpz_max);

    return index;
}



// Sends an encrypted result through a socket communication.
void sendEncryptedResult(SocketComm& socketComm, const std::vector<MPZ_Wrapper>& resultVector, int flag) {
    std::string resultString = encryptedResultToString(resultVector);
    addOperationFlag(flag, resultString);
    socketComm.send_message(resultString);
}

// Starts an Alice server and waits for Bob to connect.
void startAliceServerAndListen(SocketComm& socketComm,  int portNumber) {
    if (!socketComm.init_as_server(portNumber)) {
        std::cerr << "Error initializing server." << std::endl;
    }
    std::cout << "Waiting for Bob to connect..." << std::endl;
    
    if (!socketComm.listen_for_client()) {
        std::cerr << "Error listening for client." << std::endl;
    }
    std::cout << "Bob connected!" << std::endl;
}

// Sends a size message through socket communication.
void sendSizeMessage(std::vector<std::string> allScanpaths, SocketComm& socketComm, int flag) {
    std::string sizeMessage = std::to_string(allScanpaths.size());
    addOperationFlag(flag, sizeMessage);

    if (!socketComm.send_message(sizeMessage)) {
        std::cerr << "Error sending message." << std::endl;
    }
}

// Encrypts a vector of integers using Paillier encryption.
std::vector<MPZ_Wrapper> encryptIntVector(const Paillier& paillier, const std::vector<int>& intVector) {
    std::vector<MPZ_Wrapper> encryptedResult(intVector.size());

    for (size_t i = 0; i < intVector.size(); ++i) {
        mpz_t plaintext;
        mpz_t ciphertext;
        mpz_init_set_si(plaintext, intVector[i]);
        mpz_init(ciphertext);
        paillier.encryptMessage(ciphertext, plaintext);
        encryptedResult[i] = ciphertext;

        mpz_clear(plaintext);
        mpz_clear(ciphertext);
    }

    return encryptedResult;
}

// Encrypts and flattens a difference matrix using Paillier encryption.
std::vector<MPZ_Wrapper> encryptAndFlattenDifferenceMatrix(const Paillier& paillier, const std::string& aliceScanpath) {
    // Create alphabet mapping and map Alice's scanpath
    std::map<char, int> alphabet = createAlphabetMapping();
    std::vector<int> mappedResult = mapStringUsingAlphabet(alphabet, aliceScanpath);

    // Compute the difference matrix and then flatten it
    std::vector<std::vector<int>> differenceMatrix = computeAlphabetDifferenceMatrix(mappedResult, alphabet);
    std::vector<int> flattenedMatrix = flattenMatrix(differenceMatrix);

    // Encrypt the flattened matrix and return
    return encryptIntVector(paillier, flattenedMatrix);
}
