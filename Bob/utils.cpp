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
#include "PublicPaillier.hpp"
#include <gmp.h>
#include <iomanip>  // Include this header for setprecision


// Function to print the allocation metadata of an mpz_t variable
void print_mpz_alloc(const char* name, const mpz_t& mpzVar) {
    std::cout << name << " allocation:" << std::endl;
    std::cout << "  _mp_alloc = " << mpzVar->_mp_alloc << std::endl;
    std::cout << "  _mp_size  = " << mpzVar->_mp_size << std::endl;
    // If you also want to print whether the _mp_d pointer is NULL or not
    std::cout << "  _mp_d pointer is " << (mpzVar->_mp_d == NULL ? "NULL" : "non-NULL") << std::endl;
}


// Default constructor: initializes the mpz_t value
MPZ_Wrapper::MPZ_Wrapper() {
    mpz_init(value);
}

// Constructor from a GMP mpz_t type: initializes and sets value
MPZ_Wrapper::MPZ_Wrapper(const mpz_t& other) {
    mpz_init_set(value, other);
}

// Copy constructor: initializes and sets value from another MPZ_Wrapper
MPZ_Wrapper::MPZ_Wrapper(const MPZ_Wrapper& other) {
    mpz_init_set(value, other.value);
}

// Move constructor: initializes and swaps values
MPZ_Wrapper::MPZ_Wrapper(MPZ_Wrapper&& other) noexcept {
    mpz_init(value);
    mpz_swap(value, other.value);
}

// Destructor: clears the mpz_t value to deallocate memory
MPZ_Wrapper::~MPZ_Wrapper() {
    mpz_clear(value);
}

// Copy assignment operator from mpz_t type: sets value if different
MPZ_Wrapper& MPZ_Wrapper::operator=(const mpz_t& other) {
    if (mpz_cmp(value, other) != 0) {
        mpz_set(value, other);
    }
    return *this;
}

// Copy assignment operator from another MPZ_Wrapper: sets value if different
MPZ_Wrapper& MPZ_Wrapper::operator=(const MPZ_Wrapper& other) {
    if (this != &other) {
        mpz_set(value, other.value);
    }
    return *this;
}

// Move assignment operator: swaps values and clears the other
MPZ_Wrapper& MPZ_Wrapper::operator=(MPZ_Wrapper&& other) noexcept {
    if (this != &other) {
        mpz_swap(value, other.value);
    }
    mpz_clear(other.value); // Clear other.value to prevent memory leak
    return *this;
}




// Constructor for initializing the matrix with specific dimensions
Matrix::Matrix(unsigned int r, unsigned int c) : rows(r), cols(c), data(r, std::vector<MPZ_Wrapper>(c)) {
    // Initialization of mpz_t values is now handled by MPZ_Wrapper's constructor.
}

// Copy constructor for deep copying another matrix
Matrix::Matrix(const Matrix& other) : rows(other.rows), cols(other.cols), data(other.data) {}

// Assignment operator for deep copying another matrix
Matrix& Matrix::operator=(const Matrix& other) {
    if(this == &other) return *this;  // Check for self-assignment

    rows = other.rows;
    cols = other.cols;
    data = other.data;

    return *this;
}

// Set a specific element in the matrix
void Matrix::set(int row, int col, const mpz_t& value) {
    if(row < 0 || row >= rows || col < 0 || col >= cols) {
        throw std::out_of_range("Index out of range"); // Check for valid index
    }
    data[row][col] = value;
}

// Get a specific element from the matrix
void Matrix::get(int row, int col, mpz_t& out_value) const {
    if(row < 0 || row >= rows || col < 0 || col >= cols) {
        throw std::out_of_range("Index out of range"); // Check for valid index
    }
    mpz_set(out_value, data[row][col].value);
}

// Get the number of rows in the matrix
int Matrix::getRows() const {
    return rows;
}

// Get the number of columns in the matrix
int Matrix::getCols() const {
    return cols;
}

// Destructor for the matrix
Matrix::~Matrix() {
    // The destructor does not need to do anything special here, because
    // MPZ_Wrapper takes care of cleaning up its own resources and
    // std::vector cleans up its resources automatically.
}






// Constructor initializing all member variables
OnlineStats::OnlineStats() : n(0), mean(0), M2(0), min(0), max(0), totalSum(0) {}

// Update the statistics with a new value
void OnlineStats::update(std::chrono::duration<double> x_duration) {
    double x = x_duration.count() * 1e6;  // converting duration to microseconds
    n += 1;  // Incrementing the count of values
    double delta = x - mean;  // Difference from the current mean
    mean += delta / n;  // Updating the mean
    M2 += delta * (x - mean);  // Updating the sum of squares of differences
    totalSum += x;  // Updating the total sum

    // Updating min and max values
    if (n == 1) {
        min = max = x;  // Set the first value as both min and max
    } else {
        if (x < min) min = x;  // Update min if the new value is smaller
        if (x > max) max = x;  // Update max if the new value is larger
    }
}

// Get the count of values processed
int OnlineStats::getCount() const {
    return n;
}

// Get the mean of the values
double OnlineStats::getMean() const {
    return mean;
}

// Get the variance of the values
double OnlineStats::getVariance() const {
    if (n < 2) return 0;  // Variance is not defined for n < 2
    return M2 / n;  // Variance is the mean of the squared deviations
}

// Get the standard deviation of the values
double OnlineStats::getStandardDeviation() const {
    return std::sqrt(getVariance());  // Standard deviation is the square root of variance
}

// Get the minimum value processed
double OnlineStats::getMin() const {
    return min;
}

// Get the maximum value processed
double OnlineStats::getMax() const {
    return max;
}

// Get the total sum of all values processed
double OnlineStats::getTotalSum() const {
    return totalSum;
}

// Reset all statistics to initial state
void OnlineStats::reset() {
    n = 0;
    mean = 0;
    M2 = 0;
    min = 0;
    max = 0;
    totalSum = 0;
}





// Constructor: Initializes and starts the timer
Timer::Timer() : start(std::chrono::high_resolution_clock::now()) {
    // 'start' is set to the current time point using high resolution clock
}

// Destructor: Currently empty as there's no dynamic allocation or cleanup required
Timer::~Timer() {
}

// Returns the duration from the start of the timer to the current moment
std::chrono::duration<double> Timer::getDuration() const {
    // Capturing the current time point
    auto end = std::chrono::high_resolution_clock::now();

    // Calculating the duration from start to end
    duration = end - start;  // 'duration' is the difference between end and start

    // Returning the duration as a double representing seconds
    return duration;
}




// appendVectors function:
// Merges two vectors of MPZ_Wrapper elements into one.
// The elements of vec1 are appended at the beginning of vec2.

std::vector<MPZ_Wrapper> appendVectors(const std::vector<MPZ_Wrapper>& vec1,
                                       std::vector<MPZ_Wrapper> vec2) {
    // Inserting all elements of vec1 at the beginning of vec2
    vec2.insert(vec2.begin(), vec1.begin(), vec1.end());

    // Returning the merged vector
    return vec2;
}

// convertStringToMPZ_Wrapper function:
// Converts a string containing comma-separated values into a vector of MPZ_Wrapper objects.
// Each value in the string is converted to an MPZ_Wrapper object.

std::vector<MPZ_Wrapper> convertStringToMPZ_Wrapper(const std::string& str) {
    // Vector to store the MPZ_Wrapper objects
    std::vector<MPZ_Wrapper> result;

    // Stringstream used for parsing the input string
    std::stringstream ss(str);

    // String used to temporarily store each comma-separated item
    std::string item;

    // Parsing the string, splitting by commas
    while (std::getline(ss, item, ',')) {
        // Creating an MPZ_Wrapper object
        MPZ_Wrapper wrapper;

        // Setting the value of the wrapper using the item string
        mpz_set_str(wrapper.value, item.c_str(), 10);  // Base 10 for decimal numbers

        // Adding the wrapper to the result vector
        result.push_back(wrapper);
    }

    // Clearing the stringstream and temporary string
    ss.str(std::string());
    ss.clear(); // Clears the error state flags of the stringstream, not the content.
    item.clear();
    
    // Returning the vector of MPZ_Wrapper objects
    return result;
}



// randomlyPermutateVector function:
// Randomly permutes the elements of a vector of MPZ_Wrapper objects.
// This is typically used for shuffling the elements in a secure manner.

void randomlyPermutateVector(std::vector<MPZ_Wrapper>& vec) {
    // Looping backwards through the vector
    for (size_t i = vec.size() - 1; i > 0; i--) {
        // Generating a random index up to i (inclusive)
        size_t j = random_index_from_urandom(i + 1); // Function to generate a random index

        // Swapping the element at i with the element at the random index j
        std::swap(vec[i], vec[j]);
    }
}

// encryptScanpath function:
// Encrypts a string (representing a scanpath) using the PublicPaillier encryption scheme.
// Each character in the string is mapped to an integer, and each integer is then encrypted.

std::vector<MPZ_Wrapper> encryptScanpath(const PublicPaillier& paillier, const std::string& bobScanpath) {
    // Creating a mapping of the alphabet to integers
    std::map<char, int> alphabet = createAlphabetMapping();

    // Mapping the scanpath string to a vector of integers
    std::vector<int> mappedResult = mapString(alphabet, bobScanpath);

    // Vector to store the encrypted result
    std::vector<MPZ_Wrapper> encryptedMappedResult(mappedResult.size());

    // Encrypting each integer in the mapped result
    for(size_t i = 0; i < mappedResult.size(); ++i) {
        mpz_t plaintext;
        mpz_t ciphertext;

        // Initializing GMP integers and setting plaintext to the mapped integer
        mpz_init_set_si(plaintext, mappedResult[i]);
        mpz_init(ciphertext);

        // Encrypting the plaintext
        paillier.encrypt(ciphertext, plaintext);

        // Storing the encrypted value in the result vector
        encryptedMappedResult[i] = ciphertext;

        // Clearing the GMP integers to free memory
        mpz_clear(plaintext);
        mpz_clear(ciphertext);
    }

    // Returning the vector of encrypted values
    return encryptedMappedResult;
}




// unflattenVector function:
// Converts a flattened 1D vector into a 2D matrix-like vector of specified width.
// The height of the resulting matrix is automatically determined.
std::vector<std::vector<MPZ_Wrapper>> unflattenVector(const std::vector<MPZ_Wrapper>& flattenedVector, int width) {
    // Calculating the height of the matrix based on the size of the flattened vector and the width
    int height = flattenedVector.size() / width;

    // Creating a 2D vector (matrix) with the calculated dimensions
    std::vector<std::vector<MPZ_Wrapper>> matrix(height, std::vector<MPZ_Wrapper>(width));

    // Populating the matrix with elements from the flattened vector
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            matrix[i][j] = flattenedVector[i * width + j];
        }
    }

    // Returning the 2D matrix
    return matrix;
}





// encryptedResultToString function:
// Converts a vector of encrypted MPZ_Wrapper objects into a string.
// Each encrypted value is converted to a string and separated by commas.

std::string encryptedResultToString(const std::vector<MPZ_Wrapper>& encryptedMappedResult) {
    // Using a string stream to build the resulting string
    std::ostringstream os;

    // Iterating over each encrypted MPZ_Wrapper object in the vector
    for (size_t i = 0; i < encryptedMappedResult.size(); i++) {
        // Buffer to store the string representation of the encrypted value
        char buffer[8192];  // Assuming the encrypted value can fit into this buffer; adjust size if needed

        // Converting the encrypted value to a string using GMP's sprintf function
        gmp_sprintf(buffer, "%Zd", encryptedMappedResult[i].value);

        // Appending the string representation to the string stream
        os << buffer;

        // Adding a comma separator between values, except after the last value
        if (i != encryptedMappedResult.size() - 1) {
            os << ", ";
        }
    }

    // Returning the resulting comma-separated string
    return os.str();
}


// sendEncryptedResult function:
// Sends an encrypted result vector through socket communication.
// The result vector is first converted to a string, and an operation flag is added before sending.

void sendEncryptedResult(SocketComm& socketComm, const std::vector<MPZ_Wrapper>& resultVector, int flag) {
    // Converting the vector of encrypted MPZ_Wrapper objects to a string
    std::string resultString = encryptedResultToString(resultVector);

    // Adding an operation flag to the string. The function 'addOperationFlag' presumably modifies the string to include the flag.
    addOperationFlag(flag, resultString);

    // Sending the modified string (containing the encrypted result and the flag) through the socket
    socketComm.send_message(resultString);
}





// createAlphabetMapping function:
// Creates and returns a mapping of each letter in the alphabet to a unique integer.
std::map<char, int> createAlphabetMapping() {
    // Creating a map to store the character-to-integer mapping
    std::map<char, int> alphabet;

    // String containing all uppercase and lowercase letters of the English alphabet
    std::string allLetters = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";

    // Looping over each character in the allLetters string
    for (int i = 0; i < allLetters.length(); ++i) {
        // Inserting a pair of the letter and its position into the map
        // The position is zero-based, so 'A' maps to 0, 'B' to 1, and so on
        alphabet.insert(std::make_pair(allLetters[i], i));
    }

    // Clearing the allLetters string (optional, as it's going out of scope)
    allLetters.clear();

    // Returning the completed map
    return alphabet;
}




// mapString function:
// Maps each character in the given string to an integer based on the provided alphabet mapping.
// Characters not present in the alphabet mapping will be skipped.

std::vector<int> mapString(const std::map<char, int>& alphabet, const std::string& input) {
    // Vector to store the mapped integers
    std::vector<int> mappedString;

    // Looping through each character in the input string
    for (char c : input) {
        // Searching for the character in the alphabet map
        auto it = alphabet.find(c);

        // If the character is found in the map
        if (it != alphabet.end()) {
            // Adding the corresponding integer value to the mappedString vector
            mappedString.push_back(it->second);
        }
        // Note: Characters not found in the alphabet map will be skipped.
        // If handling of such characters is needed, additional code should be added here.
    }

    // Returning the vector containing the mapped integers
    return mappedString;
}



// loadScanpathsFromFile function:
// Reads a file and loads each line (representing a scanpath) into a vector of strings.

std::vector<std::string> loadScanpathsFromFile(const std::string& filename) {
    std::vector<std::string> scanpaths;  // Vector to store the scanpaths
    std::ifstream file(filename);  // Opening the file
    std::string line;

    // Reading each line from the file and adding it to the vector
    while (std::getline(file, line)) {
        scanpaths.push_back(line);
    }

    // Returning the vector of scanpaths
    return scanpaths;
}

// getNumberOfScanpaths function:
// Parses a string to extract the number of scanpaths. The number is expected to be after a comma.

int getNumberOfScanpaths(const std::string& input) {
    std::size_t commaPos = input.find(',');  // Find the position of the comma
    if (commaPos == std::string::npos) {  // Check if a comma was found
        // Handle the error - maybe return a special value, or throw an exception
        std::cerr << "No comma found in input string." << std::endl;
        return -1;  // or some other error-indicating value
    }
    std::string numberString = input.substr(commaPos + 1);  // Get the substring after the comma
    return std::stoi(numberString);  // Convert substring to integer and return
}

// addOperationFlag function:
// Adds an operation flag (as an integer) at the beginning of a string, separated by a comma.

void addOperationFlag(int number, std::string& str) {
    str = std::to_string(number) + "," + str;  // Prepending the number and a comma to the string
}

// getOperationFlag function:
// Extracts the operation flag (as an integer) from the beginning of a string.

int getOperationFlag(const std::string input) {
    std::stringstream ss(input);
    int firstValue;
    ss >> firstValue;  // Parse the first integer value
    return firstValue;  // Returning the extracted operation flag
}







// log function:
// Buffers a log message with contextual information and writes it to a file when the buffer is full.

void log(const std::string& message_subject, const std::string& message,
         const std::string& log_file_path, const std::string& code_name,
         int iteration_id, const std::string& phase) {
    
    static std::vector<std::string> log_buffer;  // Buffer to hold log messages
    
    // Obtaining current time point and microseconds since epoch
    auto now = std::chrono::system_clock::now();
    auto epoch_time = now.time_since_epoch();
    auto microsecs_since_epoch = std::chrono::duration_cast<std::chrono::microseconds>(epoch_time).count();
    auto now_as_time_t = std::chrono::system_clock::to_time_t(now);
    auto utc_time = gmtime(&now_as_time_t);
    
    // Formatting the log message
    std::ostringstream oss;
    oss << "UTC Time: " << std::put_time(utc_time, "%Y-%m-%d %H:%M:%S")
        << ", " << microsecs_since_epoch << " ms since epoch"
        << ", Code: " << code_name
        << ", Iteration: " << iteration_id
        << ", Phase: " << phase
        << ", Subject: " << message_subject
        << " - " << message << "\n";
    
    // Adding the formatted message to the buffer
    log_buffer.push_back(oss.str());
    
    // Writing logs to file if buffer size threshold is met
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

// flush_remaining_logs function:
// Writes any remaining logs in the buffer to the specified log file.

void flush_remaining_logs(const std::string& log_file_path) {
    static std::vector<std::string> log_buffer;  // Buffer for log messages
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

// buildLogMessage function:
// Constructs a formatted log message with statistics.

std::string buildLogMessage(const std::string& title, const OnlineStats& stats) {
    std::ostringstream logMessage;
    logMessage << title << " Time: "
               << "Mean: " << std::fixed << std::setprecision(2) << stats.getMean()
               << ", Min: " << stats.getMin()
               << ", Max: " << stats.getMax()
               << ", Standard Deviation: " << stats.getStandardDeviation()
               << ", Total Time: " << stats.getTotalSum()
               << ", Count: " << stats.getCount();

    return logMessage.str();
}

// durationToString function:
// Converts a std::chrono::duration to a string with specified precision.

std::string durationToString(const std::chrono::duration<double>& duration, int precisionVal) {
    std::ostringstream oss;
    auto microsecs = std::chrono::duration_cast<std::chrono::microseconds>(duration);
    oss << std::fixed << std::setprecision(precisionVal) << microsecs.count();
    return oss.str();
}


// random_from_urandom function:
// Generates a random number using the /dev/urandom device and ensures it is less than max_value.

void random_from_urandom(mpz_t result, const mpz_t max_value) {
    // Determine the number of bytes needed to represent max_value
    size_t bytes = (mpz_sizeinbase(max_value, 2) + 7) / 8;
    unsigned char *buf = new unsigned char[bytes];
    
    // Open /dev/urandom for reading
    std::ifstream urandom("/dev/urandom", std::ios::in|std::ios::binary);
    if (!urandom) {
        throw std::runtime_error("Failed to open /dev/urandom");
    }
    // Read random bytes from /dev/urandom
    urandom.read(reinterpret_cast<char*>(buf), bytes);
    urandom.close();

    // Convert the random bytes into an mpz_t number
    mpz_import(result, bytes, 1, 1, 0, 0, buf);
    delete[] buf;

    // Ensure the result is within the range [0, max_value)
    mpz_mod(result, result, max_value);
}



// random_from_urandom_bits function:
// Generates a random number of a specified bit length using the /dev/urandom device.

void random_from_urandom_bits(mpz_t result, unsigned int bits) {
    unsigned int size_in_bytes = (bits + 7) / 8; // Calculate bytes from bits
    std::vector<unsigned char> buffer(size_in_bytes);
    std::ifstream urandom("/dev/urandom", std::ios::in | std::ios::binary);

    if (!urandom) {
        throw std::runtime_error("Failed to open /dev/urandom");
    }

    // Read random bytes
    urandom.read(reinterpret_cast<char*>(buffer.data()), size_in_bytes);
    urandom.close();

    // Convert the random bytes into an mpz_t number
    mpz_import(result, size_in_bytes, 1, sizeof(buffer[0]), 0, 0, buffer.data());

    // Ensure the result fits within the specified bit length
    mpz_t two_power_bits;
    mpz_init(two_power_bits);
    mpz_ui_pow_ui(two_power_bits, 2, bits); // 2^bits
    mpz_mod(result, result, two_power_bits);

    mpz_clear(two_power_bits);
}

// random_double_from_urandom function:
// Generates a random double precision floating-point number using /dev/urandom.

double randomDoubleFromUrandom() {
    mpz_t result;
    mpz_init(result);
    const unsigned int bits = 52; // Precision for double
    random_from_urandom_bits(result, bits);

    // Setting up max_value as 2^bits
    mpz_t max_value;
    mpz_init(max_value);
    mpz_ui_pow_ui(max_value, 2, bits);

    // Converting mpz_t values to floating-point
    mpf_t float_result, float_max_value;
    mpf_init(float_result);
    mpf_init(float_max_value);
    mpf_set_z(float_result, result);
    mpf_set_z(float_max_value, max_value);

    // Calculating the double precision floating-point result
    mpf_div(float_result, float_result, float_max_value);
    double secure_random_double = mpf_get_d(float_result);

    // Cleaning up
    mpz_clear(result);
    mpz_clear(max_value);
    mpf_clear(float_result);
    mpf_clear(float_max_value);

    return secure_random_double;
}




// random_index_from_urandom function:
// Generates a random index (as size_t) less than a specified max_value using /dev/urandom.

size_t random_index_from_urandom(size_t max_value) {
    mpz_t result, mpz_max;
    mpz_init(result);
    mpz_init_set_ui(mpz_max, max_value);

    // Generate a random number
    random_from_urandom(result, mpz_max);

    // Convert the mpz_t number to size_t
    size_t index = mpz_get_ui(result);

    // Cleaning up
    mpz_clear(result);
    mpz_clear(mpz_max);

    return index;
}
