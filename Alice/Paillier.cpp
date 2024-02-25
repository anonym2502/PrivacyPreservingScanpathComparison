#include "Paillier.hpp"
#include <ctime>
#include <cstdlib>
#include "utils.hpp"
#include <fstream>

// Constructor: Initializes the Paillier cryptosystem keys.
Paillier::Paillier(unsigned int bit_length) {
    initializeEncryptionKeys(bit_length); // Initialize keys with given bit length.
    precomputeKeyValues(); // Precompute additional values for the cryptosystem.
}

// Destructor: Cleans up and frees allocated resources.
Paillier::~Paillier() {
    clearVariables(); // Clears and deallocates all mpz_t members.
}

// Encrypts plaintext using Paillier encryption.
void Paillier::encryptMessage(mpz_t ciphertext, const mpz_t plaintext) const {
    // Initialize variables for adjusted plaintext and r^n.
    mpz_t adjusted_plaintext, r_to_n;
    mpz_init(adjusted_plaintext);
    mpz_init(r_to_n);

    // Adjust plaintext for negative values within range [-n/2, n/2).
    if (mpz_cmp(plaintext, *n_div_2) < 0) {
        mpz_add(adjusted_plaintext, plaintext, *n);
    } else {
        mpz_set(adjusted_plaintext, plaintext);
    }
    
    // Initialize variables for random number r and gcd result.
    mpz_t r, gcd_result;
    mpz_init(r);
    mpz_init(gcd_result);

    // Generate random r co-prime with n.
    bool isRValid;
    do {
        isRValid = true;
        do {
            random_from_urandom(r, *n);
            mpz_gcd(gcd_result, r, *n);
        } while (mpz_cmp_ui(gcd_result, 1) != 0);

        // Ensure r is not zero.
        if (mpz_cmp_ui(r, 0) == 0) {
            isRValid = false;
        }
    } while (!isRValid);

    // Calculate ciphertext: (nm + 1) * r^n mod n^2.
    mpz_mul(ciphertext, adjusted_plaintext, *n);
    mpz_add_ui(ciphertext, ciphertext, 1);
    mpz_powm_sec(r_to_n, r, *n, *n_square);
    mpz_mul(ciphertext, ciphertext, r_to_n);
    mpz_mod(ciphertext, ciphertext, *n_square);

    // Clear all temporary variables.
    mpz_clear(r);
    mpz_clear(r_to_n);
    mpz_clear(adjusted_plaintext);
    mpz_clear(gcd_result);
}

// Decrypts ciphertext using Paillier decryption.
void Paillier::decryptMessage(mpz_t plaintext, const mpz_t ciphertext) {
    mpz_t temp;
    mpz_init(temp);

    // Compute L(c^lambda mod n^2) * lambda_inv mod n.
    mpz_powm_sec(temp, ciphertext, *lambda, *n_square);
    applyLFunction(temp, temp);
    mpz_mul(plaintext, temp, *lambda_inv);
    mpz_mod(plaintext, plaintext, *n);

    // Adjust for values greater than n/2.
    mpz_divexact_ui(*n_div_2, *n, 2);
    if (mpz_cmp(plaintext, *n_div_2) > 0) {
        mpz_sub(plaintext, plaintext, *n);
    }

    // Clear temporary variable.
    mpz_clear(temp);
}

// Subtracts two ciphertexts in the encrypted domain.
void Paillier::subtractCiphertexts(mpz_t result, const mpz_t ciphertext_a, const mpz_t ciphertext_b) {
    mpz_t inv_b;
    mpz_init(inv_b);

    // Compute inverse of ciphertext_b and multiply with ciphertext_a.
    mpz_invert(inv_b, ciphertext_b, *n_square);
    mpz_mul(result, ciphertext_a, inv_b);
    mpz_mod(result, result, *n_square);

    // Clear temporary variable.
    mpz_clear(inv_b);
}

// Adds two ciphertexts in the encrypted domain.
void Paillier::addCiphertexts(mpz_t result, const mpz_t ciphertext1, const mpz_t ciphertext2) {
    // Multiply ciphertexts and take modulo n^2.
    mpz_mul(result, ciphertext1, ciphertext2);
    mpz_mod(result, result, *n_square);
}

// Multiplies a ciphertext by a scalar in the encrypted domain.
void Paillier::multiplyByScalar(mpz_t result, const mpz_t ciphertext, const mpz_t scalar) {
    mpz_t adjusted_scalar;
    mpz_init(adjusted_scalar);

    // Adjust scalar to be within range [0, n-1].
    mpz_mod(adjusted_scalar, scalar, *n);

    // Raise ciphertext to the power of adjusted scalar mod n^2.
    mpz_powm_sec(result, ciphertext, adjusted_scalar, *n_square);

    // Clear temporary variable.
    mpz_clear(adjusted_scalar);
}

// Adds a scalar to a ciphertext in the encrypted domain.
void Paillier::addScalar(mpz_t result, const mpz_t ciphertext, const mpz_t scalar) {
    mpz_t temp;
    mpz_init(temp);

    // Raise g to the power of scalar mod n^2 and multiply with ciphertext.
    mpz_powm_sec(temp, *g, scalar, *n_square);
    mpz_mul(result, ciphertext, temp);
    mpz_mod(result, result, *n_square);

    // Clear temporary variable.
    mpz_clear(temp);
}

// Retrieves public key components n and g.
void Paillier::getPublicKey(mpz_t out_n, mpz_t out_g) const {
    // Assign public key components to output variables.
    mpz_realloc2(out_n, 8);
    mpz_set(out_n, *n);
    mpz_realloc2(out_g, 8);
    mpz_set(out_g, *g);
}

// Generates a random prime number of specified bit length.
void Paillier::generateRandomPrime(mpz_t prime, unsigned int bits) {
    // Generate random number and find the next prime number.
    random_from_urandom_bits(prime, bits);
    mpz_nextprime(prime, prime);
}

// Helper function L used in the decryption process.
void Paillier::applyLFunction(mpz_t result, const mpz_t x) {
    // L function: (x - 1) / n.
    mpz_sub_ui(result, x, 1);
    mpz_divexact(result, result, *n);
}

// Updates the key materials for the Paillier cryptosystem.
void Paillier::updateEncryptionKey(unsigned int bit_length) {
    // Clear and reallocate key materials.
    clearVariables();
    initializeEncryptionKeys(bit_length);
    precomputeKeyValues();
}

// Precomputes values needed for the cryptosystem.
void Paillier::precomputeKeyValues() {
    // Compute and store n divided by 2.
    mpz_init(*n_div_2);
    mpz_divexact_ui(*n_div_2, *n, 2);
}


void Paillier::initializeEncryptionKeys(unsigned int bit_length) {
    
    // Allocate and initialize mpz_t members using malloc
    p = (mpz_t*)malloc(sizeof(mpz_t));
    q = (mpz_t*)malloc(sizeof(mpz_t));
    n = (mpz_t*)malloc(sizeof(mpz_t));
    n_square = (mpz_t*)malloc(sizeof(mpz_t));
    g = (mpz_t*)malloc(sizeof(mpz_t));
    lambda = (mpz_t*)malloc(sizeof(mpz_t));
    lambda_inv = (mpz_t*)malloc(sizeof(mpz_t));
    n_div_2 = (mpz_t*)malloc(sizeof(mpz_t));
    
    mpz_init(*p);
    mpz_init(*q);
    mpz_init(*n);
    mpz_init(*n_square);
    mpz_init(*g);
    mpz_init(*lambda);
    mpz_init(*lambda_inv);
    
    // Generate two random primes, p and q.
    generateRandomPrime(*p, bit_length / 2);
    generateRandomPrime(*q, bit_length / 2);
    
    // Compute n, n^2, g, and other necessary values.
    mpz_mul(*n, *p, *q);
    mpz_mul(*n_square, *n, *n);
    
    mpz_t p1, q1;
    mpz_init(p1);
    mpz_init(q1);
    mpz_sub_ui(p1, *p, 1);
    mpz_sub_ui(q1, *q, 1);
    mpz_lcm(*lambda, p1, q1);

    mpz_add_ui(*g, *n, 1);
    mpz_invert(*lambda_inv, *lambda, *n);

    precomputeKeyValues();

    // Clear temporary variables.
    mpz_clear(p1);
    mpz_clear(q1);
}

void Paillier::clearVariables(){
    // Clear and free all mpz_t members.
    mpz_clears(*p, *q, *n, *n_square, *g, *lambda, *lambda_inv,*n_div_2, NULL);
    
    free(p);
    free(q);
    free(n);
    free(n_square);
    free(g);
    free(lambda);
    free(lambda_inv);
    free(n_div_2);
}






