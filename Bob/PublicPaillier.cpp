#include "PublicPaillier.hpp"
#include <ctime>
#include <cstdlib>
#include <stdexcept>
#include <memory>
#include <fstream>
#include "utils.hpp"

// Constructor for PublicPaillier: initializes public parameters n and g, and their dependent values
PublicPaillier::PublicPaillier(const mpz_t public_n, const mpz_t public_g)
{
    mpz_init_set(n, public_n);     // Initialize and set n
    mpz_init_set(g, public_g);     // Initialize and set g
    mpz_init(n_square);            // Initialize n_square
    mpz_init(n_half);              // Initialize n_half
    mpz_init(n_minus_1);           // Initialize n_minus_1

    mpz_mul(n_square, n, n);       // Compute n^2
    mpz_tdiv_q_ui(n_half, n, 2);   // Compute n/2
    mpz_sub_ui(n_minus_1, n, 1);   // Compute n-1

    // Random number generation setup (commented out)
    //gmp_randinit_default(state.get());
    //gmp_randseed_ui(state.get(), std::time(nullptr));
}

// Destructor for PublicPaillier: clears all mpz_t variables to free memory
PublicPaillier::~PublicPaillier() {
    mpz_clears(n, n_half, n_minus_1, n_square, g, NULL);
    // Random number state cleanup (commented out)
    //gmp_randclear(state);
    //delete state;
    //free(state);
}

// encrypt function: Encrypts plaintext using Paillier encryption scheme
void PublicPaillier::encrypt(mpz_t ciphertext, const mpz_t plaintext) const {
    mpz_t r;  // Random number for encryption
    mpz_inits(r, NULL);

    // Ensure plaintext is in range [-n/2, n/2]
    if (mpz_cmp(plaintext, n_half) > 0 || mpz_cmp_si(plaintext, -mpz_get_si(n_half)) < 0) {
        mpz_clears(r, NULL);
        throw std::runtime_error("Plaintext out of range");
    }

    // Random number generation for encryption
    mpz_t gcd_result;  // For checking GCD with n
    mpz_init(gcd_result);
    do {
        random_from_urandom(r, n);  // Generate a random number less than n
        mpz_gcd(gcd_result, r, n);  // Compute GCD of r and n
    } while (mpz_cmp_ui(gcd_result, 1) != 0);  // Ensure GCD is 1
    
    mpz_add_ui(r, r, 1);  // Adjust r to be in [1, n-1]

    // Paillier encryption computations
    mpz_mul(ciphertext, plaintext, n);  // nm
    mpz_add_ui(ciphertext, ciphertext, 1);  // nm + 1
    mpz_powm_sec(r, r, n, n_square);  // r^n mod n^2
    mpz_mul(ciphertext, ciphertext, r);
    mpz_mod(ciphertext, ciphertext, n_square);

    mpz_clears(r, NULL);
}


// Performs homomorphic addition of two ciphertexts in the Paillier cryptosystem.
void PublicPaillier::add(mpz_t result, const mpz_t ciphertext1, const mpz_t ciphertext2) const {
    // Multiply the two ciphertexts together
    // Since the Paillier cryptosystem is homomorphic with respect to multiplication,
    // this operation effectively adds the corresponding plaintexts.
    mpz_mul(result, ciphertext1, ciphertext2);
    mpz_mod(result, result, n_square);
}


// Performs homomorphic scalar addition of a ciphertext with a scalar value.
void PublicPaillier::scalarAdd(mpz_t result, const mpz_t ciphertext, const mpz_t scalar) const {
    mpz_t g_to_scalar;  // Variable to store g raised to the power of scalar
    mpz_init(g_to_scalar);

    // Compute g^scalar mod n^2. This is part of the Paillier encryption of the scalar.
    // In the Paillier cryptosystem, encrypting a plaintext value 'x' is done by computing
    // g^x mod n^2. Here, 'scalar' is treated as a plaintext value.
    mpz_powm_sec(g_to_scalar, g, scalar, n_square);

    // Multiply the original ciphertext by the computed g^scalar.
    // This is the homomorphic equivalent of adding the scalar value to the encrypted plaintext.
    // The operation is performed modulo n^2 to maintain consistency with the Paillier scheme.
    mpz_mul(result, ciphertext, g_to_scalar);
    mpz_mod(result, result, n_square);

    // Clear the temporary variable used for g^scalar
    mpz_clear(g_to_scalar);
}


// Performs homomorphic scalar multiplication of a ciphertext by a scalar value.
void PublicPaillier::scalarMultiply(mpz_t result, const mpz_t ciphertext, const mpz_t scalar) const {
    mpz_t adjusted_scalar;  // Variable to store the adjusted scalar value
    mpz_init(adjusted_scalar);

    // Adjust the scalar to be within the range [0, n-1].
    // This is necessary to ensure that the scalar multiplication stays within
    // the bounds of the Paillier cryptosystem's group.
    mpz_mod(adjusted_scalar, scalar, n);

    // Perform the scalar multiplication in a secure manner using the 'mpz_powm_sec' function.
    // This function raises 'ciphertext' to the power of 'adjusted_scalar' modulo n^2.
    // In the Paillier cryptosystem, this operation corresponds to multiplying the plaintext
    // by the scalar before encryption.
    mpz_powm_sec(result, ciphertext, adjusted_scalar, n_square);

    // Clear the temporary variable used for adjusted scalar
    mpz_clear(adjusted_scalar);
}

    
// Generates a random number with its multiplicative inverse modulo n.
void PublicPaillier::randomWithInverse(mpz_t result, mpz_t inverse) {
    mpz_t n_limit;
    mpz_init(n_limit);

    // Set a limit for the random number (It can be maximized according to requirements and Paillier scheme's defined range)
    mpz_tdiv_q_ui(n_limit, n_half, 2000000000);
    do {
        // Generate a random number less than the limit
        random_from_urandom(result, n_limit);
        // Ensure result is in the range [1, n-1]
        mpz_add_ui(result, result, 1);
    } while (mpz_invert(inverse, result, n) == 0);  // Repeat until a number with a multiplicative inverse modulo n is found

    mpz_clear(n_limit);
}
// Generates a random number with its multiplicative inverse modulo n, with a specific limit.
void PublicPaillier::randomWithInverseMultiplication(mpz_t result, mpz_t inverse) {
    mpz_t n_limit;
    mpz_init(n_limit);

    // Set a limit for the random number (It can be maximized according to requirements and Paillier scheme's defined range)
    mpz_set_str(n_limit, "50000", 10);
    do {
        // Generate a random number less than the limit
        random_from_urandom(result, n_limit);
        // Ensure result is in the range [1, n-1]
        mpz_add_ui(result, result, 1);
    } while (mpz_invert(inverse, result, n) == 0);  // Repeat until a number with a multiplicative inverse modulo n is found

    mpz_clear(n_limit);
}
// Generates a random number with its additive inverse modulo n.
void PublicPaillier::randomWithAdditiveInverse(mpz_t result, mpz_t additive_inverse) {
    mpz_t n_limit;
    mpz_init(n_limit);

    // Set a limit for the random number (It can be maximized according to requirements and Paillier scheme's defined range)
    mpz_tdiv_q_ui(n_limit, n_half, 2000000000);

    do {
        // Generate a random number less than the limit
        random_from_urandom(result, n_limit);
        // Ensure result is in the range [1, n-1]
        mpz_add_ui(result, result, 1);
    } while (mpz_cmp_ui(result, 0) == 0);  // Repeat until a non-zero number is found

    // Compute the additive inverse of result modulo n
    mpz_sub(additive_inverse, n, result);

    mpz_clear(n_limit);
}
// Generates a random number and its shifted multiplicative inverse modulo n.
void PublicPaillier::randomWithShiftedInverse(mpz_t result, mpz_t shifted_inverse) {
    mpz_t shifted, n_limit;
    mpz_init(shifted);
    mpz_init(n_limit);

    // Set a limit for the random number (It can be maximized according to requirements and Paillier scheme's defined range)
    mpz_set_str(n_limit, "50000", 10);

    do {
        // Generate a random number less than the limit
        random_from_urandom(result, n_limit);
        // Shifted value
        mpz_add_ui(shifted, result, 1);
    } while (mpz_invert(shifted_inverse, shifted, n) == 0);  // Repeat until a shifted number with a multiplicative inverse modulo n is found

    mpz_clears(shifted, n_limit, NULL);
}





void PublicPaillier::updatePublicParams(mpz_t& new_public_n, mpz_t& new_public_g) {
    // Clear the old values
    mpz_clears(n, n_half, n_minus_1, n_square, g, NULL);
    // Initialize and set the new values
    mpz_init_set(n, new_public_n);
    mpz_init_set(g, new_public_g);
    mpz_init(n_square);
    mpz_init(n_half);
    mpz_init(n_minus_1);

    mpz_mul(n_square, n, n);
    mpz_tdiv_q_ui(n_half, n, 2);
    mpz_sub_ui(n_minus_1, n, 1);

    // Reinitialize the random state with a new seed
    //gmp_randseed_ui(state.get(), std::time(nullptr));
}
