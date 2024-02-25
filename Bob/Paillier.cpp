#include "Paillier.hpp"
#include <ctime>
#include <cstdlib>
#include "utils.hpp"


Paillier::Paillier(unsigned int bit_length) {
    // Allocate and initialize mpz_t members using malloc
    p = (mpz_t*)malloc(sizeof(mpz_t));
    q = (mpz_t*)malloc(sizeof(mpz_t));
    n = (mpz_t*)malloc(sizeof(mpz_t));
    n_square = (mpz_t*)malloc(sizeof(mpz_t));
    g = (mpz_t*)malloc(sizeof(mpz_t));
    lambda = (mpz_t*)malloc(sizeof(mpz_t));
    lambda_inv = (mpz_t*)malloc(sizeof(mpz_t));

    mpz_init(*p);
    mpz_init(*q);
    mpz_init(*n);
    mpz_init(*n_square);
    mpz_init(*g);
    mpz_init(*lambda);
    mpz_init(*lambda_inv);

    // Initialize random state
    state = reinterpret_cast<__gmp_randstate_struct*>(malloc(sizeof(__gmp_randstate_struct)));
    gmp_randinit_default(state);
    gmp_randseed_ui(state, std::time(nullptr));

    // Generate p and q
    generate_random_prime(*p, bit_length / 2);
    generate_random_prime(*q, bit_length / 2);

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


    mpz_clear(p1);
    mpz_clear(q1);
}


Paillier::~Paillier() {
    int a=1;
    /*
    mpz_clear(*p);
    mpz_clear(*q);
    mpz_clear(*n);
    mpz_clear(*n_square);
    mpz_clear(*g);
    mpz_clear(*lambda);
    mpz_clear(*lambda_inv);

    free(p);
    free(q);
    free(n);
    free(n_square);
    free(g);
    free(lambda);
    free(lambda_inv);

    gmp_randclear(state);
    free(state);
     */
}

void Paillier::encrypt(mpz_t ciphertext, const mpz_t plaintext) const{
    mpz_t adjusted_plaintext, n_div_2;
    mpz_init(adjusted_plaintext);
    mpz_init(n_div_2);

    // Compute n/2 for comparison
    mpz_divexact_ui(n_div_2, *n, 2);

    // If plaintext is negative and in range [-n/2, n/2)
    if (mpz_cmp(plaintext, n_div_2) < 0) {
        mpz_add(adjusted_plaintext, plaintext, *n);
    } else {
        mpz_set(adjusted_plaintext, plaintext);
    }
    
    mpz_t r, n_minus_1;
    
    mpz_init(r);
    mpz_init(n_minus_1);

    mpz_sub_ui(n_minus_1, *n, 1);
    mpz_urandomm(r, state, n_minus_1);
    mpz_add_ui(r, r, 1);

    mpz_powm_sec(ciphertext, *g, adjusted_plaintext, *n_square);
    mpz_powm_sec(r, r, *n, *n_square);
    mpz_mul(ciphertext, ciphertext, r);
    mpz_mod(ciphertext, ciphertext, *n_square);

    mpz_clear(r);
    mpz_clear(n_minus_1);
    mpz_clear(adjusted_plaintext);
    mpz_clear(n_div_2);
}



void Paillier::decrypt(mpz_t plaintext, const mpz_t ciphertext) {
    mpz_t temp, n_div_2;

    mpz_init(temp);
    mpz_init(n_div_2);

    mpz_powm_sec(temp, ciphertext, *lambda, *n_square);
    L(temp, temp);
    mpz_mul(plaintext, temp, *lambda_inv);
    mpz_mod(plaintext, plaintext, *n);

    mpz_divexact_ui(n_div_2, *n, 2);
    if (mpz_cmp(plaintext, n_div_2) > 0) {
        mpz_sub(plaintext, plaintext, *n);
    }


    mpz_clear(temp);
    mpz_clear(n_div_2);
}

void Paillier::subtract(mpz_t result, const mpz_t ciphertext_a, const mpz_t ciphertext_b) {
    mpz_t inv_b;
    mpz_init(inv_b);

    mpz_invert(inv_b, ciphertext_b, *n_square);
    mpz_mul(result, ciphertext_a, inv_b);
    mpz_mod(result, result, *n_square);

    mpz_clear(inv_b);
}

void Paillier::add(mpz_t result, const mpz_t ciphertext1, const mpz_t ciphertext2) {
    mpz_mul(result, ciphertext1, ciphertext2);
    mpz_mod(result, result, *n_square);
}



void Paillier::scalar_multiply(mpz_t result, const mpz_t ciphertext, const mpz_t scalar) {
    mpz_t adjusted_scalar;
    mpz_init(adjusted_scalar);

    // Adjust the scalar to be within [0, n-1]
    mpz_mod(adjusted_scalar, scalar, *n);

    mpz_powm_sec(result, ciphertext, adjusted_scalar, *n_square);

    mpz_clear(adjusted_scalar);
}

void Paillier::scalar_add(mpz_t result, const mpz_t ciphertext, const mpz_t scalar) {
    mpz_t temp;
    mpz_init(temp);

    mpz_powm_sec(temp, *g, scalar, *n_square);  // g raised to the power of scalar mod n^2
    mpz_mul(result, ciphertext, temp);          // Multiply with the ciphertext
    mpz_mod(result, result, *n_square);         // Take mod n^2 to stay in the encrypted domain

    mpz_clear(temp);
}


void Paillier::get_public_key(mpz_t out_n, mpz_t out_g) const {
    mpz_realloc2(out_n, 8);

    mpz_set(out_n, *n);
    mpz_realloc2(out_g, 8);
    mpz_set(out_g, *g);
}

void Paillier::generate_random_prime(mpz_t prime, unsigned int bits) {
    mpz_urandomb(prime, state, bits);
    mpz_nextprime(prime, prime);
}

void Paillier::L(mpz_t result, const mpz_t x) {
    mpz_sub_ui(result, x, 1);
    mpz_divexact(result, result, *n);
}

/*
//  paillier.cpp
//  PaillierCrypto
//  Created by Süleyman Özdel on 6.09.2023.

#include "paillier.h"
#include <ctime>
#include <cstdlib>

Paillier::Paillier(unsigned int bit_length, mpz_t out_n, mpz_t out_g) {
    // Initialize random state
    gmp_randinit_default(state);
    gmp_randseed_ui(state, std::time(nullptr));

    mpz_inits(p, q, n, n_square, g, lambda, lambda_inv, NULL);

    // Generate p and q
    generate_random_prime(p, bit_length / 2);
    generate_random_prime(q, bit_length / 2);

    mpz_mul(n, p, q);
    mpz_mul(n_square, n, n);

    mpz_t p1, q1;
    mpz_inits(p1, q1, NULL);
    mpz_sub_ui(p1, p, 1);
    mpz_sub_ui(q1, q, 1);
    mpz_lcm(lambda, p1, q1);

    mpz_add_ui(g, n, 1);

    mpz_invert(lambda_inv, lambda, n);

    get_public_key(out_n, out_g);
    mpz_clears(p1, q1, NULL);
}

Paillier::~Paillier() {
    mpz_clears(p, q, n, n_square, g, lambda, lambda_inv, NULL);
    gmp_randclear(state);
}



void Paillier::encrypt(mpz_t ciphertext, const mpz_t plaintext) {
    mpz_t r, n_minus_1;
    mpz_inits(r, n_minus_1, NULL);

    mpz_sub_ui(n_minus_1, n, 1);
    mpz_urandomm(r, state, n_minus_1);
    mpz_add_ui(r, r, 1);  // Ensure r is in [1, n-1]

    mpz_powm_sec(ciphertext, g, plaintext, n_square);
    mpz_powm_sec(r, r, n, n_square);
    mpz_mul(ciphertext, ciphertext, r);
    mpz_mod(ciphertext, ciphertext, n_square);

    mpz_clears(r, n_minus_1, NULL);
}

void Paillier::decrypt(mpz_t plaintext, const mpz_t ciphertext) {
    mpz_t temp, n_div_2;
    mpz_inits(temp, n_div_2, NULL);

    mpz_powm_sec(temp, ciphertext, lambda, n_square);
    L(temp, temp);
    mpz_mul(plaintext, temp, lambda_inv);  // Use the precomputed inverse here
    mpz_mod(plaintext, plaintext, n);

    // Adjust for negative numbers
    mpz_divexact_ui(n_div_2, n, 2);
    if (mpz_cmp(plaintext, n_div_2) > 0) {
        mpz_sub(plaintext, plaintext, n);
    }

    mpz_clears(temp, n_div_2, NULL);
}


void Paillier::subtract(mpz_t result, const mpz_t ciphertext_a, const mpz_t ciphertext_b) {
    mpz_t inv_b;
    mpz_init(inv_b);
    

    // Compute the modular inverse of ciphertext_b mod n^2
    mpz_invert(inv_b, ciphertext_b, n_square);

    // Multiply ciphertext_a with the inverse of ciphertext_b
    mpz_mul(result, ciphertext_a, inv_b);
    mpz_mod(result, result, n_square);

    mpz_clear(inv_b);
}

void Paillier::add(mpz_t result, const mpz_t ciphertext1, const mpz_t ciphertext2) {
    mpz_mul(result, ciphertext1, ciphertext2);
    mpz_mod(result, result, n_square);
}

void Paillier::scalar_multiply(mpz_t result, const mpz_t ciphertext, const mpz_t scalar) {
    mpz_t adjusted_scalar;
    mpz_init(adjusted_scalar);

    // Adjust the scalar to be within [0, n-1]
    mpz_mod(adjusted_scalar, scalar, n);

    mpz_powm_sec(result, ciphertext, adjusted_scalar, n_square);

    mpz_clear(adjusted_scalar);
}

void Paillier::scalar_add(mpz_t result, const mpz_t ciphertext, const mpz_t scalar) {
    mpz_t temp;
    mpz_init(temp);

    mpz_powm_sec(temp, g, scalar, n_square);  // g raised to the power of scalar mod n^2
    mpz_mul(result, ciphertext, temp);       // Multiply with the ciphertext
    mpz_mod(result, result, n_square);       // Take mod n^2 to stay in the encrypted domain

    mpz_clear(temp);
}

void Paillier::get_public_key(mpz_t out_n, mpz_t out_g) const {
    mpz_set(out_n, n);
    mpz_set(out_g, g);
}




void Paillier::generate_random_prime(mpz_t prime, unsigned int bits) {
    mpz_urandomb(prime, state, bits);
    mpz_nextprime(prime, prime);
}

void Paillier::L(mpz_t result, const mpz_t x) {
    mpz_sub_ui(result, x, 1);
    mpz_divexact(result, result, n);
}
*/
