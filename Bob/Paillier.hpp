//  paillier.h
//  PaillierCrypto
//  Created by Süleyman Özdel on 6.09.2023.

#ifndef PAILLIER_hpp
#define PAILLIER_hpp

#include <gmp.h>

class Paillier {
public:
    Paillier(unsigned int bit_length);
    ~Paillier();

    void encrypt(mpz_t ciphertext, const mpz_t plaintext) const;
    void decrypt(mpz_t plaintext, const mpz_t ciphertext);
    void subtract(mpz_t result, const mpz_t ciphertext_a, const mpz_t ciphertext_b);
    void add(mpz_t result, const mpz_t ciphertext1, const mpz_t ciphertext2);
    void scalar_multiply(mpz_t result, const mpz_t ciphertext, const mpz_t scalar);
    void scalar_add(mpz_t result, const mpz_t ciphertext, const mpz_t scalar);
    void get_public_key(mpz_t out_n, mpz_t out_g) const;

private:
    void generate_random_prime(mpz_t prime, unsigned int bits);
    void L(mpz_t result, const mpz_t x);

    mpz_t *p;
    mpz_t *q;
    mpz_t *n;
    mpz_t *n_square;
    mpz_t *g;
    mpz_t *lambda;
    mpz_t *lambda_inv;
    __gmp_randstate_struct* state;
    
};

#endif // PAILLIER_H
