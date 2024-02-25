//  Paillier.h
//  PaillierCrypto
//  Created by Süleyman Özdel on 6.09.2023.

#ifndef PAILLIER_hpp
#define PAILLIER_hpp

#include <gmp.h>

class Paillier {
public:
    // Constructor: Initializes the Paillier cryptosystem with a specified key size.
    Paillier(unsigned int bit_length);

    // Destructor: Cleans up and deallocates resources.
    ~Paillier();

    // Encrypts a plaintext message using the Paillier cryptosystem.
    void encryptMessage(mpz_t ciphertext, const mpz_t plaintext) const;

    // Decrypts a ciphertext message using the Paillier cryptosystem.
    void decryptMessage(mpz_t plaintext, const mpz_t ciphertext);

    // Performs subtraction on two ciphertexts in the encrypted domain.
    void subtractCiphertexts(mpz_t result, const mpz_t ciphertext_a, const mpz_t ciphertext_b);

    // Performs addition on two ciphertexts in the encrypted domain.
    void addCiphertexts(mpz_t result, const mpz_t ciphertext1, const mpz_t ciphertext2);

    // Multiplies a ciphertext by a scalar in the encrypted domain.
    void multiplyByScalar(mpz_t result, const mpz_t ciphertext, const mpz_t scalar);

    // Adds a scalar to a ciphertext in the encrypted domain.
    void addScalar(mpz_t result, const mpz_t ciphertext, const mpz_t scalar);

    // Gets the public key components (n and g) of the Paillier cryptosystem.
    void getPublicKey(mpz_t out_n, mpz_t out_g) const;

    // Updates the key materials of the Paillier cryptosystem.
    void updateEncryptionKey(unsigned int bit_length);


private:
    // Generates a random prime number of a specified bit length.
    void generateRandomPrime(mpz_t prime, unsigned int bits);

    // Helper function L used in the decryption process.
    void applyLFunction(mpz_t result, const mpz_t x);

    // Precomputes values necessary for encryption and decryption.
    void precomputeKeyValues();

    // Initializes the keys for the Paillier cryptosystem.
    void initializeEncryptionKeys(unsigned int bit_length);

    // Clears and deallocates all key materials.
    void clearVariables();

    // Key components and precomputed values.
    mpz_t *p;
    mpz_t *q;
    mpz_t *n;
    mpz_t *n_square;
    mpz_t *g;
    mpz_t *lambda;
    mpz_t *lambda_inv;
    mpz_t *n_div_2; // Pre-computed n/2
};

#endif // PAILLIER_H
