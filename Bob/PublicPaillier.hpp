#ifndef PUBLICPAILLIER_HPP
#define PUBLICPAILLIER_HPP
#include <memory>
#include <gmp.h>

class PublicPaillier {
public:
    mpz_t n_square;

    PublicPaillier(const mpz_t public_n, const mpz_t public_g);
    ~PublicPaillier();

    void encrypt(mpz_t ciphertext, const mpz_t plaintext) const;
    void add(mpz_t result, const mpz_t ciphertext1, const mpz_t ciphertext2) const;
    void scalarMultiply(mpz_t result, const mpz_t ciphertext, const mpz_t scalar) const;
    void scalarAdd(mpz_t result, const mpz_t ciphertext, const mpz_t scalar) const;
    void randomWithInverse(mpz_t result, mpz_t inverse);
    void randomWithInverseMultiplication(mpz_t result, mpz_t inverse);
    void randomWithShiftedInverse(mpz_t result,mpz_t inverse);
    void randomWithAdditiveInverse(mpz_t result, mpz_t additive_inverse);
    //void random_with_additive_inverse(mpz_t result);

    void updatePublicParams(mpz_t& new_public_n, mpz_t& new_public_g);

private:
    mpz_t n_half;
    mpz_t n_minus_1;
    mpz_t n;
    mpz_t g;
    
    //__gmp_randstate_struct* state;
};

#endif // PUBLICPAILLIER_H
