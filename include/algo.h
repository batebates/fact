#include "gmp.h"
#ifndef ALGO_H_INCLUDED
#define ALGO_H_INCLUDED

void FirstPrimeDiviseur(mpz_t r,mpz_t n,int parameter);

void DiviseurSucc(mpz_t r,mpz_t n,int parameter);

void BrentFactor(mpz_t r,mpz_t n,int parameter);

void Pollardrho(mpz_t r,mpz_t n,int parameter);

void Pollardp(mpz_t r,mpz_t n,int parameter);

#endif
