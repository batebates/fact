/**
@file algo.c
@brief Algorithm functions that determine a factor of a number
@author Baptiste Decrand, Olivier ADJONYO, Dai Chi Do
@version 1.0
@date 30/01/2019
*/
#include "../include/algo.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "gmp.h"



/**
  This function research a factor of n in a list of prime number
@param r this is the back value that contains a factor of n
@param n Number to factorize.
*/
void FirstPrimeDiviseur(mpz_t r,mpz_t n,int parameter){
  printf("bal1se\n");
  fflush(stdout);
  char ch[2048];
  int iq;
  mpz_t q;
  mpz_init(q);
  FILE *fp;
  fp = fopen("primes50.txt", "r");
  if (fp == NULL)
  {
      perror("Error while opening the file.\n");
      exit(EXIT_FAILURE);
  }
  while (fscanf(fp, " %1023s", ch) == 1) {
        iq = atoi(ch);
        printf("bal3se\n");
        fflush(stdout);
        mpz_set_ui(q,iq);
        printf("bal4se\n");
        fflush(stdout);
        if(mpz_divisible_p(n,q)){
          printf("bal5se\n");
          fflush(stdout);
          mpz_set(r,q);
          break;
        }
  }
  mpz_clear(q);
  printf("bal2se\n");
  fflush(stdout);
  fclose(fp);
}
/**
 function that incrément r parameter
@param r return value of r+1
@param n Le nombre à factoriser.
*/
void DiviseurSucc(mpz_t r,mpz_t n,int parameter){
    if (mpz_cmp_ui(r,parameter)>0) {
      mpz_set_ui(r,0);
    }else{
      mpz_add_ui(r,r,1);
    }

}
/**
 Algorithme de Brent, la fonction recherche un facteur de n en utilisant cet algorithme
 @param r this is the return value that contains a factor of n or 0 if it didn't find it
 @param n Number to factorize.
*/
void BrentFactor(mpz_t r,mpz_t n,int rhomax){
  mpz_t a,b,t;
  long i=0;
  int flag =0;
  mpz_init(a);
  mpz_init(b);
  mpz_init(t);
  mpz_set_ui(r,2);
  mpz_set_ui(a,2);
  mpz_set_ui(b,2);

  while(flag==0  && i<rhomax){
    i++;
    mpz_pow_ui(a,a,2);
    mpz_add_ui(a,a,1);
    mpz_mod(a,a,n);
    mpz_sub(t,a,b);
    mpz_abs(t,t);
    mpz_gcd(r,t,n);
    if(mpz_cmp_ui(r,1)>0 && mpz_cmp(r,n)<0){
      flag=1;
    }
    if(i && (!(i&(i-1)))){
      mpz_set(b,a);
    }
  }

  if(flag==0 || mpz_cmp(r,n)==0){
    mpz_set_ui(r,0);
  }
  mpz_clear(a);
  mpz_clear(b);
  mpz_clear(t);
}
/**
 Algorithme de Pollard rho, la fonction recherche un facteur de n en utilisant cet algorithme
 @param r this is the return value that contains a factor of n or 0 if it didn't find it
 @param n Number to factorize.
*/
void Pollardrho(mpz_t r,mpz_t n,int rhomax){
  mpz_t a,b,t;
  long i=0;
  int flag =0;
  mpz_init(a);
  mpz_init(b);
  mpz_init(t);
  mpz_set_ui(r,2);
  mpz_set_ui(a,2);
  mpz_set_ui(b,2);

  while(flag==0 && mpz_cmp(r,n)!=0 && i<rhomax){
    i++;
    /* f(x) = x^2+1)%N*/
    /* a = f(a)*/
    mpz_pow_ui(a,a,2);
    mpz_add_ui(a,a,1);
    mpz_mod(a,a,n);

    /* b = f(f(b))*/
    mpz_pow_ui(b,b,2);
    mpz_add_ui(b,b,1);
    mpz_mod(b,b,n);

    mpz_pow_ui(b,b,2);
    mpz_add_ui(b,b,1);
    mpz_mod(b,b,n);


    /*r = gcd((a-b),n)*/
    mpz_sub(t,a,b);
    mpz_gcd(r,t,n);

    if(mpz_cmp_ui(r,1)>0){
      flag=1;
    }
  }

  if(flag==0 || mpz_cmp(r,n)==0){
    mpz_set_ui(r,0);
  }
  /*Clear*/
  mpz_clear(a);
  mpz_clear(b);
  mpz_clear(t);
}

/**
 Algorithme de Pollard p-1, la fonction recherche un facteur de n en utilisant cet algorithme
 @param r this is the return value that contains a factor of n or 0 if it didn't find it
 @param n Number to factorize.
*/
void Pollardp(mpz_t r,mpz_t n,int parameter){
  /*Initiation random*/
  mpz_t a;
  gmp_randstate_t rstate;

  unsigned long seed = rand();
  gmp_randinit_default(rstate);
  gmp_randseed_ui(rstate, seed);
  mpz_init(a);
  mpz_urandomm(a,rstate,n);

  /*Initiation*/
  mpz_t q,g;
  mpz_init(q);
  mpz_init(g);
  mpz_set_ui(q,1);

  int iq,e,flag,c;
  flag=0;
  /*Selection du seuil de friabilité de B prime number between 0 and ..*/
  char ch[2048];
  FILE *fp;
  fp = fopen("primes1.txt", "r");
  if (fp == NULL)
  {
      perror("Error while opening the file.\n");
      exit(EXIT_FAILURE);
  }
  c=0;
  while (fscanf(fp, " %1023s", ch) == 1 ) {
        c++;
        iq = atoi(ch);
      	mpz_set_ui(q,iq);
      	e = 6/mpz_sizeinbase(q,10);
        mpz_pow_ui(q,q,e);
        mpz_powm(a,a,q,n);
        mpz_sub_ui (a, a,1);
        mpz_gcd(g,a,n);
        if(mpz_cmp_ui(g,1)!=0 && mpz_cmp(g,n)!=0){
          mpz_set(r,g);
          flag = 1;
          break;
        }
  }
  if(flag==0 ){
    mpz_set_ui(r,0);
  }
  fclose(fp);
  mpz_clear(q);
  mpz_clear(a);
  mpz_clear(g);
}
