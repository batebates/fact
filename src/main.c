/**
@file main.c
@brief Main programm that manage user options, define and use factorisation functions
@author Baptiste Decrand, Olivier ADJONYO, Dai Chi Do
@version 1.0
@date 30/01/2019
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "gmp.h"
#include "../include/algo.h"
#include <time.h>
#include <getopt.h>
#define TEST_MAX 10000
#define SIZE_NUMBER 20
static int verbose_flag;

/**
 Fonction factorisant un nombre en utilisant un des algorithmes choisis
@param n Le nombre à factoriser
@param choix Choix de l'algorithme à utiliser 0: Pollard rho, 1:Division successives, 2: Brent, 3: Pollard p-1
@param parameter Iteration number or band larger
*/
int Factorisation(mpz_t n,int choix,int parameter){
  int retour =1;
  void (*NextFactor[])(mpz_t,mpz_t,int) = {Pollardrho,DiviseurSucc,BrentFactor,Pollardp};
  mpz_t z,p,t;
  int e,i,flag;
  mpz_init(z);
  mpz_init(t);
  mpz_init(p);
  mpz_sqrt(z,n);
  mpz_set_ui(p,2);
  mpz_set(t,n);
  flag=0;
  while(mpz_cmp(p,z)<0){
    e=0;
    while(mpz_divisible_p(t,p)){
      mpz_divexact(t,t,p);
      e++;
    }
    if(e>0){
      for(i=0;i<e;i++){
        gmp_printf("%Zux",p);
        fflush(stdout);
      }
    }
    mpz_sqrt(z,t);
    NextFactor[choix](p,t,parameter);
    if(mpz_cmp_ui(p,0)==0){
      flag=1;
      break;
    }
  }
  if(flag==1){
      switch (mpz_probab_prime_p(t,20)) {
        case 1:
          retour = 10;
          gmp_printf("%Zu\n",t);
          break;
        case 0:
          retour = 100;
          gmp_printf("%Zu\n",t);
          break;
        default:
          gmp_printf("%Zu\n",t);
          if(mpz_sizeinbase (p, 10)>=14){
            fflush(stdout);
          }
          break;
      }
  }else{
    gmp_printf("%Zu\n",t);
  }
  mpz_clear(z);
  mpz_clear(t);
  mpz_clear(p);

  return retour;
}

/**
 Fonction factorisant un nombre en utilisant les algorithmes en fonction de la taille des facteurs
@param n Le nombre à factoriser
*/
int FactorisationOpti(mpz_t n,int parameter){
  int retour =1;
  void (*NextFactor[])(mpz_t,mpz_t,int) = {Pollardrho,DiviseurSucc,BrentFactor,Pollardp,FirstPrimeDiviseur};
  mpz_t z,p,t;
  size_t size;
  int e,i,flag;
  mpz_init(z);
  mpz_init(t);
  mpz_init(p);
  mpz_sqrt(z,n);
  mpz_set_ui(p,2);
  mpz_set(t,n);
  flag=0;
  char ch[2048];
  int iq;
  FILE *fp;
  fp = fopen("primes1.txt", "r");
  if (fp == NULL)
  {
      perror("Error while opening the file.\n");
      exit(EXIT_FAILURE);
  }
  while (fscanf(fp, " %1023s", ch) == 1) {
        iq = atoi(ch);
        mpz_set_ui(p,iq);
        e=0;
        while(mpz_divisible_p(t,p)){
          mpz_divexact(t,t,p);
          e++;
        }
        if(e>0){
          for(i=0;i<e;i++){
            gmp_printf("%Zux",p);
            fflush(stdout);
          }
        }
  }
  fclose(fp);
  while(mpz_cmp(p,z)<0){
    e=0;
    while(mpz_divisible_p(t,p)){
      mpz_divexact(t,t,p);
      e++;
    }
    if(e>0){
      for(i=0;i<e;i++){
        gmp_printf("%Zux",p);
        fflush(stdout);
      }
    }

    mpz_sqrt(z,t);
    size = mpz_sizeinbase (t, 10);
    if(size<=14){

      NextFactor[1](p,t,parameter);
    }else{
      //printf("oui");
      i=0;
      NextFactor[2](p,t,parameter);
      while (mpz_cmp_ui(p,0)==0 && i<1) {

        NextFactor[3](p,t,parameter);
        i++;
      }
    }

    if(mpz_cmp_ui(p,0)==0){
      flag=1;
      break;
    }


  }
  if(mpz_sizeinbase(t,10)<15 && flag==1){
    Factorisation(t,1,10000000);
  }else if(flag==1){

      switch (mpz_probab_prime_p(t,20)) {
        case 1:
          retour = 10;
          gmp_printf("%Zu\n",t);
          break;
        case 0:
          retour = 100;
          gmp_printf("%Zu\n",t);

          break;
        default:
          gmp_printf("%Zu\n",t);
          break;
      }
  }else{
    gmp_printf("%Zu\n",t);
  }

  mpz_clear(z);
  mpz_clear(t);
  mpz_clear(p);
  return retour;
}


/**
 Fonction principal qui va dans un premier temps gérer les options, puis appeler les fonctions de Factorisation adéquates
@param argc nombre d'arguments
@param argv Tableau des arguments
*/
int main(int argc, char* argv[]){
  /*Initiation des variables Aléatoires*/
  srand(time(NULL));
  mpz_t random;
  gmp_randstate_t rstate;
  unsigned long seed = rand();
  gmp_randinit_default(rstate);
  gmp_randseed_ui(rstate, seed);
  mpz_init(random);

  /*Gestions des Options*/
  int c,choix;
  choix = -1;
  int hFlag,fFlag,iFlag,parameter;
  iFlag = 0;
  hFlag = 0;
  fFlag = 0;
  char filename[50];
  while (1)
    {
      static struct option long_options[] =
        {
          {"verbose", no_argument,       &verbose_flag, 1},
          {"rho",     no_argument,       0, 'r'},
          {"p-1",  no_argument,       0, 'p'},
          {"division",  no_argument,       0, 'd'},
          {"help",  no_argument,       0, 'h'},
          {"verbose",  no_argument,       0, 'v'},
          {"brent",     no_argument,       0, 'b'},
          {"parameter", required_argument,       0, 'b'},
          {"file",    required_argument, 0, 'f'},
          {0, 0, 0, 0}
        };
      int option_index = 0;

      c = getopt_long (argc, argv, "rpdhvbi:f:",
                       long_options, &option_index);
      if (c == -1)
        break;

      switch (c)
        {
        case 0:
          /* If this option set a flag, do nothing else now. */
          if (long_options[option_index].flag != 0)
            break;
          printf ("option %s", long_options[option_index].name);
          if (optarg)
            printf (" with arg %s", optarg);
          printf ("\n");
          break;

        case 'r':
          choix = 0;
          break;

        case 'p':
          choix = 3;
          break;

        case 'd':
          choix = 1;
          break;
        case 'b':
          choix = 2;
          break;

        case 'h':
          hFlag = 1;
          break;

        case 'f':
          fFlag = 1;
          strcpy(filename,optarg);
          break;
        case 'i':
          iFlag = 1;
          parameter=atoi(optarg);
          break;
        case 'v':
          verbose_flag = 1;
          break;
        case '?':
          puts("Erreur");
          break;
        default:
          abort ();
        }
    }
  if (verbose_flag)
    puts ("verbose flag is set");
  if (optind < argc)
    {
      printf ("non-option ARGV-elements: ");
      while (optind < argc)
        printf ("%s ", argv[optind++]);
      putchar ('\n');
    }
  if ((parameter > 1000000 || parameter < 1) && iFlag==1) {
    printf("Erreur paramètre incorrect non compris entre 1 et 1000000\n");
    return 0;
  }
  if(hFlag || fFlag!=1){
    printf("Use: fact -f [FILE] [OPTION]\n\n");
    printf("  By default optimized factorization depending of size number\n");
    printf("  -d Factorization by successive divisions\n");
    printf("  -p Factorization by the method p-1 of Pollard \n");
    printf("  -r Factorization by the method rho of Pollard \n");
    printf("  -b Factorization by the method rho of Pollard Brent variant\n");
    printf("  -i Parameter that will set borne larger and number of iterations\n\n");
    printf("Example:  \n");
    printf("  ./fact -f 10chiffres.txt\n");
    printf("  ./fact -f 70chiffres.txt -p\n");
    return 0;
  }
  if(iFlag==0){
    parameter=1000000;
  }
  /*Initialization */
  mpz_t n,r;
	mpz_init(n);
  mpz_init(r);
  mpz_ui_pow_ui(r,10,SIZE_NUMBER);

  /*File opening*/
  FILE *fp;
  char *line = NULL;
  size_t len = 0;
  ssize_t read;
  /*fp = fopen(filename, "w");
  if (fp == NULL)
    exit(EXIT_FAILURE);

  for(i=0;i<TEST_MAX;i++){
    mpz_urandomm(random,rstate,r);
    mpz_out_str(fp,10,random);
    fprintf(fp, "\n");
  }
  fclose(fp);*/
  fp = fopen(filename, "r");
  if (fp == NULL)
    exit(EXIT_FAILURE);
  /*Timer Set*/
  clock_t begin = clock();
  int notPrime=0;
  int maybePrime=0;
  int factParfaites=0;
  int resultat;
  /*Read the file*/
  while ((read = getline(&line, &len, fp)) != -1){
    mpz_set_str(n, line,10);
    printf("N = %s",line);
    printf("Factorization p = 1x");
    /*Call of Factorization functions depending of the choice of the user*/
    if(choix!=-1){
      resultat = Factorisation(n,choix,parameter);
      switch (resultat) {
        case 100:
            notPrime++;
            printf("Warning!: Last factor is not Prime !!!\n");
            break;
        case 10:
            maybePrime++;
            printf("Warning!: The last factor was fund by default because the algorithm didn't find any other factors. it is not necessary prime! \n");
            break;
        default:
          factParfaites++;

      }
    }else{

      resultat = FactorisationOpti(n,parameter);

      switch (resultat) {
        case 100:
            notPrime++;
            printf("Warning!: Last factor is not Prime !!!\n");
            break;
        case 10:
            maybePrime++;
            printf("Warning!: The last factor was fund by default because the algorithm didn't find any other factors. it is not necessary prime! \n");
            break;
        default:
          factParfaites++;
      }
    }
    printf("\n\n");

    mpz_set_ui(r,0);
  }
  clock_t end = clock();
  /*Print of Execution Time*/
  double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
  printf("Execution Time %fs\n",time_spent);
  printf("Factorisations parfaites %d\n",factParfaites);
  printf("Dernier facteur pas forcement Premier %d\n",maybePrime);
  printf("Dernier facteur non premier %d\n",notPrime);
  fclose(fp);
  gmp_randclear(rstate);
  mpz_clear(n);
  mpz_clear(r);
  mpz_clear(random);

  return 0;
}
