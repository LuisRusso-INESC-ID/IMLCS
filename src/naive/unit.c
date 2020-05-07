/* BSD 2-Clause License */

/* Copyright (c) 2020, Luís M. S. Russo */
/* All rights reserved. */

/* Redistribution and use in source and binary forms, with or without */
/* modification, are permitted provided that the following conditions are met: */

/* 1. Redistributions of source code must retain the above copyright notice, this */
/*    list of conditions and the following disclaimer. */

/* 2. Redistributions in binary form must reproduce the above copyright notice, */
/*    this list of conditions and the following disclaimer in the documentation */
/*    and/or other materials provided with the distribution. */

/* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" */
/* AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE */
/* IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE */
/* DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE */
/* FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL */
/* DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR */
/* SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER */
/* CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, */
/* OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE */
/* OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. */

/**
 *  \brief     Unit tests
 *  \details   Testing the dynamic MLCS algorithm
 *  \author    Luís M. S. Russo
 *  \version   0.1.0-alpha
 *  \date      04-05-2020
 *  \copyright BSD 2-Clause License
 */

#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <bsd/stdlib.h>

#include "string.h"
#include "naivemlcs.h"

int
commandShell(void)
{
  int resets = 0;
  int count = 0;
  char LC = 'I'; /* Last change command */
  char C = 'A';
  int t;
  string *A = NULL;
  int dim;
  int sigma; /* Alphabet size */
  struct timespec start, stop;
  int cpu_time_used = 0;

  clock_gettime(CLOCK_MONOTONIC, &start);
  while('X' != C && cpu_time_used < TIME_LIMIT){
    C = getchar();
    switch(C){
    case 'K': /* Define number of strings */
      resets++;
      if(NULL != A){
        char *S[dim]; /* Array with the strings */

        for(int j = 0; j < dim; j++)
          S[j] = stringGet(A[j]);
	/* printf("%d\n", */
	naiveMLCS(dim, S, sigma)
		 /* ) */
	;
        for(int j = 0; j < dim; j++)
          free(S[j]);

        for(int j = 0; j < dim; j++)
          stringFree(A[j]);
        free(A);
      }
      scanf("%d", &dim);
      scanf("%d", &sigma);
      A = calloc(dim, sizeof(string));
      for(int j = 0; j < dim; j++)
	A[j] = stringAlloc(sigma); /* Alphabet is not relevant */

      break;
    case 'I': /* Insert letter */
      count++;
      scanf("%d", &t);
      do{
	C = getchar();
      } while(' ' == C);
      stringAppend(A[t], C);

      if('D' == LC){
        LC = 'I';
        char *S[dim]; /* Array with the strings */

        for(int j = 0; j < dim; j++)
          S[j] = stringGet(A[j]);
	/* printf("%d\n", */
	naiveMLCS(dim, S, sigma)
		 /* ) */
	;
        for(int j = 0; j < dim; j++)
          free(S[j]);
      }
      break;
    case 'D': /* Insert letter */
      count++;
      scanf("%d", &t);
      stringPop(A[t]);

      if('I' == LC){
        LC = 'D';
        char *S[dim]; /* Array with the strings */

        for(int j = 0; j < dim; j++)
          S[j] = stringGet(A[j]);
	/* printf("%d\n", */
	naiveMLCS(dim, S, sigma)
		 /* ) */
	;
        for(int j = 0; j < dim; j++)
          free(S[j]);
      }
      break;
    }
    clock_gettime(CLOCK_MONOTONIC, &stop);
    cpu_time_used = (stop.tv_sec - start.tv_sec);
  }

  if(NULL != A){
    char *S[dim]; /* Array with the strings */

    for(int j = 0; j < dim; j++)
      S[j] = stringGet(A[j]);
	/* printf("%d\n", */
    naiveMLCS(dim, S, sigma)
		 /* ) */
	;
    for(int j = 0; j < dim; j++)
      free(S[j]);

    for(int j = 0; j < dim; j++)
      stringFree(A[j]);
    free(A);
  }

  return count;
}

int
main(__attribute__((unused)) int argc,
     __attribute__((unused)) char** argv)
{
  int count = 0;
  /* int minCount = 10; */
  int minCount = 1000000;

  if(1 < argc)
    sscanf(argv[1], "%d", &minCount);

  count = commandShell();
  printf("%d ", count); /* repetitions */

  return count < minCount ? 1 : 0;
}
