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

#include "ort.h"
#include "string.h"
#include "mlcs.h"
#include "naivemlcs.h"

/* Test code by comparing with random strings */
void
unitTest(int sigma, /* Alphabet size */
	 int n,     /* Average string length */
	 int dim,   /* Number of strings */
	 int ops    /* Number of operations */
	 )
{
#ifndef NDEBUG
  char *S[dim]; /* Array with the strings */
#endif /* NDEBUG */
  mlcs m = allocMLCS(dim, sigma);

  for(; 0 < ops; ops--){
    int t;
    t = arc4random_uniform(dim);  /* Choose a string */
    int option = arc4random_uniform(2); /* choose option */
    if(stringSize(pullString(m,t)) < n)
      option = 0; /* Insert */
    if(stringSize(pullString(m,t)) > 2*n)
      option = 1; /* Delete */

    switch(option){
    default:
    case 0: /* Insert */
      append(m, t, 'A' + arc4random_uniform(sigma));
      break;
    case 1: /* Delete */
      pop(m, t);
      break;
    }

#ifndef NDEBUG
    for(int j = 0; j < dim; j++)
      S[j] = stringGet(pullString(m,j));

    printf("[CHECK] %d = %d\n",
	   mlcsSize(m),
	   naiveMLCS(dim, S, sigma)
	   );

    for(int j = 0; j < dim; j++)
      printf("[CHECK] %d : %s\n", j, S[j]);

    assert(mlcsSize(m) == naiveMLCS(dim, S, sigma));

    for(int j = 0; j < dim; j++)
      free(S[j]);
#endif /* NDEBUG */
  }

  freeMLCS(m);
}

int
commandShell(void)
{
  int resets = 0;
  int count = 0;
  char C = 'A';
  mlcs m = NULL;
  int t;
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
      if(NULL != m)
	freeMLCS(m);
      scanf("%d", &dim);
      scanf("%d", &sigma);
      m = allocMLCS(dim, sigma);
      break;
    case 'I': /* Insert letter */
      count++;
      scanf("%d", &t);
      do{
	C = getchar();
      } while(' ' == C);
      append(m, t, C);
      break;
    case 'D': /* Delete letter */
      count++;
      scanf("%d", &t);
      pop(m, t);
      break;
    }
#ifndef NDEBUG
    if('I' == C || 'D' == C){
      char *S[dim]; /* Array with the strings */

      for(int j = 0; j < dim; j++)
        S[j] = stringGet(m->S[j]);

      assert(mlcsSize(m) == naiveMLCS(dim, S, sigma));

      printf("[CHECK] %d = %d\n",
             mlcsSize(m),
             naiveMLCS(dim, S, sigma)
             );

      for(int j = 0; j < dim; j++)
        printf("[CHECK] %d : %s\n", j, S[j]);

      for(int j = 0; j < dim; j++)
        free(S[j]);
    }
#endif /* NDEBUG */
    clock_gettime(CLOCK_MONOTONIC, &stop);
    cpu_time_used = (stop.tv_sec - start.tv_sec);
    /* printf("diff %d\n", cpu_time_used); */
  }

  if(NULL != m)
    freeMLCS(m);

  return count;
}

#if 1 /* Comment */

int
main(__attribute__((unused)) int argc,
     __attribute__((unused)) char** argv)
{
  int count = 0;
  /* int minCount = 10; */
  int minCount = 1000000;

  if(1 < argc)
    sscanf(argv[1], "%d", &minCount);

  adjustCut(-CUTOFF);
  /* configure(1<<30, 60, 4.0/5.0, CUTOFF); */

  count = commandShell();
  printf("%d ", count); /* repetitions */

  adjustCut(0); /* Release internal array */
  /* deconfigure(); */

  return count < minCount ? 1 : 0;
}

#endif

#if 0 /* Comment */

int
main(__attribute__((unused)) int argc,
     __attribute__((unused)) char** argv)
{
  FILE *f = stdout;

  stdout = fopen("outf", "w");
  setbuf(stdout, NULL);

  adjustCut(-CUTOFF);
  /* configure(1<<30, 60, 4.0/5.0, CUTOFF); */
  /* configure(1<<30, 60, 4.0/5.0, 0); */

  /* char **S = malloc(2*sizeof(char*)); */
  /* S[0] = strdup("ABBBA"); */
  /* S[1] = strdup("BAAB"); */
  /* naiveMLCS(2, S, 2); */

  /* printf("First test\n"); */
  /* unitTest(1, 10, 2, 10); */
  /* unitTest(1, 10, 2, 100); */
  /* unitTest(1, 10, 4, 100); */
  /* unitTest(1, 20, 4, 400); */

  /* printf("Second test\n"); */
  /* unitTest(2, 3, 2, 100); */
  /* unitTest(2, 3, 3, 100); */
  /* unitTest(2, 3, 7, 100); */
  /* unitTest(3, 10, 4, 100); */
  /* unitTest(3, 20, 4, 400); */
  /* unitTest(4, 20, 4, 400); */
  /* unitTest(4, 3, 7, 100); */
  /* unitTest(10, 40, 4, 5); */
  /* unitTest(10, 20, 4, 400); */

  /* printf("Other tests\n"); */
  /* unitTest(4, 20, 4, 80); */
  /* unitTest(20, 20, 4, 100); */

#define K 4

  char T[K][100] = {
		    "BBBABAAAAABBBACAABCBB",
		    "CAACACACBABBACBCAC",
		    "ACCBACABBACCCBABACCA",
		    "ACAAAACBBACAABCCCCCB"
  };

  mlcs m = allocMLCS(K, 3);

  for(int i = 0; i < K; i++){
    int l = strlen(T[i]);
    append(m, i, 'A');
    pop(m,i);
    for(int j = 0; j < l; j++){
      append(m, i, T[i][j]);
    }
  }

  printMLCS(m);
  pop(m, 3);
  append(m, 2, 'C');
  printf("After\n");
  printMLCS(m);

  /* append(m, 0, 'A'); */
  /* append(m, 0, 'A'); */
  /* append(m, 0, 'A'); */
  /* append(m, 0, 'A'); */
  /* append(m, 0, 'A'); */
  /* append(m, 0, 'A'); */

  /* append(m, 1, 'A'); */
  /* append(m, 1, 'A'); */
  /* append(m, 1, 'A'); */
  /* append(m, 1, 'A'); */
  /* append(m, 1, 'A'); */
  /* append(m, 1, 'A'); */
  /* append(m, 1, 'A'); */
  /* append(m, 1, 'A'); */
  /* append(m, 1, 'A'); */

  /* append(m, 2, 'A'); */
  /* append(m, 2, 'A'); */
  /* append(m, 2, 'A'); */
  /* append(m, 2, 'A'); */
  /* append(m, 2, 'A'); */
  /* append(m, 2, 'A'); */

  /* append(m, 3, 'A'); */
  /* append(m, 3, 'A'); */
  /* append(m, 3, 'A'); */
  /* append(m, 3, 'A'); */
  /* append(m, 3, 'A'); */
  /* append(m, 3, 'A'); */
  /* append(m, 3, 'A'); */
  /* append(m, 3, 'A'); */
  /* append(m, 3, 'A'); */
  /* append(m, 3, 'A'); */

  /* pop(m, 3); */

  /* mlcs m = allocMLCS(2, 2); */

  /* append(m, 0, 'A'); */
  /* append(m, 0, 'A'); */
  /* append(m, 0, 'B'); */
  /* append(m, 0, 'A'); */
  /* append(m, 0, 'A'); */

  /* append(m, 1, 'A'); */
  /* append(m, 1, 'A'); */
  /* append(m, 1, 'A'); */
  /* append(m, 1, 'B'); */
  /* append(m, 1, 'B'); */

  /* pop(m, 0); */

  /* printf("MLCS has size %d\n", m->lambda); */

  /* for(int i = 0; i < nn; i++) */
  /*   pop(m, 0); */

  /* printf("MLCS has size %d\n", m->lambda); */

  /* for(int i = 0; i < n; i++) */
  /*   pop(m, 1); */

  /* printf("MLCS has size %d\n", m->lambda); */

  /* return 0; */

 /* mlcs m = allocMLCS(3); */

 /*  n = nn = 20; */

 /*  for(int i = 0; i < n; i++) */
 /*    append(m, 1, 'a'); */

 /*  for(int i = 0; i < nn; i++){ */
 /*    append(m, 0, 'a'); */
 /*    append(m, 0, 'b'); */
 /*  } */

 /*  for(int i = 0; i < nn/2; i++){ */
 /*    append(m, 2, 'd'); */
 /*    append(m, 2, 'a'); */
 /*    append(m, 2, 'c'); */
 /*  } */

 /*  printf("MLCS has size %d\n", mlcsSize(m)); */

 /*  for(int i = 0; i < nn; i++) */
 /*    pop(m, 0); */

 /*  printf("MLCS has size %d\n", mlcsSize(m)); */

 /*  for(int i = 0; i < n; i++) */
 /*    pop(m, 1); */

 /*  printf("MLCS has size %d\n", mlcsSize(m)); */

  freeMLCS(m);
  fclose(stdout);
  stdout = f;

  /* deconfigure(); */
  adjustCut(0);

  return 0;
}

#endif	/* Comment */
