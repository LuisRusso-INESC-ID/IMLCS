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

#include <stdlib.h>
#include <assert.h>
#include <stdio.h>

#include "pointHash.h"

static int primes[] = {
  3, 5, 7, 11, 17, 29, 53, 97, 193, 389, 769, 1543, 3079, 6151, 12289, 24593,
  49157, 98317, 196613, 393241, 786433, 1572869, 3145739, 6291469,
  12582917, 25165843, 50331653, 100663319, 201326611, 402653189, 805306457,
  1610612741
};

struct pointHash{
  point *A; /* Array for storing points */
  int M; /* Size of A */
  int n;  /* Number of points inside A */
  int dim; /* Point dimension */
};

#ifndef NDEBUG

static void
printTable(pointHash h
           )
{
  for(int i = 0; i < h->M ; i++){
    printf("[%d] = ", i);
    if(NULL != h->A[i]){
      printf("%d ", *h->A[i]);
    }
    printf("; ");
  }
  printf("\n");
}
#endif /* NDEBUG */

pointHash
allocPH(int dim
        )
{
  pointHash h = NULL;

  h = calloc(1, sizeof(struct pointHash));
  h->dim = dim;
  h->M = 3;
  h->A = calloc(h->M, sizeof(point));

  return h;
}

void
freePH(pointHash h
       )
{
  for(int i = 0; i < h->M; i++){
    if(NULL != h->A[i])
      free(h->A[i]);
  }
  free(h->A);
  free(h);
}

static unsigned int
hash(point p,
     int dim,
     int M
     )
{
  unsigned int r = 0;
  unsigned int a = 31415;
  const unsigned int b = 27183;

  unsigned char *S = (unsigned char *)p;

  for(int i = 0; i < dim*4 ; i++){
    r = (a*r + *S) % M;
    S++;
    a = (a*b) % (M-1);
  }


  return r;
}

static int
findPosition(pointHash h,
             point p
             )
{
  int i = hash(p, h->dim, h->M);

  while(NULL != h->A[i]
        && !pointEquals(h->A[i], p, h->dim)){
    i++;
    i %= h->M;
  }

  return i;
}

static void
expandPH(pointHash h,
         int M /* New size */
         )
{
  int oM = h->M; /* Old table size */
  point *A = h->A;

  h->M = M;
  h->A = calloc(h->M, sizeof(point));
  h->n = 0;

  for(int i = 0; i < oM; i++){
    if(NULL != A[i])
      insertPH(h, A[i]);
  }

  free(A);
}

void
insertPH(pointHash h,
         point p /* A pointer to the point */
         )
{
#ifndef NDEBUG
  int ni = h->n;
#endif /* NDEBUG */

  if(2*(h->n+1) > h->M){ /* Expand the table */
    int i;
    for(i = 0; primes[i] <= h->M; i++)
      ;
    expandPH(h, primes[i]);
  }

  h->n++;
  assert(1 + ni == h->n && "Lost point on insert.");

  h->A[findPosition(h, p)] = p;
}

int
containsPH(pointHash h,
           point p
           )
{
  /* printf("Searching %d\n", *p); */
  /* printTable(h); */

  return NULL != h->A[findPosition(h, p)];
}

void
deletePH(pointHash h,
         point p /* A pointer to the point */
         )
{
  /* printf("Deleting %d\n", *p); */
  /* printTable(h); */

  int i;
#ifndef NDEBUG
  int ni = h->n;
#endif /* NDEBUG */

  if(3-1 < h->n && 8*(h->n-1) < h->M){ /* Reduce table */
    for(i = 0; primes[i] < h->M; i++)
      ;
    expandPH(h, primes[i-1]);
  }

  i = findPosition(h, p);
  h->A[i] = NULL; /* Emptied position */
  h->n--;
  i++;
  i %= h->M;
  while(NULL != h->A[i]){
    point t = h->A[i];
    h->A[i] = NULL; /* Emptied position */
    h->n--;
    insertPH(h, t);
    i++;
    i %= h->M;
  }

  /* printTable(h); */
  assert(!containsPH(h, p) && "Hash delete failed.");
  assert(ni - 1 == h->n && "Failed Hash delete.");
}

#if 0 /* Comment */

int
main(__attribute__((unused)) int argc,
     __attribute__((unused)) char** argv)
{
  pointHash h;
  point p;
  const int limit = 1000;


  h = allocPH(1);

  for(int i = 0; i < limit; i++){
    p = calloc(1, sizeof(int));
    *p = i;
    insertPH(h, p);
  }

  p = calloc(1, sizeof(int));
  for(int i = 0; i < limit; i+=2){
    *p = i;
    deletePH(h, p);
  }

  for(int i = 1; i < limit; i+=2){
    *p = i;
    assert(containsPH(h, p) && "Ups lost a point.");
    *p = i-1;
    assert(!containsPH(h, p) && "Ups ghost point.");
  }

  freePH(h);

  return 0;
}

#endif
