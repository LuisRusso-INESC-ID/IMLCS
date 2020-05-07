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
 *  \brief     Dynamic Multiple Longest Common Sub-sequence
 *  \details   This file implements the dynamic MLCS algorithm
 *  \author    Luís M. S. Russo
 *  \version   0.1.0-alpha
 *  \date      04-05-2020
 *  \copyright BSD 2-Clause License
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>

#include "ort.h"
#include "string.h"
#include "mlcs.h"
#include "naivemlcs.h"
#include "pointQueue.h"
#include "pointHash.h"

/* A structure for storing a multiple longest common sub-string. */
#ifdef NDEBUG
struct mlcs{
  int dim; /* The number of strings considered. */
  int zeros;  /* Number of strings of size 0 */
  string *S; /* Array that stores the current strings.  */
  int lambda; /* The size of the mlcs */
  ort *PF;    /* Array with pareto fronts. */
  int pfA;    /* Size of the PF array */
};
#endif /* NDEBUG */

#ifndef NDEBUG
static void
checkMLCS(mlcs m)
{
  point *T;
  int n;
  int A[m->dim];
  point p = A;
  int B[m->dim];
  point z = B;

  for(int k = 0; k < m->dim; k++)
    z[k] = -1;

  for(int i = 1; i <= m->lambda; i++){
    T = collect(m->PF[i], z, &n);
    assert(0 < n && "Empty layer");
    for(int j = 0; j < n; j++){
      assert(0 < countQ(m->PF[i-1], T[j]) && "Unjustified point in layer");

      for(int k = 0; k < m->dim; k++)
	p[k] = T[j][k]+1;

      assert(1 == countQ(m->PF[i], p) && "Non minimal point in layer");
      free(T[j]);
    }
    free(T);
  }
}
#endif /* NDEBUG */


void
printMLCS(mlcs m)
{ /* Print the data structure for debuging. */
  point *T;
  int n;
  int A[m->dim];
  point p = A;

  for(int k = 0; k < m->dim; k++)
    p[k] = -1;

  for(int i = 0; i <= m->lambda; i++){
    printf("\n");
    printf("@ %d : ", i);
    n = 0;
    T = collect(m->PF[i], p, &n);
    for(int j = 0; j < n; j++){
      printf("(");
      for(int k = 0; k < m->dim; k++)
	printf("%d, ", T[j][k]);
      printf(") ");
      free(T[j]);
    }
    printf("\n");
    free(T);
  }
}

/* Allocs an MLCS of dimension dim. The first string is set to S.
  All other strings are kept empty. */
mlcs
allocMLCS(int dim, /* The number of strings. */
	  int sigma /* Alphabet size */
	  )
{
  assert(1 < dim && "Error: MLCS should contain at least 2 strings.");
  mlcs r = malloc(sizeof(struct mlcs));

  r->dim = dim;
  r->zeros = dim;
  r->lambda = 0;
  r->pfA = 3;
  r->PF = calloc(r->pfA, sizeof(ort));

  r->S = calloc(dim, sizeof(string));
  for(int i = 0; i < dim; i++)
    r->S[i] = stringAlloc(sigma);

  point p = malloc(dim*sizeof(int));
  for(int i = 0; i < dim; i++)
    p[i] = -1;

  r->PF[0] = allocORT(dim);
  insert(r->PF[0], p);
  free(p); /* Can free it after insert */

  return r;
}

void
freeMLCS(mlcs m /* The structure */
	 )
{
  assert(NULL != m && "Error: freeing NULL MLCS.");

  for(int i = 0; i < m->dim; i++)
    stringFree(m->S[i]);

  free(m->S);
  m->S = NULL;

  for(int j = 0; j < m->pfA; j++){
    if(NULL != m->PF[j])
      freeORT(m->PF[j]);
    m->PF[j] = NULL;
  }
  free(m->PF);

  m->dim = 0;
  m->lambda = 0;
  m->pfA = 0;

  free(m);
}

string
pullString(mlcs m,
	   int t /* Which string do you want */
          )
{
  return m->S[t];
}

void
append(mlcs m, /* The MLCS data struct */
       int j,  /* Which string */
       char c  /* Which letter */
       )
{
  string S = m->S[j];
  if(0 == stringSize(S))
    m->zeros--; /* Another string gets a size */

  if(0 == m->zeros){
    if(1+m->lambda == m->pfA){
      m->pfA *= 2;
      m->PF = realloc(m->PF, m->pfA*sizeof(ort));
      bzero(&(m->PF[m->lambda+1]),
	    (m->pfA/2)*sizeof(ort)
	    );
    }

    if(NULL == m->PF[m->lambda+1])
      m->PF[m->lambda+1] = allocORT(m->dim);

    int dim = m->dim;
    /* candidate point */
    point p = malloc(dim*sizeof(int));
    /* Lower point */
    point z = malloc(dim*sizeof(int));
    /* further right point */
    point f = malloc(dim*sizeof(int));
    for(int i = 0; i < dim; i++){
      z[i] = -1;
      f[i] = stringLast(m->S[i], c);
    }
    z[j] = stringLast(m->S[j], c);
    f[j] = stringEnd(m->S[j]);

    for(int r = 0; r <= m->lambda; r++){
      int n;
      /* point *T = collect(m->PF[r], z, &n); */
      point *T = rangeCollect(m->PF[r], z, f, &n);
      int excludeP;

      if(0 < n){
	ort tempP = allocORT(dim); /* Temporary points */
	for(int k = 0; k < n; k++){ /* Process points */
	  excludeP = 0;
	  for(int i = 0; !excludeP && i < dim; i++){
	    p[i] = stringCeil(m->S[i], c, T[k][i]+1);
	    if(i == j){
	      p[i] = stringEnd(m->S[i]);
	      excludeP = T[k][i] == p[i];
	    } else
	      excludeP = -2 == p[i];
	  }

	  /* if(!excludeP){ */
	  /*   printf("\n [CHECK] j: %d >> ( ", j); */
	  /*   for(int i = 0; i < dim; i++) */
	  /*     printf("%d < %d, ", T[k][i], p[i]); */
	  /*   printf(")\n"); */
	  /* } */

	  free(T[k]);
	  if(!excludeP){
	    /* Check to see if it is dominated on the PF */
            plusPlus(p, dim);
	    excludeP = 0 < countQ(m->PF[r+1], p);
	    if(!excludeP) /* Avoid duplicate inserts */
	      excludeP = 0 < countQ(tempP, p);
            minusMinus(p, dim);
	    if(!excludeP)
	      insert(tempP, p);
	  }
	}

	free(T);
	T = collect(tempP, z, &n);
	for(int k = 0; k < n; k++){ /* Process points */
	  /* Check to see if it is a minima */
	  plusPlus(T[k], dim);
	  excludeP = 1 < countQ(tempP, T[k]);
	  minusMinus(T[k], dim);
	  if(!excludeP)
	    insert(m->PF[r+1], T[k]);
	  free(T[k]);
	}
	freeORT(tempP);
      }
      free(T);
    }

    if(0 < weightORT(m->PF[m->lambda+1]))
       m->lambda++;

    free(f);
    free(z);
    free(p);
  }

  stringAppend(S, c);

#ifndef NDEBUG
  printMLCS(m);
  checkMLCS(m);
#endif /* NDEBUG */
}

#ifndef NDEBUG
static void
gdbBreak(void) {}
#endif /* NDEBUG */

static void
uncover(pointQueue M, /* For removing non-minima. */
	pointHash CleanM,
	mlcs m,
        char c,
        int r,
        int j, /* Important only when 0 == r */
        int *pprev, /* The point of previous pos */
        int *baseTop, /* The point that is moving */
	int *baseBot /* The far away point */
        )
{
  int dim = m->dim;

  if(1 == r){
    if(baseBot[j] > baseTop[j]){
      point p = malloc(dim*sizeof(int));
      memcpy(p, baseTop, dim*sizeof(int));

      p[j] = stringCeil(m->S[j], c, p[j]+1);

      plusPlus(p, dim);
      int insertQ = 1 == countQ(m->PF[r], p);
      minusMinus(p, dim);

      if(insertQ){
	insert(m->PF[r], p);
	pushPQ(M, p);
	insertPH(CleanM, p);
	/* printf("MIN REMOVAL >>>>>> Pushed %d %d\n", p[0], p[1]); */
      } else
	free(p);
    }
  } else {
    pointHash H = allocPH(dim); /* Register of considered points */

    for(int i = 0; i < dim; i++){
      int h = baseBot[i];
      baseBot[i] = baseTop[i]; /* Set one coordinate to Top. */
      /* slice : Part of the previous frontier */
      int n; /* slice size */
      /* point *slice = dominatedCollect(m->PF[r-1], baseBot, &n); */
      point *slice = rangeCollect(m->PF[r-1], pprev, baseBot, &n);

      for(int k = 0; k < n; k++){
	if(!containsPH(H, slice[k])){
	  insertPH(H, slice[k]);

	  point p = malloc(dim*sizeof(int));
	  for(int l = 0; l < dim; l++)
	    p[l] = stringCeil(m->S[l], c, slice[k][l]+1);

	  /* First check All equal */
	  int insertQ = 1;
	  for(int l = 0; insertQ && l < dim; l++)
	    insertQ = baseTop[l] == p[l];

	  /* Then check forward move */
	  insertQ = !insertQ;
	  /* for(int l = 0; insertQ && l  < dim; l++) */
	  /*   insertQ = baseTop[l] <= p[l]; */

	  /* Now check Non-dominated */
	  plusPlus(p, dim);
	  insertQ = insertQ && (1 == countQ(m->PF[r], p));
	  minusMinus(p, dim);
	  /* insertQ = insertQ && !containsQ(m->PF[r], p); */

	  if(insertQ){

	  /*   if(3 == r && */
	  /*      10 == p[0] && */
	  /*      8 == p[1] && */
	  /*      7 == p[2] && */
	  /*      7 == p[3] */
	  /*      ) gdbBreak(); */

	    insert(m->PF[r], p);
	    pushPQ(M, p);
	    insertPH(CleanM, p);
	    /* printf("MIN REMOVAL >>>>>> Pushed %d %d\n", p[0], p[1]); */
	  } else
	    free(p);
	} else
	  free(slice[k]);
      }
      if(NULL != slice)
	free(slice);
      baseBot[i] = h;
    }

    freePH(H); /* Frees all the points involved */
  }
}

static void
cleanCovered(pointQueue M,
	     int dim,
	     ort t,
	     point p
	     )
{
  while(!isEmptyPQ(M)){
    memcpy(p, topPQ(M), dim*sizeof(int));
    /* printf("CLEAN REMOVAL >>>>>> Pushed %d %d\n", p[0], p[1]); */
    plusPlus(p, dim);
    if(1 < countQ(t, p))
      delete(t, topPQ(M));
    popPQ(M);
  }
}

/* Remove the first letter of a given string. */
void
pop(mlcs m, /* The MLCS data struct */
    int j  /* Which string */
    )
{
  string S = m->S[j];
  if(1 == stringSize(S))
    m->zeros++;

  if(0 < m->lambda){
    int dim = m->dim;
    int n ;
    point *T = NULL;

    /* Secondary point queue for cleanning non minima */
    pointQueue M = allocPQ();
    expandPQ(M, stringSigma(S));
    /* Hash to store points for cleaning */
    pointHash CleanM = allocPH(dim);

    /* The main Queue of the algorithm */
    pointQueue Q = allocPQ();

    point p = malloc(dim*sizeof(int));
    /* Which letter are you removing ? */
    char c = stringFstLetter(S);
    int valid = 1;
    for(int i = 0; valid && i < dim; i++){
      p[i] = stringCeil(m->S[i], c, 0);
      valid = 0 <= p[i];
    }
    if(valid){
      pushPQ(Q, p);
      markPQ(Q);
      p = malloc(dim*sizeof(int));
    }
    int r = 1;
    /* Queue Load Complete */

    while(!isEmptyPQ(Q)){ /* There is stuff in the queue. */
      if(markingTruePQ(Q)){ /* Frontier transition. */
	cleanCovered(M, dim, m->PF[r], p);
        r++; /* Update the pareto index */
        markPQ(Q);
      }

      /* 0. Uncover hidden points. */
      memcpy(p, topPQ(Q), dim*sizeof(int));
      char c = stringLetter(S, p[j]);

      /* Point with large indexes */
      int *plI = malloc(dim*sizeof(int));
      /* previous positions of p */
      int *pprev = malloc(dim*sizeof(int));
      for(int i = 0; i < dim; i++){
        plI[i] = stringLast(m->S[i], c);
	pprev[i] = stringFloor(m->S[i], c, p[i]-1);
      }
      uncover(M, CleanM, m, c, r, j, pprev, p, plI);
      free(pprev);
      free(plI);

      /* 1. Collect points that dominate current. */
      n = 0;
      if(r+1 < m->pfA && NULL != m->PF[r+1])
        T = collect(m->PF[r+1], topPQ(Q), &n);

      /* 2. remove the point. */
      /* printMLCS(m); */
      delete(m->PF[r], topPQ(Q));
      /* printMLCS(m); */

      /* 3. filter candidates */
      for(int i = 0; i < n; i++){
        if(0 == countQ(m->PF[r], T[i])){
          pushPQ(Q, T[i]); /* Mark for deletion */
          pushPQ(M, T[i]); /* Verify non-minima later */

	  /* if(3 == r && */
	  /*    10 == T[i][0] && */
	  /*    8 == T[i][1] && */
	  /*    7 == T[i][2] && */
	  /*    7 == T[i][3] */
	  /*    ) gdbBreak(); */

          insert(m->PF[r], T[i]);
          /* printMLCS(m); */
        } else /* Free elements that do not go into Q */
          free(T[i]);
      }
      if(NULL != T){
        free(T);
        T = NULL;
      }

      /* 4. process Q */
      free(topPQ(Q));
      popPQ(Q);
    }
    cleanCovered(M, dim, m->PF[r], p);
    free(p);

    freePQ(Q);

    freePH(CleanM); /* Frees all the points involved */
    freePQ(M);
    if(NULL != T){
      free(T);
      T = NULL;
    }

    if(0 == weightORT(m->PF[m->lambda]))
      m->lambda--;
  }

  stringPop(S); /* Remove the letter */

#ifndef NDEBUG
  printMLCS(m);
  checkMLCS(m);
#endif /* NDEBUG */
}

int
mlcsSize(mlcs m
	 )
{
  return m->lambda;

}
