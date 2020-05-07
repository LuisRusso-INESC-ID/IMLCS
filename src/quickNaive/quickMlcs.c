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
 *  \brief     Quick Multiple Longest Common Sub-sequence
 *  \details   This file implements the dynamic MLCS algorithm
 *  \author    Luis Russo
 *  \version   0.1
 *  \date      19-12-2019
 *  \copyright Simplified BSD License
 */

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <bsd/stdlib.h>
#include <limits.h>

#include "point.h"

#ifndef NDEBUG
void
assertOrder(int m, int *A)
{
  for(int i = 0; i+1 < m; i++)
    assert(A[i] < A[i+1] && "Failed order call.");
}
#endif /* NDEBUG */

static void
swap(int *a, int *b)
{
  int t = *a;
  *a = *b;
  *b = t;
}

static void
partition(int Gdim, /* Global problem dimension */
	  point G, /* The list of points to consider.
		      Unsorted. */
	  int Ldim, /* Sub-problem dimension */
	  int *l, /* The left index, inclusive */
	  int *r, /* The right index, exclusive */
	  int p, /* The dividing pivot */
	  int *R, /* List of points */
	  int *iR /* indexes of R */
	  )
{
  int pl = *l;
  int pr = *r;

#ifndef NDEBUG
  assertOrder(*r-*l, R);
#endif /* NDEBUG */

  int i = *l;
  while(i < *r){
    if(G[R[iR[i]]*Gdim+Ldim-1] > p){
      (*r)--;
      swap(&iR[*r], &iR[i]);
    } else {
      if(G[R[iR[i]]*Gdim+Ldim-1] < p){
	swap(&iR[*l], &iR[i]);
	(*l)++;
      }
      i++;
    }
  }
#ifndef NDEBUG
  assertOrder(*r-*l, R);
#endif /* NDEBUG */

  for(int i = pl; i < *l; i++)
    assert(G[R[iR[i]]*Gdim+Ldim-1] <
	   G[R[iR[*l]]*Gdim+Ldim-1]
	   && "Failed order");

  for(int i = *l; i < *r; i++)
    assert(G[R[iR[i]]*Gdim+Ldim-1] ==
	   G[R[iR[*l]]*Gdim+Ldim-1]
	   && "Failed order");

  for(int i = *r; i < pr; i++)
    assert(G[R[iR[i]]*Gdim+Ldim-1] >
	   G[R[iR[*l]]*Gdim+Ldim-1]
	   && "Failed order");
}

/* Return a list of booleans that are 1
   iff the point R[i] is bellow position t */
char *
orderSelect(int Gdim, /* Global problem dimension */
	    point G, /* The list of points to consider.
			Unsorted. */
	    int Ldim, /* Sub-problem dimension */
	    int m, /* Number of points in sub-problem */
	    int t, /* Desired position */
	    int *R /* List of points in sub-problem,
		       ordered on coordinate 0. */
	    )
{
#ifndef NDEBUG
  assertOrder(m, R);
#endif /* NDEBUG */


  int l = 0; /* left index, inclusive */
  int r = m; /* Right index, exclusive */

  int *iR = malloc(m*sizeof(int));  /* Indexes to array R */
  for(int i = 0; i < m; i++)
    iR[i] = i;

  while(2 < r - l){
    /* Select median of 3 */
    swap(&iR[r-1], &iR[l+arc4random_uniform(r-0-l)]);
    swap(&iR[r-2], &iR[l+arc4random_uniform(r-1-l)]);
    swap(&iR[r-3], &iR[l+arc4random_uniform(r-2-l)]);

    if(G[R[iR[r-1]]*Gdim+Ldim-1] < G[R[iR[r-2]]*Gdim+Ldim-1])
      swap(&iR[r-1], &iR[r-2]);
    if(G[R[iR[r-2]]*Gdim+Ldim-1] < G[R[iR[r-3]]*Gdim+Ldim-1])
      swap(&iR[r-2], &iR[r-3]);
    if(G[R[iR[r-1]]*Gdim+Ldim-1] < G[R[iR[r-2]]*Gdim+Ldim-1])
      swap(&iR[r-1], &iR[r-2]);
    swap(&iR[r-1], &iR[r-2]); /* Put pivot on r-1 */

    /* Partition */
    int p = G[R[iR[r-1]]*Gdim+Ldim-1]; /* Pivot */
    int pl = l;
    int pr = r;
    partition(Gdim, G, Ldim, &pl, &pr, p, R, iR);

    /* Divide */
    if(pl <= t && t < pr){  /* Target hit */
      l = pl;
      r = pr;
      break;
    }
    else {
      if(pr <= t)
	l = pr;
      else if(pl > t)
	r = pl;
    }
  }

  /* You Now have a list with at most 2 elements */
  if(G[R[iR[r-1]]*Gdim+Ldim-1] < G[R[iR[l]]*Gdim+Ldim-1])
    swap(&iR[r-1], &iR[l]);
  if(G[R[iR[r-1]]*Gdim+Ldim-1] != G[R[iR[l]]*Gdim+Ldim-1])
    r = r-1;

  for(int i = 0; i < l; i++)
    assert(G[R[iR[i]]*Gdim+Ldim-1] <
	   G[R[iR[l]]*Gdim+Ldim-1]
	   && "Failed order");

  for(int i = l; i < r; i++)
    assert(G[R[iR[i]]*Gdim+Ldim-1] ==
	   G[R[iR[l]]*Gdim+Ldim-1]
	   && "Failed order");

  for(int i = r; i < m; i++)
    assert(G[R[iR[i]]*Gdim+Ldim-1] >
	   G[R[iR[l]]*Gdim+Ldim-1]
	   && "Failed order");

  int d = l; /* The division point */
  if(d < m-r)
    d = r;

  char *tB = NULL;

  if(0 < d && d < m){
    tB = malloc(m*sizeof(char)); /* store result */

    for(int i = 0; i < d; i++)
      tB[iR[i]] = 1;
    for(int i = d; i < m; i++)
      tB[iR[i]] = 0;
  }
  free(iR);

#ifndef NDEBUG
  assertOrder(m, R);
#endif /* NDEBUG */

  return tB;
}

static void
printCall(int Gdim,
	  point G,
	  char *B,
	  int m,
	  int *R,
	  char *killerR,
	  char *victimR
	  )
{
  for(int i = 0; i < m; i++){
    printf("i : %d \t G : ", i);
    for(int j = 0; j < Gdim; j++)
      printf("%6d,", G[R[i]*Gdim+j]);
    printf("\tB : %d", B[R[i]]);
    printf("\tK : %d", killerR[i]);
    printf("\tV : %d\n", victimR[i]);
  }
}

#ifndef NDEBUG
static void
assertMinimaR(int Gdim, /* Global problem dimension */
	      point G, /* The list of points to consider.
			  Sorted on coordinate 0. On ties sorted
			  on coordinate 1. On ties no further criteria. */
	      char *B, /* List for storing the result, it should indicate, for
			  the points in the sub-problem if they are minima or
			  not. This array is indexed as G, meaning that B[i] = 1
			  implies that the i-th point in G is a minima.  As
			  mentioned only if i = R[j] for some j is the result
			  meaningfull. */
	      /* Arguments related to sub-problem */
	      int Ldim, /* Sub-problem dimension */
	      int m, /* Number of points in sub-problem */
	      int *R, /* List of points in sub-problem,
			 ordered as G. */
	      char *killerR, /* 1 if R[i] point can kill others. */
	      char *victimR /* 1 if R[i] point can be kille. */
	      )
{
  for(int i = 0; i < m; i++){
    if(!victimR[i])
      assert(B[R[i]] && "Consistent result.");
    else {
      int count = 0;
      if(killerR[i])
	count--;

      for(int j = 0; 1 > count && j < m; j++){
	if(1 == killerR[j]){
	  int dominates = 1;
	  for(int k = 0; dominates && k < Ldim; k++){
	    dominates = G[R[i]*Gdim+k] >= G[R[j]*Gdim+k];
	  }
	  /* if(dominates) */
	  /*   printf("(%d -> %d)\n",i, j); */
	  count += dominates;
	}
	if(B[R[i]])
	  assert(count < 1 && "Consistent result.");
      }
    }
  }
}
#endif /* NDEBUG */

static void
minimaR( /* Arguments related to global problem */
	int Gdim, /* Global problem dimension */
	point G, /* The list of points to consider.
		    Sorted on coordinate 0. On ties sorted
		    on coordinate 1. On ties no further criteria. */
	char *B, /* List for storing the result, it should indicate, for
		    the points in the sub-problem if they are minima or
		    not. This array is indexed as G, meaning that B[i] = 1
		    implies that the i-th point in G is a minima.  As
		    mentioned only if i = R[j] for some j is the result
		    meaningfull. */
	/* Arguments related to sub-problem */
	int Ldim, /* Sub-problem dimension */
	int m, /* Number of points in sub-problem */
	int *R, /* List of points in sub-problem,
		   ordered as G. */
	char *killerR, /* 1 if R[i] point can kill others. */
	char *victimR /* 1 if R[i] point can be kille. */
	/* char *tB /\* List of booleans, with tB[i] = 1 iff R[i] is larger than */
	/* 	   median in the next dimension. *\/ */
	 )
{
  int j;

#ifndef NDEBUG
  assertOrder(m, R);
#endif /* NDEBUG */

  /* printf("minimaR call begin Ldim: %d .\n", Ldim); */
  /* printCall(Gdim, G, B, m, R, killerR, victimR); */

  if(1 < m){
    switch(Ldim){
    case 1: /* One dimentional solution */
      {
	int i=0;
	while(i < m && 0 == killerR[i])
	  i++;
	if(1 == killerR[i])
	  i++;
	while(i < m){
	  if(1 == victimR[i])
	    B[R[i]] = 0;
	  i++;
	}
      }
      break;

    case 2: /* 2D solution */
      {
	int smin = INT_MAX; /* Max on coord 2 */
	for(int i = 0; i < m; i++){
	  if(1 == killerR[i]
	     && G[R[i]*Gdim + 1] < smin)
	    smin = G[R[i]*Gdim + 1];
	  else if(1 == victimR[i]
		  && G[R[i]*Gdim + 1] >= smin)
	    B[R[i]] = 0;
	}
      }
      break;

    default :
      {
	int sm = m/2;
	char *lowtB = orderSelect(Gdim, G, Ldim, m, sm, R); /* Divide the array R */
#ifndef NDEBUG
	assertOrder(m, R);
#endif /* NDEBUG */

	if(NULL != lowtB){ /* Valid partition */
	  int *sR = malloc((m+2)*sizeof(int));
	  char *skillerR = malloc((m+2)*sizeof(char));
	  char *svictimR = malloc((m+2)*sizeof(char));
	  j = 0; /* index over sR */
	  for(int i = 0; i < m; i++){
	    sR[j] = R[i];
	    skillerR[j] = killerR[i];
	    svictimR[j] = victimR[i];
	    if(1 == lowtB[i]) j++;
	  }
	  minimaR(Gdim, G, B, Ldim, j, sR,
		  skillerR, svictimR);

	  j = 0; /* index over sR */
	  for(int i = 0; i < m; i++){
	    sR[j] = R[i];
	    skillerR[j] = killerR[i];
	    svictimR[j] = victimR[i];
	    if(0 == lowtB[i]) j++;
	  }
	  minimaR(Gdim, G, B, Ldim, j, sR,
		  skillerR, svictimR);

	  free(sR);
	  free(skillerR);
	  free(svictimR);

	  j = 0;
	  for(int i = 0; i < m; i++)
	    if(1 == B[R[i]]) /* If is minima */
	      j++;
	  int *redR = malloc((j+1)*sizeof(int));
	  char *redKillerR = malloc((j+1)*sizeof(char));
	  char *redVictimR = malloc((j+1)*sizeof(char));

	  /* Contract R and tB */
	  j = 0;
	  for(int i = 0; i < m; i++){
	    redR[j] = R[i];
	    redKillerR[j] = killerR[i] &&  lowtB[i];
	    redVictimR[j] = victimR[i] && !lowtB[i];
	    if(1 == B[R[i]]) /* If is minima */
	      j++;
	  }
	  minimaR(Gdim, G, B, Ldim-1, j, redR,
		  redKillerR, redVictimR);

	  free(redR);
	  free(redKillerR);
	  free(redVictimR);
	  free(lowtB);
	} else /* On invalid partition call for lower dim */
	  minimaR(Gdim, G, B, Ldim-1, m, R, killerR, victimR);
      }
      break;
    }
  }

  /* printf("minimaR call >>[END]<< Ldim: %d .\n", Ldim); */
  /* printCall(Gdim, G, B, m, R, killerR, victimR); */

#ifndef NDEBUG
  assertMinimaR(Gdim, G, B, Ldim, m, R, killerR, victimR);
#endif /* NDEBUG */
}


static int Global_Dimension;

static int
CoordCmp(const void *p1, const void *p2)
{
  int *i1 = (int *)p1;
  int *i2 = (int *)p2;

  int r = 0;
  int i = 0;

  while(0 == r && i < Global_Dimension){
    r = i1[i] - i2[i];
    i++;
  }

  return r;
}

/* Return array of booleans as result */
char *
minima(int Gdim, /* Dimension */
       point G, /* Points */
       int *pm /* Number of points */
       )
{
  int  m = *pm;

  Global_Dimension = Gdim;
  qsort(G, m, Gdim*sizeof(int), CoordCmp);

  /* Removing duplicates */
  int t = m;
  m = 0;
  for(int i = 1; i < t; i++){
    int dup = 1;
    for(int j = 0; dup && j < Gdim; j++)
      dup = G[m*Gdim+j] == G[i*Gdim+j];

    if(!dup){
      m++;
      if(m < i)
	memcpy(&G[m*Gdim], &G[i*Gdim], Gdim*sizeof(int));
    }
  }
  m++;

  int *R = malloc(m*sizeof(int));
  for(int i = 0; i < m; i++)
    R[i] = i;

/* Algorithm assumes all points are minima.
   Recursive calls prune the booleans. */
  char *B = calloc(m, sizeof(char));
  char *killerR = calloc(m, sizeof(char));
  char *victimR = calloc(m, sizeof(char));

  for(int i=0; i<m; i++){
    B[i] = 1;
    killerR[i] = 1;
    victimR[i] = 1;
  }

  minimaR(Gdim, G, B, Gdim, m, R, killerR, victimR);

  free(R);
  free(killerR);
  free(victimR);

  *pm = m;
  return B;
}

static int
sigmaSize(int dim, /* Dimension */
	  char** S /* The strings */
	  )
{
  char maxC = 'A';

  for(int i = 0; i < dim; i++){
    for(int j = 0; '\0' != S[i][j]; j++){
      if(S[i][j] > maxC)
	maxC=S[i][j];
    }
  }

  return 1 + maxC - 'A';
}

int
naiveMLCS(int dim, /* Dimension */
	  char** S /* The strings */
	  )
{
  int len[dim]; /* String dimensions */
  int r = -1;  /* The size of the LCS */
  int sigma = sigmaSize(dim, S);
  int ***next; /* next position with the same letter. */
  int nextP[sigma]; /* Next occ of a given letter */

  next = malloc(dim*sizeof(int**));

  for(int i = 0; i < dim; i++){
    len[i] = strlen(S[i]);
    if(0 == len[i]) /* If some string is empty */
      return 0;

    next[i] = malloc(sigma*sizeof(int*));

    for(int s = 0; s < sigma; s++){
      next[i][s] = calloc(len[i], sizeof(int));
      nextP[s] = len[i];
    }

    for(int j = len[i]; 0 < j;){
      j--;
      nextP[S[i][j]-'A'] = j;
      for(int s = 0; s < sigma; s++)
	next[i][s][j] = nextP[s];
    }
  }

  int m = 1; /* number of points in currentF */
  point currentF = malloc(dim*sizeof(int));
  point c;
  point p;

  for(int i = 0; i < dim; i++)
    currentF[i] = -1;

  while(0 < m){
    r++;

    int a = 2;
    int im = 0;
    point tempP = malloc(a*dim*sizeof(int));

    for(int i = 0; i<m; i++){
      c = &currentF[i*dim];
      for(int s = 0; s < sigma; s++){
	if(im == a){
	  a *= 2;
	  tempP = realloc(tempP, a*dim*sizeof(int));
	}
	p = &tempP[im*dim]; /* Generated point */

	int valid = 1;
	for(int i = 0; i < dim; i++){
	  if(c[i]+1 < len[i] &&
	     next[i][s][c[i]+1] < len[i]){
	    p[i] = next[i][s][c[i]+1];
	  } else {
	    valid = 0;
	    break; /* Go to the exterior cycle */
	  }
	}
	if(valid){
	  /* printf("@ %d creating point (%d, %d) for letter %c\n", */
	  /* 	   r, p[0], p[1], (char)('A'+s)); */
	  im++;
	}
      }
    }

    m = 0;
    if(0 < im){
      char *filter = minima(dim, tempP, &im);
      for(int i = 0; i < im; i++){
	if(m < i)
	  memcpy(&tempP[m*dim], &tempP[i*dim],
		 dim*sizeof(int));
	if(filter[i]) m++;
      }
      free(filter);
    }
    free(currentF);
    currentF = tempP;
  }

  free(currentF);

  for(int i = 0; i < dim; i++){
    for(int s = 0; s < sigma; s++)
      free(next[i][s]);
    free(next[i]);
  }
  free(next);

  return r;
}
