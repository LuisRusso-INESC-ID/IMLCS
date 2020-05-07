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
 *  \brief     Naive Multiple Longest Common Sub-sequence
 *  \details   This file implements the dynamic MLCS algorithm
 *  \author    Luís M. S. Russo
 *  \version   0.1.0-alpha
 *  \date      04-05-2020
 *  \copyright BSD 2-Clause License
 */

#include <string.h>
#include <stdlib.h>

/* Index to coordinate function */
static int /* Returns true if some coord is 0 */
idx2coord(int dim, /* Dimensions */
	  int* n,  /* Sizes */
	  int j,   /* Index */
	  int* c   /* Coordinates */
	  )
{
  int zc = 0; /* Zero in coordiante */
  for(int i = 0; i < dim; i++){
    c[i] = j % n[i];
    j /= n[i];
    zc = zc || (0 == c[i]);
  }

  return zc;
}

static int /* Returns index from coordinates */
coord2idx(int dim, /* Dimensions */
	  int* n,  /* Sizes */
	  int* c   /* Coordinates */
	  )
{
  int r = c[dim-1];
  for(int i = dim-2; 0 <= i; i--){
    r *= n[i]; /* Create space */
    r += c[i]; /* Include info */
  }

  return r;
}

int
naiveMLCS(int dim, /* Dimension */
	  char** S, /* The strings */
	  __attribute__((unused)) int sigma
	  )
{
  int n[dim]; /* String dimensions */
  int c[dim]; /* The coordinates of an index */
  int sc[dim]; /* The smaller coordinates of an index */
  int sz = 1; /* Size of the DP table */
  int* T;     /* Famous DP table */
  int r;

  /* The strings */
  for(int i = 0; i < dim; i++){
    n[i] = strlen(S[i]) + 1;
    if(1 == n[i]) /* If some string is empty */
      return 0;
    sz *= n[i];
  }

  T = malloc(sz*sizeof(int));

  for(int j = 0; j < sz ; j++){
    /* First decode coordinates of j */
    int zc = idx2coord(dim, n, j, c);

    T[j] = 0; /* Default with some zero coordinate */
    if(!zc){
      int lm = 1; /* Letter match */
      for(int i = 1; lm && i < dim; i++)
	lm = S[i][c[i]-1] == S[0][c[0]-1];
      if(lm){ /* Letter matched */
	for(int i = 0; i < dim; i++)
	  sc[i] = c[i]-1;
	T[j] = T[coord2idx(dim, n, sc)] + 1;
      } else {
	for(int i = 0; i < dim; i++){
	  c[i]--;
	  int jj = coord2idx(dim, n, c);
	  c[i]++;

	  if(T[j] < T[jj])
	    T[j] = T[jj];
	}
      }
    }
  }
  r = T[sz-1];
  free(T);
  return r;
}
