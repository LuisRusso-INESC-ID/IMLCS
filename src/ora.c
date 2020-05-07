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
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include "ora.h"

struct ora{
  int w; /* Weight, means number of elements stored,
	  Counting multiplicities. */
  int d; /* Number of dimensions */
  int o; /* Number of occupied positions */
  int a; /* Number of alloced points */
  int *A; /* array storing the points */
};

void
dsRnORA(ora R,
	FILE *f)
{
  fprintf(f, "m%x [label=%c0x%x | w: %d| d: %d| o: %d| a: %d| A: ",
	  (int)R, '"', (int)R,
	  R->w,
	  R->d,
	  R->o,
	  R->a
	  );

  for(int i = 0; i < R->o; i++){
    fprintf(f, "(");
    for(int j = 0; j <= R->d; j++){
      fprintf(f, "%d, ", R->A[i*(R->d+1)+j]);
    }
    fprintf(f, ") ");
  }
  fprintf(f, "%c]", '"');
}


/* Create a new ORA */
ora
allocORA(int dim /* Number of dimensions */
	 )
{
  ora R = calloc(1, sizeof(struct ora));

  R->d = dim;
  R->a = 2;
  R->A = malloc(R->a*(R->d+1)*sizeof(int));
  assert(0 < R->a && "Invalid alloc of ORA");

  return R;
}

/* Free the corresponding struct*/
void
freeORA(ora R
	)
{
  free(R->A);
#ifndef NDEBUG
  bzero(R, sizeof(struct ora));
#endif /* NDEBUG */
  free(R);
}

void
checkORA(ora R
	 )
{
  if(NULL != R){
    assert(0 >= R->w && "Weight fail on ORA");
    assert(0 < R->d && "Dimension fail on ORA");
    assert(10 >= R->d && "Dimension overflow on ORA");
    assert(0 <= R->o && "Use fail on ORA");
    assert(0 < R->a && "Allocation fail on ORA");
    assert(NULL != R->A && "Array fail on ORA");
  }
}

int
weightORA(ora R
	  )
{
  return R->w;
}

static int *
findORA(ora R, /* The orthogonal range tree */
	int* coords /* Point coordinates, in LSD order. */
	)
{
  int *found = NULL;

  for(int i=0; NULL == found && i < R->o; i++){
    found = &R->A[i*(R->d+1)];
    for(int j=0; NULL != found && j < R->d; j++){
      if(found[1+j] != coords[j])
	found = NULL;
    }
  }

  return found;
}

/* Range Queries */
int
countQORA(ora R, /* The orahogonal range tree */
          int* coords /* Point coordinates, in LSD order. */
          )
{
  int r = 0;

  for(int i=0; i < R->o; i++){
    int delta = R->A[i*(R->d+1)];
    for(int j=0; 0 < delta && j < R->d; j++){
      if(R->A[i*(R->d+1)+1+j] >= coords[j])
        delta = 0;
    }
    r += delta;
  }

  return r;
}

/* Important for application */
int
containsQORA(ora R, /* The orthogonal range tree */
             int* coords /* Point coordinates, in LSD order. */
             )
{
  return NULL != findORA(R, coords);
}

#include "point.h"

ora
buildORA(point *C, /* Array with point pointers */
	 int l, /* Left index, inclusive */
	 int r, /* Right index, inclusive */
	 int dim /* Current dimension */
	 )
{
  ora R = calloc(1, sizeof(struct ora));

  R->d = dim;
  R->o = r-l+1;
  R->a = 2*R->o; /* Leave some room */
  if(0 == R->a)
    R->a = 2;
  assert(0 < R->a && "Invalid alloc of ORA");
  R->A = malloc(R->a*(R->d+1)*sizeof(int));

  for(int i = l; i <= r; i++){
    memcpy(&R->A[(i-l)*(R->d+1)], C[i], (dim+1)*sizeof(int));
    R->w += C[i][0];
  }

  return R;
}

static int GLOBAL_dim;
static int *GLOBAL_Array;

static int
lastCCompare(const void *pp, const void *pq)
{
  int p = *(int *)pp;
  int q = *(int *)pq;

  return GLOBAL_Array[(1+p)*(GLOBAL_dim+1)-1] -
    GLOBAL_Array[(1+q)*(GLOBAL_dim+1)-1];
}

void
teleportORA(ora R, /* Tree node to traverse */
	    int *C, /* Array for storing points and multipls */
	    int *i, /* at the end it is the size */
	    int dim, /* Current dimension */
	    int mxdim, /* Maximum dimension */
	    point p,   /* Current common point coords from
			  mxdim to dim */
	    int cp /* Copy this node */
            )
{
  assert(NULL != R && "Teleporting empty struct");

  if(cp){ /* Means you want the nodes copied */
    /* First sort the nodes */
    if(1+dim == mxdim){ /* Needs sorting */
      int *P = malloc(R->o*sizeof(int));
      for(int j = 0; j < R->o; j++)
	P[j] = j;

      GLOBAL_dim = dim+1;
      GLOBAL_Array = R->A;
      qsort(P, R->o, sizeof(int), lastCCompare);

      for(int j = 0; j < R->o; j++){
	memcpy(&C[(*i)*(mxdim+1)], &(R->A[P[j]*(dim+2)]), (dim+2)*sizeof(int));
	(*i)++;
      }
      free(P);
    } else { /* Does not need sorting */
      for(int j = 0; j < R->o; j++){
	memcpy(&C[(*i)*(mxdim+1)], &(R->A[j*(dim+2)]), (dim+2)*sizeof(int));
	memcpy(&C[(*i)*(mxdim+1)+(dim+2)], &p[dim+1], (mxdim-dim-1)*sizeof(int));
	(*i)++;
      }
    }
  }
  freeORA(R);
}

/* Returns an array with the points that dominate
   the coordinates. */
void
collectORA(ora R, /* The orthogonal range array */
	   int *C,   /* Array for storing points. */
	   int *coords, /* Point coordinates, in LSD order. */
	   int maxdim,
	   point hp, /* Temporary coords */
	   int *n
           )
{
  point p; /* Current pointer */

  for(int i=0; i < R->o; i++){
    p = &R->A[i*(R->d+1)];
    for(int j=0; NULL != p && j < R->d; j++){
      if(p[1+j] <= coords[j])
	p = NULL;
    }
    if(NULL != p){
      for(int k = 0; k < p[0]; k++){
	assert(1 == p[0] && "Missed point rep.");
	memcpy(&(C[(*n)*maxdim]), &p[1], R->d*sizeof(int));
	memcpy(&(C[(*n)*maxdim+R->d]), &hp[R->d], (maxdim-R->d)*sizeof(int));
	(*n)++;
      }
    }
  }
}

void
dominatedCollectORA(ora R, /* The orthogonal range array */
                    int *C,   /* Array for storing points. */
                    int *coords, /* Point coordinates, in LSD order. */
                    int maxdim,
                    point hp, /* Temporary coords */
                    int *n
                    )
{
  point p; /* Current pointer */

  for(int i=0; i < R->o; i++){
    p = &R->A[i*(R->d+1)];
    for(int j=0; NULL != p && j < R->d; j++){
      if(p[1+j] >= coords[j])
	p = NULL;
    }
    if(NULL != p){
      for(int k = 0; k < p[0]; k++){
	assert(1 == p[0] && "Missed point rep.");
	memcpy(&(C[(*n)*maxdim]), &p[1], R->d*sizeof(int));
	memcpy(&(C[(*n)*maxdim+R->d]), &hp[R->d], (maxdim-R->d)*sizeof(int));
	(*n)++;
      }
    }
  }
}

void
rangeCollectORA(ora R, /* The orthogonal range array */
		int *C,   /* Array for storing points. */
		int *minCoords, /* Point coordinates, in LSD order. */
		int *maxCoords, /* Point coordinates, in LSD order. */
		int maxdim,
		point hp, /* Temporary coords */
		int *n
		)
{
  point p; /* Current pointer */

  for(int i=0; i < R->o; i++){
    p = &R->A[i*(R->d+1)];
    for(int j=0; NULL != p && j < R->d; j++){
      if(p[1+j] >= maxCoords[j])
	p = NULL;
      if(NULL != p && p[1+j] < minCoords[j])
	p = NULL;
    }
    if(NULL != p){
      for(int k = 0; k < p[0]; k++){
	assert(1 == p[0] && "Missed point rep.");
	memcpy(&(C[(*n)*maxdim]), &p[1], R->d*sizeof(int));
	memcpy(&(C[(*n)*maxdim+R->d]), &hp[R->d], (maxdim-R->d)*sizeof(int));
	(*n)++;
      }
    }
  }
}


void
insertORA(ora R, /* The orthogonal range tree */
          point p,
	  int w
          )
{
  assert(0 < R->a && "Invalid alloc of ORA");
  int *P = findORA(R, p); /* Position in array */
  assert(0 < R->a && "Invalid alloc of ORA");
  if(NULL != P){ /* Found point */
    P[0]+=w;
  } else { /* Adding point at the end. */
    assert(0 < R->a && "Invalid alloc of ORA");
    if(R->o == R->a){ /* Exausted memory */
      R->a *= 2;	    /* Double space */
      assert(0 < R->a && "Invalid alloc of ORA");
      R->A = realloc(R->A, R->a*(R->d+1)*sizeof(int));
    }

    R->A[R->o*(R->d+1)] = w; /* w points */
    memcpy(&R->A[R->o*(R->d+1)+1], p, R->d*sizeof(int));
    R->o++;
  }

  R->w+=w; /* Account for the new point */
}

int
deleteORA(ora R, /* The orthogonal range tree */
          point p
          )
{
  int deletedQ = 0;
  int *P = findORA(R, p); /* Position in array */
  assert(NULL != P && "Deleting non-existing point.");

  R->w--;
  P[0]--;
  if(0 == P[0]){ /* Pop element out of position */
    R->o--;
    if(&R->A[R->o*(R->d+1)] != P)
      memcpy(P, &R->A[R->o*(R->d+1)], (R->d+1)*sizeof(int));

    if(0 < R->o && 2*R->o < R->a){
      R->a /= 2;
      R->A = realloc(R->A, R->a*(R->d+1)*sizeof(int));
    }
  }
  assert(0 < R->a && "Invalid alloc of ORA");

  if(0 == R->w){
    free(R->A);
    bzero(R, sizeof(struct ora));
    free(R);
    deletedQ = 1;
  }

  return deletedQ;
}
