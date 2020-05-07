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
 *  \brief     Orthogonal Range Tree implementation
 *  \details   This structure supports the dynamic multiple longest common
sub-sequence algorithm. The binary search trees are weight balanced.
 *  \author    Luís M. S. Russo
 *  \version   0.1.0-alpha
 *  \date      17-12-2019
 *  \copyright BSD 2-Clause License
 */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <bsd/stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "point.h"
#include "ort.h"
#include "ora.h"

/**
The data structure is organized as follows:

 1) The BST is internal, meaning that nodes store actual values.
 2) We avoid key repetition and store a coordenate value only once.

*/

struct node{
  /* There are no virtual nodes, but only the
     v in the last dimension owns the point. */
  /* point v;    /\* Pointer to coordinates *\/ */
  int w;      /* Weight */
  int v;      /* coordinate in the tree */
  node left;  /* Sub-tree */
  node right; /* Sub-tree */
  node equal; /* Sub-tree, lower dimension */
  node lower; /* Sub-tree, lower dimension */
  /* node hook;  /\* Usually the parent, but it may */
  /*       	 also be the pointer to prev tree. *\/ */
};

struct ort{ /* Orthogonal range tree */
  int d;    /* Number of dimensions */
  node root; /* Tree root */
  int *cache; /* The cache for balancing and iterators */
  int n;	/* Number of points in tree. */
  int ca;	/* Cache size. */
};

/* Basic tree configuration */
static int mxd; /* Maximum allowed depth */
static double *T; /* A table containing balance constants */
static int cutoff; /* The cut-off value. If the number of elements
                      in a sub-tree is smaller then use an array. */

#ifndef NDEBUG
void
gdbBreak(void)
{
}
#endif /* NDEBUG */

void
configure(int mxn, /* Maximum number of elements in the tree. */
          int defmxd, /* Maximum depth of the tree */
          double alpha, /* Initial alpha value */
          int confcut
          )
{
  /* Verify alpha value */
  assert(3.0 / 4.0 <= alpha && alpha < 1 && "Invalid alpha value");

  mxd = defmxd;

  double epsilon; /* The other parameter */

  epsilon = log(mxn) / log(1/alpha);
  epsilon = log(epsilon) / log(mxd);

  assert(0 <= epsilon && "Invalid epsilon value");

  T = (double*)malloc((1+mxd)*sizeof(double));
  T[0] = 2; /* Simplifies some tests */
  int d; /* depth counter */
  for(d=0; d < mxd; d++){
    T[d+1]  = pow(d+1, epsilon);
    T[d+1] -= pow(d, epsilon);
    T[d+1]  = pow(alpha, T[d+1]);
    /* printf("T[%d] = %f\n", d+1, T[d+1]); */
  }

  cutoff = confcut;
}

void
deconfigure(void)
{
  free(T);
}

double
adjustCut(int param)
{
  double R=CUTOFF; /* The result */

  static int cutOff = 1; /* This is the original cutOff Value */
  static int cutAl = 0; /* size of alloced array */
  static double *cA = NULL;

  if(0 < param){
    if(cutAl <= param){
      cutAl *= 2;
      if(0 == cutAl)
	cutAl = 3;
      if(cutAl <= param)
	cutAl = param+1;
      cA = realloc(cA, cutAl*sizeof(double));

      double Choose[cutOff+1+cutAl][cutAl];
      for(int i = 0; i <= cutOff+cutAl; i++){
	for(int j = 0; j < cutAl; j++){
	  Choose[i][j] = 0;
	}
	Choose[i][0] = 1;
	if(0 < i){
	  for(int j = 1; j < cutAl; j++){
	    Choose[i][j] = Choose[i-1][j-1] + Choose[i-1][j];
	  }
	}
      }

      cA[0] = 1; /* Non-sense value */
      cA[1] = cutOff;
      cA[2] = cutOff;

      for(int i = 3; i < cutAl; i++){
	cA[i] = Choose[cutOff+i-1][i-1]/i;
      }
    }
    R = cA[param];
  }
  if(0 == param){
    free(cA);
  }
  if(0 > param){
    cutOff = -param;
  }

  return R;
}

int
weightORT(ort rt
	  )
{
  int r = 0;
  if(NULL != rt->root)
    r = abs(rt->root->w);

  return r;
}

/* Create a new ORT */
ort
allocORT(int dim /* Number of dimensions */
	 )
{
  ort r = calloc(1, sizeof(struct ort));

  r->d = dim;

  return r;
}

/* The proper weight of the node */
static int
properW(node t)
{
  int r = 0;

  if(NULL != t){
    r = t->w;
    if(r <= 0) /* Means ORT struct */
      r *= -1;
    else { /* Actual tree */
      if(NULL != t->left)
	r -= abs(t->left->w);
      if(NULL != t->right)
	r -= abs(t->right->w);
      }
  }

  return r;
}

/* This function should return the number of points that
are strictly dominated by the point, i.e., coords are < */

/* Recursive version. */

static int
countR(node t,     /* The orthogonal range tree */
       int* coords, /* Point coordinates, in LSD order. */
       int dim    /* Current dimension */
       )
{
  int r = 0;

  while(NULL != t && 0 < t->w){
    if(t->v < coords[dim]){
      if(0 == dim){
	r += t->w;
	if(NULL != t->right)
	  r -= abs(t->right->w);
      } else { /* 0 < dim */
	if(NULL != t->left){
	  if(0 < t->left->w)
	    r += countR(t->left->lower, coords, dim-1);
	  else
	    r += countR(t->left, coords, dim-1);
	}
	r += countR(t->equal, coords, dim-1);
      }
      t = t->right;
    } else
      t = t->left;
  }

  /* Handling cut-off structure */
  if(NULL != t && 0 >= t->w){
    t->w *= -1;
    r += countQORA((ora)t, coords);
    t->w *= -1;
  }

  return r;
}

int
countQ(ort rt, /* The orthogonal range tree */
       int* coords /* Point coordinates, in LSD order. */
       )
{
  return countR(rt->root, coords, rt->d-1);
}

int
containsQ(ort rt, /* The orthogonal range tree */
	  int* coords /* Point coordinates, in LSD order. */
	  )
{
  /* dotShow(rt->root); */
  node t = rt->root;
  int dim = rt->d-1;
  int r = 1;

  while(0 <= dim && r){
    if(NULL != t && 0 >= t->w)
      return containsQORA((ora)t, coords);
    while(NULL != t && t->v != coords[dim]){
      if(t->v < coords[dim])
	t = t->right;
      else
	t = t->left;
      if(NULL != t && 0 >= t->w)
	return containsQORA((ora)t, coords);
    }

    if(NULL == t || 0 == properW(t))
      r = 0;
    else { /* Coordinates matched */
      dim--;
      t = t->equal;
    }
  }

  return r;
}

static void
dsRn(node t,
     FILE *f)
{
  if(NULL != t && 0 < t->w){
    fprintf(f, "m%x [label=%c0x%x| w: %d| v: %d| <l> left|<r> right|<e> equal|<o> lower%c];\n",
	    (int)t, '"', (int)t, t->w, t->v,'"');
    dsRn(t->left, f);
    dsRn(t->right, f);
    dsRn(t->equal, f);
    dsRn(t->lower, f);
 } else if(NULL != t && 0 >= t->w)
    dsRnORA((ora)t, f);
}

static void dsRp(node t,
		 FILE *f)
{
  if(NULL != t && 0 < t->w){
    if(NULL !=t->left){
      fprintf(f, "m%x:l -> m%x;\n", (int)t, (int)t->left);
      dsRp(t->left, f);
    }
    if(NULL !=t->right){
      fprintf(f, "m%x:r -> m%x;\n", (int)t, (int)t->right);
      dsRp(t->right, f);
    }
    if(NULL !=t->equal){
      fprintf(f, "m%x:e -> m%x;\n", (int)t, (int)t->equal);
      dsRp(t->equal, f);
    }
    if(NULL !=t->lower){
      fprintf(f, "m%x:o -> m%x;\n", (int)t, (int)t->lower);
      dsRp(t->lower, f);
    }
  }
}

void dotShow(node t
	     )
{
  static int fn = 1; /* File number */
  static char fname[40];
  FILE *f;

  if(NULL != t && 0 < t->w){

    sprintf(fname, "OrthRangeT%.4d.dot", fn++);
    f = fopen(fname, "w");

    fprintf(f, "digraph g { ");
    fprintf(f, "rankdir=LR; ");
    fprintf(f, "node[shape=record]; \n");

    dsRn(t, f); /* Recursive call for names */
    dsRp(t, f); /* Recursive call for links */
    fprintf(f, "}\n");

    fclose(f);
  }
}

#ifndef NDEBUG
static void
checkTree(node t, int dim)
{
  if(NULL != t){
    if(0 >= t->w){
      checkORA((ora)t);
    } else { /* 0 < t->w */
      if(0 < dim){
	if(0 < properW(t)){
	  assert(NULL != t->equal && "Check lower dim.");
	  assert(properW(t) == abs(t->equal->w) && "Equal weight check.");
	  checkTree(t->equal, dim-1);
	}
	assert(NULL != t->lower && "Check lower dim.");
	assert(t->w == abs(t->lower->w) && "Lower weight check.");
	checkTree(t->lower, dim-1);
      } else {
	assert(NULL == t->equal && "Check lower dim.");
	assert(NULL == t->lower && "Check lower dim.");
      }

      if(NULL != t->left){
	/* Incomplete verification but ok */
	if(0 < t->left->w)
	  assert(t->left->v < t->v && "Order verification");
	checkTree(t->left, dim);
      }

      if(NULL != t->right){
	/* Incomplete verification but ok */
	if(0 < t->right->w)
	  assert(t->right->v > t->v && "Order verification");
	checkTree(t->right, dim);
      }
    }
  }
}
#endif /* NDEBUG */

/* Recursive version. */

static void
collectR(node t,     /* The orthogonal range tree. */
	 int *C,   /* Array for storing points. */
	 int *coords, /* Point coordinates, in LSD order. */
	 int dim,    /* Current dimension */
	 int maxdim,
	 point hp, /* Temporary coords */
	 int *n
	 )
{
  if(NULL != t){
    if(0 >= t->w){ /* Means ORA struct */
      t->w *= -1;
      collectORA((ora)t, C, coords, maxdim, hp, n);
      t->w *= -1;
    } else { /* Means inside tree. */
      if(t->v > coords[dim]){
	collectR(t->left, C, coords, dim, maxdim, hp, n);
	hp[dim] = t->v;
	if(0 < dim)
	  collectR(t->equal, C, coords, dim-1, maxdim, hp, n);
	else if(0 < properW(t)){
	  assert(1 == properW(t) && "Missed point rep.");
	  memcpy(&C[((*n)++)*maxdim], hp, maxdim*sizeof(int));
	}
      }
      collectR(t->right, C, coords, dim, maxdim, hp, n);
    }
  }
}

point *
collect(ort rt, /* The orthogonal range tree */
	int* coords, /* Point coordinate */
	int* n /* Number of points */
	)
{
  int hp[rt->d];   /* temp memory */
  point *R = NULL; /* The result */
  *n = 0;
  collectR(rt->root, rt->cache, coords,
	   rt->d-1, rt->d, hp, n);

  R = malloc(*n*sizeof(point));
  for(int i = 0; i<*n; i++){
    R[i] = malloc(rt->d*sizeof(int));
    memcpy(R[i], &(rt->cache[i*rt->d]), rt->d*sizeof(int));
  }

  return R;
}

static void
dominatedCollectR(node t,     /* The orthogonal range tree. */
		  int *C,   /* Array for storing points. */
		  int *coords, /* Point coordinates, in LSD order. */
		  int dim,    /* Current dimension */
		  int maxdim,
		  point hp, /* Temporary coords */
		  int *n
		  )
{
  if(NULL != t){
    if(0 >= t->w){ /* Means ORA struct */
      t->w *= -1;
      dominatedCollectORA((ora)t, C, coords, maxdim, hp, n);
      t->w *= -1;
    } else { /* Means inside tree. */
      if(t->v < coords[dim]){
	dominatedCollectR(t->right, C, coords, dim, maxdim, hp, n);
	hp[dim] = t->v;
	if(0 < dim)
	  dominatedCollectR(t->equal, C, coords, dim-1, maxdim, hp, n);
	else if(0 < properW(t)){
	  assert(1 == properW(t) && "Missed point rep.");
	  memcpy(&C[((*n)++)*maxdim], hp, maxdim*sizeof(int));
	}
      }
      dominatedCollectR(t->left, C, coords, dim, maxdim, hp, n);
    }
  }
}

point *
dominatedCollect(ort rt, /* The orthogonal range tree */
		 int* coords, /* Point coordinate */
		 int* n /* Number of points */
		 )
{
  int hp[rt->d];   /* temp memory */
  point *R = NULL; /* The result */
  *n = 0;
  dominatedCollectR(rt->root, rt->cache, coords,
	   rt->d-1, rt->d, hp, n);

  if(0 < *n){
    R = malloc(*n*sizeof(point));
    for(int i = 0; i<*n; i++){
      R[i] = malloc(rt->d*sizeof(int));
      memcpy(R[i], &(rt->cache[i*rt->d]), rt->d*sizeof(int));
    }
  }

  return R;
}

static void
rangeCollectR(node t,     /* The orthogonal range tree. */
	      int *C,   /* Array for storing points. */
	      int *minCoords, /* Point coordinates, in LSD order. */
	      int *maxCoords, /* Point coordinates, in LSD order. */
	      int dim,    /* Current dimension */
	      int maxdim,
	      point hp, /* Temporary coords */
	      int *n
	      )
{
  if(NULL != t){
    if(0 >= t->w){ /* Means ORA struct */
      t->w *= -1;
      rangeCollectORA((ora)t, C, minCoords, maxCoords,
		      maxdim, hp, n);
      t->w *= -1;
    } else { /* Means inside tree. */
      if(t->v < maxCoords[dim]){
	rangeCollectR(t->right, C, minCoords, maxCoords,
		      dim, maxdim, hp, n);
	if(minCoords[dim] <= t->v){
	  hp[dim] = t->v;
	  if(0 < dim)
	    rangeCollectR(t->equal, C, minCoords, maxCoords,
			  dim-1, maxdim, hp, n);
	  else if(0 < properW(t)){
	    assert(1 == properW(t) && "Missed point rep.");
	    memcpy(&C[((*n)++)*maxdim], hp, maxdim*sizeof(int));
	  }
	}
      }
      if(minCoords[dim] <= t->v)
	rangeCollectR(t->left, C,  minCoords, maxCoords,
		      dim, maxdim, hp, n);
      
    }
  }
}

point *
rangeCollect(ort rt, /* The orthogonal range tree */
	     int* minCoords, /* Point coordinate */
	     int* maxCoords, /* Point coordinate */
	     int* n /* Number of points */
	     )
{
  int hp[rt->d];   /* temp memory */
  point *R = NULL; /* The result */
  *n = 0;
  rangeCollectR(rt->root, rt->cache, minCoords, maxCoords,
		rt->d-1, rt->d, hp, n);

  if(0 < *n){
    R = malloc(*n*sizeof(point));
    for(int i = 0; i<*n; i++){
      R[i] = malloc(rt->d*sizeof(int));
      memcpy(R[i], &(rt->cache[i*rt->d]), rt->d*sizeof(int));
    }
  }

  return R;
}

/* Rebalance */
/* Swap elements i and j in array. */
static void
swap(point* A,
     int i,
     int j
     )
{
  point t;
  t = A[i];
  A[i] = A[j];
  A[j] = t;
}

/* QuickSort 3-way partition, assumes the pivot is at *r. */
/* See Sedgewick Algorithms 4th ed, pag 299 */
static void
partition(point *A, /* Array with point pointers */
          int *l, /* Left index, inclusive */
          int *r, /* Right index, inclusive */
          int p,  /* The dividing pivot */
          int dim /* Current dimension */
          )
{ /* Use *l for exclusive end of less than elements */
  /* Use *r for exclusive end of greater than elements */
  int i = *l; /* simple counter */
  while(i <= *r){
    if(getCoord(A[i], dim) > p)
      swap(A, i, (*r)--);
    else {
      if(getCoord(A[i], dim) < p)
        swap(A, i, (*l)++);
      i++;
    }
  }

  /* At the end [*l, *r] is the interval of pivot positions */
#ifndef NDEBUG
  i = *l;
  while(i < *r){
    assert(getCoord(A[i], dim) == p && "Failed pivot interval in partition.");
    i++;
  }
#endif /* NDEBUG */
}

/* Sum up multiplicities */
static int
msum(point *A,
     int l,
     int r
     )
{
  int sum = 0;
  while(l <= r){
    sum += A[l][0];
    l++;
  }
  return sum;
}

static int *
accSum(point *A, /* The array with point multiplicities */
       int n
       )
{
  int *R = malloc((n+1)*sizeof(int));
  R[0] = 0;
  for(int i = 1; i <= n; i++){
    R[i]  = R[i-1];
    R[i] += A[i-1][0];
  }

  return R;
}

/* int */
/* myRandom(int ubound) */
/* { */
/*   return ubound/4; */
/*   return ubound/2; /\* Good for debugging. *\/ */
/* } */

/* Returns the median inside the interval [*l, *r]  and the
array A properly partitioned. */
static void
median(point *A, /* Array with points and mults */
       int *l, /* Left index, inclusive */
       int *r, /* Right index, inclusive */
       int dim, /* Current dimension */
       int *Acc /* Array of accumulated values */
       )
{
  if(NULL == Acc){
    int lsum = 0;
    int rsum = msum(A, *l, *r);
    assert(0<rsum && "Failed median call.");
    int mid = rsum/2;

    while(*r - *l + 1 > 2){ /* The interval has at least 3 elements */
      /* Start by selecting median of 3 */
      swap(A, *r-0, *l+arc4random_uniform(*r-0-*l+1));
      swap(A, *r-1, *l+arc4random_uniform(*r-1-*l+1));
      swap(A, *r-2, *l+arc4random_uniform(*r-2-*l+1));

      if(pointcmp(A[*r-2], A[*r-1], dim) > 0)
	swap(A, *r-2, *r-1);
      if(pointcmp(A[*r-1], A[*r-0], dim) > 0)
	swap(A, *r-1, *r-0);
      if(pointcmp(A[*r-2], A[*r-1], dim) > 0)
	swap(A, *r-2, *r-1);
      /* Median of 3 is at *r-1 */

      int pl = *l;
      int pr = *r;
      partition(A, &pl, &pr, getCoord(A[*r-1], dim), dim);
      int addl = msum(A, *l, pl-1);
      int subr = msum(A, pr+1, *r);

      if(lsum + addl <= mid && mid <= rsum - subr){
	*l = pl;
	*r = pr;
	break; /* End cycle */
      } else {
	if(rsum - subr < mid){
	  *l = pr+1;
	  lsum = rsum - subr;
	} else if(mid < lsum + addl){
	  *r = pl-1;
	  rsum = lsum + addl;
	}
      }
    }

    /* These are for lists of size 2 */
    if(pointcmp(A[*l], A[*r], dim) > 0)
      swap(A, *l, *r);
    if(pointcmp(A[*l], A[*r], dim) != 0){
      if(lsum + A[*l][0] < mid)
	*l = *r;
      else
	*r = *l;
    }
  } else { /* No need for partition call, A is sorted */
    int ol = *l; /* Original l */
    int or = *r; /* Original r */

    /* The invariant is *l is included
       and *r is excluded */
    (*r)++;
    int mid = Acc[*r]/2;

    while(*l + 1 < *r){ /* Binary search */
      int m = (*l + *r)/2;
      if(Acc[m+1] <= mid)
	*l = m;
      else
	*r = m;
    } /* Now *l contains the target */
    *r = *l;
    while(ol < *l &&
	  getCoord(A[*l], dim) == getCoord(A[*l-1], dim))
      (*l)--;
    while(*r < or &&
	  getCoord(A[*l], dim) == getCoord(A[*r+1], dim))
      (*r)++;
  }

  assert(getCoord(A[*l], dim) ==
	 getCoord(A[*r], dim) &&
	 "Failed median value.");
}

static node
buildBalanced(point *C, /* Array with point pointers */
              int l, /* Left index, inclusive */
              int r, /* Right index, inclusive */
              int dim, /* Current dimension */
	      int mxdim, /* Max dimension */
	      int *Acc	 /* Array of accumulated values */
              )
{
  node root;

  /* if(r - l + 1 <= cutoff){ */
  if(r - l + 1 <= adjustCut(dim+1)){
    root = (node)buildORA(C, l, r, dim+1);
    root->w *= -1; /* Signal that it is an ora struct */
  } else {

    int ml = l; /* Median left position. */
    int mr = r; /* Median right position. */
    median(C, &ml, &mr, dim+1, Acc); /* dim + 1 offset mult */

    root = calloc(1, sizeof(struct node));
    root->v = getCoord(C[ml], dim+1); /* Node representant */
    root->w = 0;
    for(int i = l; i <= r; i++)
      root->w += C[i][0];

    /* WARNING: Do not change recursion order, needs extra partitions. */
    if(l < ml)
      root->left = buildBalanced(C, l, ml-1, dim, mxdim, Acc);
    if(mr < r)
      root->right = buildBalanced(C, mr+1, r, dim, mxdim, Acc);
    if(0 < dim){
      /* The equal points. */
      root->equal = buildBalanced(C, ml, mr, dim-1, mxdim, NULL);
      root->lower = buildBalanced(C, l, r, dim-1, mxdim, NULL);
    }

#ifndef NDEBUG
    /* dotShow(root); */
    checkTree(root, dim);
#endif /* NDEBUG */
  }

  /* dotShow(root); */
  return root;
}

/* Copies the points in sub-tree to cache and deletes it. */
static void
teleportR(node t, /* Tree node to traverse */
	  int *C, /* Array for storing points and multipls */
	  int *i, /* at the end it is the size */
	  int dim, /* Current dimension */
	  int mxdim, /* Maximum dimension */
	  point p,   /* Current common point coords from
			mxdim to dim */
	  int cp /* Copy this node */
	  )
{
  if(0 >= t->w){ /* In ORA struct */
    t->w *= -1;
    teleportORA((ora)t, C, i, dim, mxdim, p, cp);
    /* No need to sign t->w. It is free by now. */
  } else { /* In tree */
    int wt = properW(t);

    if(NULL != t->left)
      teleportR(t->left, C, i, dim, mxdim, p, cp);

    p[dim] = t->v;
    if(0 < dim){
      if(NULL != t->equal) /* When there are equal points */
	teleportR(t->equal, C, i, dim-1, mxdim, p, cp);
      teleportR(t->lower, C, i, dim-1, mxdim, p, 0); /* Call for free */
    } else if(cp && 0 < wt){
      memcpy(&C[(*i)*(mxdim+1)+1], p, mxdim*sizeof(int));
      C[(*i)*(mxdim+1)] = wt;
      (*i)++;
    }

    if(NULL != t->right)
      teleportR(t->right, C, i, dim, mxdim, p, cp);

#ifdef NDEBUG
    bzero(t, sizeof(struct node));
#endif /* NDEBUG */
    free(t); /* Delete over here. Avoids code duplication. */
  }
}

static point *
teleport(node t, /* Tree node to traverse */
	 int *C, /* Array for storing points and multipls */
	 int *i, /* at the end it is the size */
	 int dim, /* Current dimension */
	 int cp /* Copy this node */
	 )
{
  assert(NULL != t && "Calling teleport on empty tree.");
#ifndef NDEBUG
  int oi = *i; /* Original i value */
#endif /* NDEBUG */

  point *R = NULL;
  int p[dim+1]; /* Current dimension rewrite */

  teleportR(t, C, i, dim, dim+1, p, cp);
  if(*i > 0 && cp){
    R = (point *)malloc(*i*sizeof(point));
    for(int j=0; j<*i; j++)
      R[j] = (point)&C[j*(dim+1+1)];
  }

#ifndef NDEBUG
  for(; oi+1 < *i; oi++)
    assert(pointcmp( &(R[oi][1]), &(R[oi+1][1]), dim) <= 0
	   && "Messed up teleport order.");
#endif /* NDEBUG */

  return R;
}

void
freeORT(ort rt)
{
  /* Free Points */
  int n = 0;

  if(NULL != rt->root)
    teleport(rt->root, NULL, &n, rt->d-1, 0);

  rt->d = 0;
  rt->root = NULL;
  free(rt->cache);
  rt->cache = NULL;
  rt->n = 0;
  rt->ca = 0;
  free(rt);
}

/* A function to insert a point into the ort. */
/* Recursive version. */

static node *
selectChild(node pt,
	    point p,
	    int dim
	    )
{
  node *sel = NULL;
  if(NULL != pt && 0 < pt->w){ /* Tree exists and is not ORA */
    sel = &(pt->equal);
    if(getCoord(p, dim) < pt->v)
      sel = &(pt->left);
    if(pt->v < getCoord(p, dim))
      sel = &(pt->right);
  }

  return sel;
}

static void
insertR(node *t, /* Pointer to orthogonal range tree */
        int depth, /* Checking node balance. */
        point p, /* The point info */
        int *C, /* Cache for points */
        int dim /* Current dimension */
        )
{
  int emptyt = (*t == NULL); /* Inital t was empty */
  int equalFound = 0;
  int edgeFound = 0;
  int balanced = 1;
  node *sel = selectChild(*t, p, dim);

  if(emptyt){
    /* Do not change balanced */
    /* Do not change equal found */
    sel = t;
  }else{
    edgeFound = (NULL == *sel);
    equalFound = (&((*t)->equal) == sel);
    if(!edgeFound)
      /* balanced = (1+(*sel)->w <= T[depth]*(1+(*t)->w)); */
      balanced = (1+(*sel)->w < (BALANCE_FACTOR)*(1+(*t)->w));

    while(!edgeFound && !equalFound && balanced){
      if(0 < dim){
	insertR(&((*t)->lower), depth,
                p, C, dim-1); /* For sure in this sub-tree */
	/* (*t)->lower->hook = *t; */
      }
      (*t)->w++; /* Add the new point to the total weight */
      t = sel; /* Execute move */
      depth++;
      sel = selectChild(*t, p, dim);

      edgeFound = (NULL == *sel);
      equalFound = (&((*t)->equal) == sel);
      if(!edgeFound)
	/* balanced = (1+(*sel)->w <= T[depth]*(1+(*t)->w)); */
	balanced = (1+(*sel)->w < (BALANCE_FACTOR)*(1+(*t)->w));
    }
  }

  if(!balanced){
    int n = 1;
    point *R = teleport(*t, C, &n, dim, 1);
    int i = 1; /* Remove possible duplicates */
    while(i < n &&
	  !pointEquals(&(R[i][1]), p, dim+1))
      i++;

    if(i < n){	 /* Found duplicate */
      R[i][0]++; /* Increase weight */
      i = 1;
    } else { /* Insert because it is new */
      memcpy(&C[1], p, (dim+1)*sizeof(int));
      C[0] = 1; /* One point to copy */
      i = 0;
      while(i+1 < n &&
	    pointFullcmp(&(R[i][1]), &(R[i+1][1]), dim) > 0){
	swap(R, i, i+1);
	i++;
      }
      i = 0;
    }
    int *Acc = accSum(&R[i], n-i);
    *t = buildBalanced(R, i, n-1, dim, dim+1, &(Acc[-i]));
    free(Acc);
     if(NULL != R)
       free(R);
  } else { /* Tree is balanced */
    if(!equalFound){ /* Alloc and move to new node */
      if(!emptyt){
	(*t)->w++;
	if(0 < dim){
	  insertR(&((*t)->lower), depth,
                  p, C, dim-1); /* For sure in this sub-tree */
	  /* (*t)->lower->hook = *t; */
	}
      }
      *sel = calloc(1, sizeof(struct node));
      (*sel)->v = getCoord(p, dim);
      /* (*sel)->hook = *t; */
      t = sel;
      depth++;
    } /* Now equalFound = 1 */
    (*t)->w++; /* Account for new point */
    if(0 < dim){ /* Insert into lower dim. */
      insertR(&((*t)->equal), depth, p, C, dim-1);
      /* (*t)->equal->hook = *t; */
      insertR(&((*t)->lower), depth, p, C, dim-1);
      /* (*t)->lower->hook = *t; */
    }
  }
}

void
insertORAWrap(
	      node *t,
	      point p,
	      int *C, /* Cache for points */
	      int dim /* Current dimension */
	      )
{
  if(NULL == *t)
    *t = (node)allocORA(dim+1);

  (*t)->w *= -1;
  insertORA((ora)(*t), p, 1);
  (*t)->w *= -1;

  if(-(*t)->w >= 2*adjustCut(dim+1)){
    int n = 0;
    point *R = teleport(*t, C, &n, dim, 1);

    int *Acc = accSum(R, n);
    *t = buildBalanced(R, 0, n-1, dim, dim+1, &(Acc[0]));
    free(Acc);
    if(NULL != R)
      free(R);
  }
}

static void
insertCR(node *t, /* Pointer to orthogonal range tree */
        int depth, /* Checking node balance. */
        point p, /* The point info */
        int *C, /* Cache for points */
        int dim /* Current dimension */
        )
{
  if(NULL == *t || 0 >= (*t)->w)
    insertORAWrap(t, p, C, dim);   /* End function on else */
  else {
    int equalFound = 0;
    int edgeFound = 0; /* Now means cutoff */
    int balanced = 1;
    node *sel = selectChild(*t, p, dim);

    equalFound = (&((*t)->equal) == sel);
    if(!equalFound){
      edgeFound = NULL == *sel;
      if(!edgeFound) /* Might happen on equal */
	edgeFound = 0 >= (*sel)->w;
    }
    if(NULL != *sel)
      balanced = (1+abs((*sel)->w) < (BALANCE_FACTOR)*(1+(*t)->w));

    while(!edgeFound && !equalFound && balanced){
      if(0 < dim){
	insertCR(&((*t)->lower), depth,
		 p, C, dim-1); /* For sure in this sub-tree */
      }
      (*t)->w++; /* Add the new point to the total weight */
      t = sel; /* Execute move */
      depth++;
      sel = selectChild(*t, p, dim);

      equalFound = (&((*t)->equal) == sel);
      if(!equalFound){
	edgeFound = NULL == *sel;
	if(!edgeFound) /* Might happen on equal */
	  edgeFound = 0 >= (*sel)->w;
      }
      if(NULL != *sel)
	balanced = (1+abs((*sel)->w) < (BALANCE_FACTOR)*(1+(*t)->w));
    }

    if(!balanced){
      int n = 1;
      point *R = teleport(*t, C, &n, dim, 1);
      int i = 1; /* Remove possible duplicates */
      while(i < n &&
	    !pointEquals(&(R[i][1]), p, dim+1))
	i++;

      if(i < n){	 /* Found duplicate */
	R[i][0]++; /* Increase weight */
	i = 1;
      } else { /* Insert because it is new */
	memcpy(&C[1], p, (dim+1)*sizeof(int));
	C[0] = 1; /* One point to copy */
	i = 0;
	while(i+1 < n &&
	      pointFullcmp(&(R[i][1]), &(R[i+1][1]), dim) > 0){
	  swap(R, i, i+1);
	  i++;
	}
	i = 0;
      }
      int *Acc = accSum(&R[i], n-i);
      *t = buildBalanced(R, i, n-1, dim, dim+1, &(Acc[-i]));
      free(Acc);
      if(NULL != R)
	free(R);
    } else { /* Tree is balanced */
      if(edgeFound)
	insertORAWrap(sel, p, C, dim);
      (*t)->w++;
      if(0 < dim){ /* Insert into lower dim. */
	if(equalFound)
	  insertCR(&((*t)->equal), depth, p, C, dim-1);
	insertCR(&((*t)->lower), depth, p, C, dim-1);
      }
    }
  }
}

void
insert(ort rt,
       point p
       )
{
  /* printf(" >I"); */
  /* for(int i = 0; i < rt->d ; i++) */
  /*   printf(" %d,", p[i]); */
  /* printf("\n"); */
#ifndef NDEBUG
  assert(!containsQ(rt, p) && "Inserting point that is in the tree");
  /* dotShow(rt->root); */
  checkTree(rt->root, rt->d-1);
#endif /* NDEBUG */

  if(rt->n+1 >= rt->ca){
    if(0 == rt->ca)
      rt->ca = 1;
    rt->ca *= 2;
    if(NULL != rt->cache)
      free(rt->cache);

    rt->cache = malloc(rt->ca*(1+rt->d)*sizeof(int));
  }
  rt->n++;

  if(0 >= adjustCut(1)) /* No cut-off struct */
    insertR(&(rt->root), 1, p, rt->cache, rt->d-1);
  else /* With cut-off struct */
    insertCR(&(rt->root), 1, p, rt->cache, rt->d-1);

#ifndef NDEBUG
  /* dotShow(rt->root); */
  assert(0 < rt->n && "Missed point count on insert.");
  assert(rt->n == abs(rt->root->w) && "Missed point insert.");
  if(!containsQ(rt, p)){
    gdbBreak();
    dotShow(rt->root);
  }
  assert(containsQ(rt, p) && "Insert failed");
  checkTree(rt->root, rt->d-1);
#endif /* NDEBUG */
}

/* A function to delete a point from the ort. Note the function
 assumes the point is inside the tree. */

static void
deleteR(node *t, /* Pointer to orthogonal range tree */
        int depth, /* Checking node balance. */
        point p, /* The point info */
        int *C, /* Cache for points */
        int dim /* Current dimension */
        )
{
  assert(*t != NULL && "Deleting on empty tree");

  node *sel = selectChild(*t, p, dim); /* Select next */
  int nodeFound = (&((*t)->equal) == sel);
  /* Deleting on a found node keeps the balance, so that is
     the default. */
  int balanced = 1 < (*t)->w; /* Prune trees with 1 element */
  if(!nodeFound && balanced)
    balanced = (-1+(*sel)->w < (BALANCE_FACTOR)*(-1+(*t)->w));

  /* End cycle if: */
  while(!nodeFound && balanced){
    if(0 < dim){
      deleteR(&((*t)->lower), depth,
              p, C, dim-1); /* For sure in this sub-tree */
    }
    (*t)->w--; /* remove the point from the total weight */
    t = sel; /* Execute move */
    depth++;

    sel = selectChild(*t, p, dim);
    nodeFound = (&((*t)->equal) == sel);
    balanced = 1 < (*t)->w; /* Prune trees with 1 element */
    if(!nodeFound && balanced)
      balanced = (-1+(*sel)->w < (BALANCE_FACTOR)*(-1+(*t)->w));
  }

  if(!balanced){
    int n = 0; /* Empty list */
    point *R = teleport(*t, C, &n, dim, 1);
    int i; /* Make sure p gets removed */
    /* This asserts point in tree, will segfault */
    for(i = 0;
	!pointEquals(&(R[i][1]), p, dim+1);
	i++) ;
    if(pointEquals(&(R[i][1]), p, dim+1)){
      R[i][0]--; /* Decrease weight */
      if(0 == R[i][0]){ /* Remove on empty */
	for(;i+1<n; i++)
	  R[i] = R[i+1];
	n--;
      }
    }
    *t = NULL;
    if(0 < n){
      int *Acc = accSum(R, n);
      *t = buildBalanced(R, 0, n-1, dim, dim+1, Acc);
      free(Acc);
    }
    if(NULL != R)
      free(R);
  } else { /* Tree is balanced */
    assert(nodeFound && "Error p not in the tree.");
    /* Now deal nodeFound */
    (*t)->w--; /* Remove point */
    if(0 < dim){ /* Recursive into lower dim. */
      deleteR(&((*t)->equal), depth, p, C, dim-1);
      deleteR(&((*t)->lower), depth, p, C, dim-1);
    }
    if(0 == (*t)->w){ /* Free the node */
      bzero(*t, sizeof(struct node));
      free(*t);
      *t = NULL;
    }
  }
}

static void
deleteCR(node *t, /* Pointer to orthogonal range tree */
        int depth, /* Checking node balance. */
        point p, /* The point info */
        int *C, /* Cache for points */
        int dim /* Current dimension */
        )
{
  assert(*t != NULL && "Deleting on empty tree");

  if(0 >= (*t)->w){
    (*t)->w *= -1;
    if(deleteORA((ora)*t, p))
      *t = NULL;
    else{
      assert(0 < (*t)->w && "Deleting weight to 0");
      (*t)->w *= -1;
    }

  } else {
    node *sel = selectChild(*t, p, dim); /* Select next */
    int edgeFound = 0;
    if(NULL != *sel) /* Happens on equal @ dim=0 */
      edgeFound = 0 >= (*sel)->w;
    int nodeFound = (&((*t)->equal) == sel);
    /* Deleting on a found node keeps the balance, so that is
       the default. */
    int balanced = 1; /* Prune trees with 1 element */
    if(NULL != *sel)
      balanced = (-1+abs((*sel)->w) < (BALANCE_FACTOR)*(-1+(*t)->w));
    if(adjustCut(dim+1)+1 == (*t)->w)
      balanced = 0;

    /* End cycle if: */
    while(!nodeFound && !edgeFound && balanced){
      if(0 < dim){
	deleteCR(&((*t)->lower), depth,
		 p, C, dim-1); /* For sure in this sub-tree */
      }
      (*t)->w--; /* remove the point from the total weight */
      assert(0 < (*t)->w && "Deleting weight to 0");
      t = sel; /* Execute move */
      depth++;

      sel = selectChild(*t, p, dim);
      if(NULL != *sel) /* Happens on equal @ dim=0 */
	edgeFound = 0 >= (*sel)->w;
      nodeFound = (&((*t)->equal) == sel);
      if(NULL != *sel)
	balanced = (-1+abs((*sel)->w) < (BALANCE_FACTOR)*(-1+(*t)->w));
      if(adjustCut(dim+1)+1 == (*t)->w)
	balanced = 0;
    }

    if(!balanced){
      int n = 0; /* Empty list */
      point *R = teleport(*t, C, &n, dim, 1);
      int i; /* Make sure p gets removed */
      /* This asserts point in tree, will segfault */
      for(i = 0;
	  !pointEquals(&(R[i][1]), p, dim+1);
	  i++) ;
      if(pointEquals(&(R[i][1]), p, dim+1)){
	R[i][0]--; /* Decrease weight */
	if(0 == R[i][0]){ /* Remove on empty */
	  for(;i+1<n; i++)
	    R[i] = R[i+1];
	  n--;
	}
      }
      *t = NULL;
      if(0 < n){
	int *Acc = accSum(R, n);
	*t = buildBalanced(R, 0, n-1, dim, dim+1, Acc);
	free(Acc);
      }
      if(NULL != R)
	free(R);
    } else { /* Tree is balanced */
      /* Now deal nodeFound */
      (*t)->w--; /* Remove point */
      assert(0 < (*t)->w && "Deleting weight to 0");
      if(0 < dim){ /* Recursive into lower dim. */
	if(nodeFound)
	  deleteCR(&((*t)->equal), depth, p, C, dim-1);
	deleteCR(&((*t)->lower), depth, p, C, dim-1);
      }
      if(edgeFound && !nodeFound){
	(*sel)->w *= -1;
	if(deleteORA((ora)*sel, p))
	  *sel = NULL;
	else
	  (*sel)->w *= -1;
      }
    }
  }
}

void
delete(ort rt,
       point p
       )
{
  /* printf(" >D"); */
  /* for(int i = 0; i < rt->d ; i++) */
  /*   printf(" %d,", p[i]); */
  /* printf("\n"); */
#ifndef NDEBUG
  /* dotShow(rt->root); */
  checkTree(rt->root, rt->d-1);
#endif /* NDEBUG */

  if(4*(rt->n-1) <= rt->ca){
    rt->ca /= 2;
    if(NULL != rt->cache)
      free(rt->cache);

    rt->cache = NULL;
    if(0 < rt->ca)
      rt->cache = malloc(rt->ca*(1+rt->d)*sizeof(int));
  }
  rt->n--;

  if(0 >= adjustCut(1)) /* No cut-off struct */
    deleteR(&(rt->root), 1, p, rt->cache, rt->d-1);
  else /* With cut-off struct */
    deleteCR(&(rt->root), 1, p, rt->cache, rt->d-1);

#ifndef NDEBUG
  /* dotShow(rt->root); */
  if(0 == rt->n)
    assert(NULL == rt->root && "Missed point count on delete.");
  else
    assert(rt->n == abs(rt->root->w) && "Missed point insert.");
  assert(!containsQ(rt, p) && "Delete failed");
  checkTree(rt->root, rt->d-1);
#endif /* NDEBUG */
}
