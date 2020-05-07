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
#include <stdio.h>
#include <assert.h>

#include "sortedList.h"

struct sortedList{
  int b;  /* begin, exclusive index of where the list starts */
  int e;  /* end, inclusive index of where the list ends */
  int* A; /* The list with numbers */
  int al; /* Number of alloced positions */
};

sortedList
listAlloc(void
	  )
{
  sortedList L;

  L = calloc(1, sizeof(struct sortedList));

  return L;
}

void
listFree(sortedList L
	 )
{
  free(L->A);
  free(L);
}

static void
listRes(sortedList L, /* The previous string */
	int n     /* The new size */
	)
{
  int* C = NULL;

  if(0 < n)
    C = malloc(n*sizeof(int)); /* new array. */

  int i;
  for(i = L->b+1; i <= L->e; i++)
    C[i%n] = L->A[i%L->al];

  L->al = n;
  if(NULL != L->A)
    free(L->A);
  L->A = C;
}

/* Appends the value v to the list. */
void
listAppend(sortedList L,
	   int v
	   )
{
  if(L->e-L->b == L->al){ /* resize */
    if(0 == L->al)
      L->al = 1;
    listRes(L, 2*L->al);
  }

#ifndef NDEBUG
  if(L->b < L->e){
    assert(L->A[L->e % L->al] < v && "Inserting a small value in the ordered list.");
  }
#endif /* NDEBUG */

  L->e++;
  L->A[L->e % L->al] = v;
}

/* Remove the first element from the list. */
void
listPop(sortedList L
	)
{
  L->b++;

  if(4*(L->e-L->b) <= L->al){ /* resize */
    listRes(L, 1+(L->al/2));
  }
}

/* Returns the last element from the list */
int
listLast(sortedList L
	 )
{
  int r = -2; /* Default for empty list */

  if(L->b < L->e)
    r = L->A[L->e % L->al];

  return r;
}

int
listIdx(sortedList L,
        int k
        )
{ /* Do a binary search */
  int b = L->b; /* Begin exclusive */
  int e = L->e; /* End inclusive */

  while(b+1 < e){
    int m = (e+b)/2; /* Middle element */

    if(L->A[m % L->al] < k)
      b = m;
    else
      e = m;
  }

  return e;
}

/* Returns the first element of the list that is larger than or equal to
   k.

   Returns -2, to represent +infty. */
int
listCeil(sortedList L,
	 int k
	 )
{ /* Do a binary search */
  int r = -2;

  if(L->b < L->e && L->A[L->e % L->al] >= k)
    r = L->A[listIdx(L, k) % L->al];

  return r;
}

int
listFloor(sortedList L,
	 int k
	 )
{ /* Do a binary search */
  int r = -2;

  if(L->b < L->e && L->A[(L->b+1) % L->al] <= k){
    int idx = listIdx(L, k);
    r = L->A[idx % L->al];
    if(r > k)
      r = L->A[(idx-1) % L->al];
  }

  return r;
}

int
listAccess(sortedList L,
           int i
           )
{
  int r = -2;

  if(L->b < i && i <= L->e)
    r = L->A[i % L->al];

  return r;
}

#if 0 /* Comment */

int
main(void)
{
  sortedList L = listAlloc();

  for(int i = 0; i < 100 ; i++)
    listAppend(L, 2*i);

  for(int i = 0; i < 30 ; i++)
    listPop(L);

  for(int i = 0; i < 220 ; i++){
    printf("Last: %d\n", listLast(L));
    printf("Ceil of %d is %d\n", i, listCeil(L, i));
  }

  listFree(L);

  return 0;
}

#endif /* Comment */
