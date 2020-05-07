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

#include "sortedList.h"
#include "string.h"

 /* A dynamic string representation. Letters get removed from the
    beginning and added at the end. */
struct string{
  int b; /* External begin value, these increase. Inclusive */
  int e; /* External end value. Increasing. Exclusive. */
  char* A; /* The string buffer. String starts at S[b-offset]. */
  int Sa;  /* Number of alloced bytes. */
  int sigma; /* Size of the alphabet */
  sortedList *Nxt; /* Sorted Lists, indexed by letters */
};

string
stringAlloc(int sigma
	    )
{
  string S;

  S = calloc(1, sizeof(struct string));

  S->sigma = sigma;
  S->Nxt = malloc(sigma*sizeof(sortedList));
  S->Nxt = &(S->Nxt[-'A']);
  int c = 'A';
  for(int i = 0; i < S->sigma; i++){
    S->Nxt[c] = listAlloc();
    c++;
  }

  return S;
}

void
stringFree(string S)
{
  free(S->A);
  int c = 'A';
  for(int i = 0; i < S->sigma; i++){
    listFree(S->Nxt[c]);
    c++;
  }

  free(&(S->Nxt['A']));
  free(S);
}

void
printString(string S
	    )
{
  for(int i=S->b; i<S->e; i++)
    printf("%c", stringLetter(S, i));
  printf("\n");
}

int
stringSize(string S
	   )
{
  return S->e-S->b;
}

int
stringSigma(string S
            )
{
  return S->sigma;
}

int
stringBegin(string S
	    )
{
  return S->b;
}

int
stringEnd(string S
	    )
{
  return S->e;
}

/* String resizing function. */
static void
stringRes(string S, /* The previous string */
	  int n     /* The new size */
	  )
{
  char* C = NULL;

  if(0 < n)
    C = malloc(n*sizeof(char)); /* new array. */

  int i;
  for(i = S->b; i < S->e; i++)
    C[i%n] = S->A[i%S->Sa];

  S->Sa = n;
  if(NULL != S->A)
    free(S->A);
  S->A = C;
}

/* Append letter c to the string */
void
stringAppend(string S,
	     char c
	     )
{
  if(S->e-S->b == S->Sa){ /* resize */
    if(0 == S->Sa)
      S->Sa = 1;
    stringRes(S, 2*S->Sa);
  }

  S->A[S->e % S->Sa] = c; /* Append letter */
  listAppend(S->Nxt[(int)c], S->e);
  S->e++;
}

/* Get a letter from the string */
char
stringLetter(string S,
	     int j
	     )
{
  return S->A[j % S->Sa];
}

char
stringFstLetter(string S)
{
  return S->A[S->b % S->Sa];
}

int
stringLast(string S,
	   char c
	   )
{
  return listLast(S->Nxt[(int)c]);
}

int
stringCeil(string S,
	   char c,
	   int k
	   )
{
  return listCeil(S->Nxt[(int)c], k);
}

int
stringFloor(string S,
	    char c,
	    int k
	    )
{
  return listFloor(S->Nxt[(int)c], k);
}

int
stringIdx(string S,
	   char c,
	   int k
	   )
{
  return listIdx(S->Nxt[(int)c], k);
}

int
stringPos(string S,
          char c,
          int i
          )
{
  return listAccess(S->Nxt[(int)c], i);
}

void
stringPop(string S)
{
  listPop(S->Nxt[(int)stringFstLetter(S)]);
  S->b++;

  if(4*(stringSize(S)) <= S->Sa) /* resize */
    stringRes(S, 1+(S->Sa/2));
}

char *
stringGet(string S
	  )
{
  int n = stringSize(S);
  char *R = malloc((n+1)*sizeof(char));

  int i = 0;
  for(int j = S->b; j < S->e; j++)
    R[i++] = stringLetter(S, j);
  R[i] = '\0';

  return R;
}
