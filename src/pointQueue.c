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

#include "pointQueue.h"

struct pointQueue{
  int in; /* Where points go in */
  int out; /* Where points go out */
  int Qa;  /* Alloced size */
  int Qm;  /* The marked position */
  point *Q; /* The queue */
};

pointQueue
allocPQ(void
        )
{
  pointQueue Q = calloc(1, sizeof(struct pointQueue));

  return Q;
}

void
freePQ(pointQueue Q
       )
{
  if(NULL != Q->Q)
    free(Q->Q);
  free(Q);
}

int
isEmptyPQ(pointQueue Q
          )
{
  return Q->out == Q->in;
}

void
markPQ(pointQueue Q)
{
  Q->Qm = Q->in;
}

int
markingTruePQ(pointQueue Q
              )
{
  return Q->Qm == Q->out;
}

/* Number of free positions */
static int
freePosPQ(pointQueue Q)
{
  return Q->out + Q->Qa - Q->in;
}

/* Transfer the points from the queue to a new array.
   of size n. */
static void
transfer(pointQueue Q,
         int n
         )
{
  point *nQ = malloc(n*sizeof(point));
  for(int i = Q->out; i < Q->in; i++)
    nQ[i % n] = Q->Q[i % Q->Qa];
  free(Q->Q);
  Q->Q = nQ;
  Q->Qa = n;
}

void
expandPQ(pointQueue Q,
         int n
         )
{
  if(n >= freePosPQ(Q)){ /* Increase size */
      int newSize = Q->in + n - Q->out;
      if(0 == newSize)
        newSize = 1;
      newSize *= 2;
      transfer(Q, newSize);
  }
}


/* Stores a pointer */
void
pushPQ(pointQueue Q,
      point p
     )
{
  if(0 == freePosPQ(Q))
    expandPQ(Q, 1);

  Q->Q[Q->in++ % Q->Qa] = p;
}

void
popPQ(pointQueue Q
     )
{
  Q->out++;

  if(Q->Qa > 4*(Q->in-Q->out)){
    int newSize = Q->Qa/2;
    if(0 == newSize)
      newSize = 1;
    transfer(Q, newSize);
  }

  while(Q->out >= Q->Qa){
    Q->out -= Q->Qa;
    Q->in -= Q->Qa;
    Q->Qm -= Q->Qa;
  }
}

point
topPQ(pointQueue Q
      )
{
  return Q->Q[Q->out];
}
