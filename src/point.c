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

#include "point.h"

static unsigned int sizeNode;

void
setSize(unsigned int n
        )
{
  sizeNode = n;
}

int
pointcmp(point p,
         point q,
         int dim
         )
{
  return p[dim] - q[dim];
}

int
pointFullcmp(point p,
	     point q,
	     int dim
	     )
{
  int r = 0;
  while(0 == r && 0 <= dim){
    r = pointcmp(p, q, dim);
    dim--;
  }

  return r;
}

int
pointEquals(point p,
	    point q,
	    int dim
	    )
{
  int r = 1;

  for(int i = 0; i < dim && r; i++)
    r = p[i] == q[i];

  return r;
}

int
getCoord(point p, /* Point */
         int d    /* dimension */
         )
{
  return p[d];
}

void
plusPlus(point p, /* Point */
         int dim  /* dimension */
         )
{
  for(int i = 0; i < dim; i++)
    p[i]++;
}

void
minusMinus(point p, /* Point */
           int dim  /* dimension */
           )
{
  for(int i = 0; i < dim; i++)
    p[i]--;
}

int
pointSum(point p,
	 int dim
	 )
{
  int r = 0;
  for(int i = 0; i < dim; i++)
    r += p[i];

  return r;
}
