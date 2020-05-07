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

#ifndef _ORA_H
#define _ORA_H

typedef struct ora *ora;

void
dsRnORA(ora R,
	FILE *f
	);

/* Create a new ORA */
ora
allocORA(int dim /* Number of dimensions */
	 );

/* Free the corresponding struct*/
void
freeORA(ora R
	);

void
checkORA(ora R
	 );

int
weightORA(ora R
	  );

/* Range Queries */
int
countQORA(ora R, /* The orahogonal range tree */
          int* coords /* Point coordinates, in LSD order. */
          );

/* Deprecated Not in use */
int
containsQORA(ora R, /* The orthogonal range tree */
             int* coords /* Point coordinates, in LSD order. */
             );

#include "point.h"

ora
buildORA(point *C, /* Array with point pointers */
	 int l, /* Left index, inclusive */
	 int r, /* Right index, inclusive */
	 int dim /* Current dimension */
	 );

void
teleportORA(ora R, /* Tree node to traverse */
            int *C, /* Array for storing points and multipls */
            int *i, /* at the end it is the size */
            int dim, /* Current dimension */
            int mxdim, /* Maximum dimension */
            point p,   /* Current common point coords from
			  mxdim to dim */
	    int cp
            );

/* Copy to array the points that dominate
   the coordinates. */
void
collectORA(ora R, /* The orthogonal range array */
	   int *C,   /* Array for storing points. */
	   int *coords, /* Point coordinates, in LSD order. */
	   int maxdim,
	   point hp, /* Temporary coords */
	   int *n
           );

void
rangeCollectORA(ora R, /* The orthogonal range array */
		int *C,   /* Array for storing points. */
		int *minCoords, /* Point coordinates, in LSD order. */
		int *maxCoords, /* Point coordinates, in LSD order. */
		int maxdim,
		point hp, /* Temporary coords */
		int *n
		);

/* Copy to array the points dominated by
   the coordinates. */
void
dominatedCollectORA(ora R, /* The orthogonal range array */
	   int *C,   /* Array for storing points. */
	   int *coords, /* Point coordinates, in LSD order. */
	   int maxdim,
	   point hp, /* Temporary coords */
	   int *n
           );

void
insertORA(ora R, /* The orthogonal range tree */
          point p,
	  int w
          );

int
deleteORA(ora R, /* The orthogonal range tree */
          point p
          );

#endif /* _ORA_H */
