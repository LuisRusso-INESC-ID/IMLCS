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

#ifndef _ORT_H
#define _ORT_H

typedef struct ort *ort;
typedef struct node *node;

void
configure(int mxn, /* Maximum number of elements in the tree. */
          int mxd, /* Maximum depth of the tree */
          double alpha, /* Initial alpha value */
          int confcut
          );

/* Free the corresponding structure */
void
deconfigure(void);

double
adjustCut(int param
	  );

/* Create a new ORT */
ort
allocORT(int dim /* Number of dimensions */
	 );

/* Free the corresponding struct*/
void
freeORT(ort rt
	);

/* Returns the root weight */
int
weightORT(ort rt
	  );

/* Range Queries */
int
countQ(ort rt, /* The orthogonal range tree */
       int* coords /* Point coordinates, in LSD order. */
       );

/* Important for application */
int
containsQ(ort rt, /* The orthogonal range tree */
	  int* coords /* Point coordinates, in LSD order. */
	  );


#include "point.h"

/* Returns an array with the points that dominate
   the coordinates. */
point *
collect(ort rt, /* The orthogonal range tree */
	int* coords, /* Point coordinate */
	int* n /* Number of points */
	);

/* Returns the points that are dominated by the coords */
point *
dominatedCollect(ort rt, /* The orthogonal range tree */
                 int* coords, /* Point coordinate */
                 int* n /* Number of points */
                 );

point *
rangeCollect(ort rt, /* The orthogonal range tree */
	     int* minCoords, /* Point coordinate */
	     int* maxCoords, /* Point coordinate */
	     int* n /* Number of points */
	     );

void
insert(ort rt, /* The orthogonal range tree */
       point p
       );

void
delete(ort rt, /* The orthogonal range tree */
       point p
       );

#endif /* _ORT_H */
