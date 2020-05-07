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

#ifndef _MLCS_H
#define _MLCS_H

typedef struct mlcs *mlcs;

#ifndef NDEBUG
struct mlcs{
  int dim; /* The number of strings considered. */
  int zeros;  /* Number of strings of size 0 */
  string* S; /* Array that stores the current strings.  */
  int lambda; /* The size of the mlcs */
  ort *PF;    /* Array with pareto fronts. */
  int pfA;    /* Size of the PF array */
};
#endif /* NDEBUG */

void
printMLCS(mlcs m);

mlcs
allocMLCS(int dim, /* The number of strings. */
	  int sigma /* Alphabet size */
	  );

void
freeMLCS(mlcs m /* The structure */
	 );

string
pullString(mlcs m,
	   int t
	   );

void
append(mlcs m, /* The MLCS data struct */
       int i,  /* Which string */
       char c  /* Which letter */
       );

void
pop(mlcs m, /* The MLCS data struct */
    int i  /* Which string */
    );

/* Returns the size of the MLCS */
int
mlcsSize(mlcs m
	 );

#endif /* _MLCS_H */
