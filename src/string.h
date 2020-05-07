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

#ifndef _STRING_H_MLCS
#define _STRING_H_MLCS

typedef struct string *string;

/* Allocs a string data structure. Assuming
   that the underlying alphabet has size sigma
   and starts at 'A'.
 */
string
stringAlloc(int sigma
	    );

void
stringFree(string S
	   );

void
printString(string S
	    );

int
stringSize(string S
	   );

int
stringSigma(string S
            );

int
stringBegin(string S
	    );

int
stringEnd(string S
	    );

void
stringAppend(string S,
	     char c
	     );

/* Get a letter from the string */
char
stringLetter(string S,
	     int j
	     );

char
stringFstLetter(string S
		);

/* What is the position of the last occ of letter c */
int
stringLast(string S,
	   char c
	   );

/* What is the next occ of letter c */
int
stringCeil(string S,
	   char c,
	   int k
	   );

/* What is the prev occ of letter c */
int
stringFloor(string S,
	    char c,
	    int k
	    );

/* Index of the next occ of letter c */
int
stringIdx(string S,
	   char c,
	   int k
	   );

/* Returns a string position where letter c occurs.
   The position is the i-th in its list */
int
stringPos(string S,
          char c,
          int i
          );

void
stringPop(string S
	  );

char *
stringGet(string S
	  );

#endif /* _STRING_H */
