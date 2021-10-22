/*


Bayes Aligner: A bayesian sequence alignment program
(Zhu, J.;Liu, J.S.; Lawrence, C.E. Bayesian adaptive sequence alignment algorithms. Bioinformatics, 14:25-39, 1998)
Please acknowledge the program authors on any publication of scientific 
results based in part on use of the program and cite the above article in 
which the program was described.

Copyright (C) 2006   Health Research Inc. 
                                                                    
HEALTH RESEARCH INCORPORATED (HRI),                   
ONE UNIVERSITY PLACE, RENSSELAER, NY 12144-3455.                  
                                                                   
Email:  Bayes-Block@wadsworth.org
                                       
This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.



*/

/* defines.h - generic codes and constants for afn biosequence programs. */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include <float.h>
#include <math.h>

/* VALUES */
#define NIL    		-1
#define FALSE		0	
#define TRUE		1	

#define ILLEGAL		-1.0
#define BELL		((char) 7)	

#define FILE_BEGIN	0
#define FILE_CURRENT	1
#define FILE_END	2

#define Boolean		char

/* CONSTANTS */
#define MAX_INTEGER		LONG_MAX

/* MACROS - standard macro definitions and static types */
#define	MEW(x,n,t)	(( (x=(t*) malloc(((n)*sizeof(t))))==NULL) ? \
			 (t*) (fprintf(stderr,"Out of Memory."),exit(1),0):x)

#define	NEW(x,n,t)	(( (x=(t*) calloc(n,sizeof(t)))==NULL) ? \
			 (t*) (fprintf(stderr,"Out of Memory."),exit(1),0):x)

#define	NEWP(x,n,t)	(( (x=(t**) calloc(n,sizeof(t*)))==NULL) ? \
			(t**) (fprintf(stderr,"Out of Memory."),exit(1),0):x)

#define	NEWPP(x,n,t)	(( (x=(t***) calloc(n,sizeof(t**)))==NULL) ? \
			(t***) (fprintf(stderr,"Out of Memory."),exit(1),0):x)

#define	NEWP3(x,n,t)	(( (x=(t****) calloc(n,sizeof(t***)))==NULL) ? \
			(t****) (fprintf(stderr,"Out of Memory."),exit(1),0):x)

#define	GETCHAR(m,C)	do{ fprintf(stderr,"%s ", m); \
			  if(fscanf(stdin,"%c",(C)) == 1) { \
                	    while(getchar()!='\n') if(feof(stdin)) exit(1);\
			    break;\
 			  } while(getchar()!='\n') if(feof(stdin)) exit(1);\
			} while(TRUE);

#define	GETINT(m,i)	do{ fprintf(stderr,"%s ",m); \
			  if(fscanf(stdin,"%d",(i)) == 1) { \
                	    while(getchar()!='\n') if(feof(stdin)) exit(1);\
			    break;\
 			  } while(getchar()!='\n') if(feof(stdin)) exit(1);\
			} while(TRUE);

#define print_error(str) for(fprintf(stderr,"%s\n",str); TRUE; exit(1))
#define DIGIT2INT(c)    ((int)(c - 48))
#define MIN(t,x,y)	(((t)(x) < (t)(y)) ? (t)(x) : (t)(y))
#define MAX(t,x,y)	(((t)(x) > (t)(y)) ? (t)(x) : (t)(y))
#define SUM(x)		(((x) * (x+1)) / 2)

