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

/********************************************************************/
/* $Id: Sankoff.h,v 1.17 1997/09/03 19:46:49 junzhu Exp junzhu $    */
/* Date: 1997/07/03                                                 */
/* Author: Jun Zhu                                                  */
/* Description: header file for Sankoff (Bayes aligner) program     */
/*  modified 6/12/98 to set change block limit from 20 to user defined*/
/********************************************************************/

#ifndef SANKOFF_H
#define SANKOFF_H

#include "common.h"


#define BASE_BLOCKSIZE 15
#define MAX_NUM_MATRIX_B 8
#define MAX_NUM_MATRIX_P 27
#define DEFAULT_NUM_BLOCKS 20
#define MAX_FILE_NAME_LENGTH 200
#define MAX_STRING_LENGTH 1000
#define MAX_NUMBER_LENGTH 200
#define MAX_NUMBER_ARGS 50
#define SPARSE_CUTOFF_LIMIT .005


/* structure which contains all the option flags */
typedef struct {
   char best_alignment;             /* calculate the best alignment? */
   char output_sequence;            /* output best aligned sequence? */
   char back_sampling;              /* back sampling?                */
   char back_sampling_sequence;     /* output sampled sequences?     */
   char back_sampling_cutoff;       /* sample all the alignment above*/
                                    /* a cutoff value?               */
   char debug;                      /* output detail information at  */
                                    /* each step?                    */
   char exact_posterior;            /* calculate the exact posterior?*/
   char protein_sequence;           /* protein or DNA                */
   char min_memory;                 /* hack to reduce memory requirements */
   char back_sampling_proportional;   /* backsample propotional to probability*/
   char profile;                    /* print out profile ?           */
   char full_profile;               /* full or sparse profile?       */
} Sankoff_Option;

typedef struct link_list *pt_save_loop;

typedef struct link_list{
   short startb_x;
   short startb_y;
   int last_step;
   int Nblocks;
   int k;
   short **start_x;                 /*------------------------------------*/
   short **start_y;                 /* terminals of matching blocks       */
   short **end_x;                   /*                                    */
   short **end_y;                   /*------------------------------------*/
   double PP;                       /* posterior probability of alignment */
   int done;                        /* has it been processed ?            */
   int id;                          /* debugging id                       */
   pt_save_loop next;
} save_loop;   

/* structure which contains information related to relation matrice  */
/* such as BLOSUM and PAM matrices.                                  */
enum Matrix_type {BLOSUM, PAM};
typedef struct {
    int    matrix_num;               /* indicate which matrices to be used */
    int    index_start;              /* index starting number              */
    int    index_step;               /* index step size                    */
    int    Nmatrix;                  /* how many matrices used             */
    double ***pp_related;            /* relation matrices for protein or   */
    float  ***ff_related;            /* DNA. pp_related is in probability  */
                                     /* ff_related in bit unit.            */
    char   **matrix_name;            /* name of the relation matrix        */ 
} Sankoff_relation_matrix;
 
/* main structure for Bayes aligner */
typedef struct {
   int Sankoff_max_blocks;          /* upper limit on Sankoff_blocks (csc)*/
   int Sankoff_blocks;              /* maximum number of blocks           */
   int Sankoff_blocksize;           /* base block size to calculate number*/
                                    /* of blocks.                         */
   FILE* queryfile_ptr;             /* file pointer for query file        */
   char comment[2][50];             /* space for store the comment(header)*/
   int backsampling_blocks;         /* maximum blocks for back sampling   */
                                    /* procedures                         */
   int backsampling_number;         /* number of back sampling            */
   float backsampling_cutoff_value; /* cutoff for tracing all possible    */
                                    /* alignments                         */ 
   double **profile;                /* sampling result or marginal joint  */
                                    /* posterior distribution             */

   float ****SV;                    /* best matching score with t matching*/
   float ****SW;                    /* blocks (See Sankoff, 1972)         */

   float *S_END;                    /* best matching score at the very end*/

   short **start_x;                 /*------------------------------------*/
   short **start_y;                 /* terminals of matching blocks       */
   short **end_x;                   /*                                    */
   short **end_y;                   /*------------------------------------*/

   pt_save_loop top;                  /* top of save struct for min_bayes   */
   
   double PM;                       /* margial probability P(R1,R2)       */
   double *PK;                      /* posterior probability of k matching*/
                                    /* blocks                             */
   double *PSI;                     /* posterior probability of matrix psi*/
   double **PKSI;                   /* posterior probability of k matching*/
                                    /* blocks using matrix psi            */

   double ***WD;                    /*------------------------------------*/
   double ***WR;                    /* Sum of all probabilitoes with t    */
   double ***WC;                    /* matching blocks end with going     */
                                    /* down, right or central             */
                                    /*------------------------------------*/

   double **W;                      /* WC+WD+WR                           */

   double ***ND;                    /*------------------------------------*/
   double ***NR;                    /* Number of all prossible alignments */
   double ***NC;                    /*  with t matching blocks and last   */
                                    /* step going down, right or central. */
                                    /*------------------------------------*/
      
   double *N;                       /* NC+ND+NR                           */

   double ***BWD;                   /*------------------------------------*/
   double ***BWR;                   /* backward sum of all probabilities  */
   double ***BWC;                   /* with t matching blocks and last    */
                                    /* step going down, right or central  */
                                    /* corresponding to forward direction */ 
                                    /*------------------------------------*/

   int    **factor;                 /* factor of decrease for W           */
                                    /* to prevent overflow, W is mutiple  */
                                    /* 1.0E-300 whenever it above some    */
                                    /* threathod (1.0e100)                */

   float  *prior_N;                 /* prior probability of residue i in  */
                                    /* query sequence included in an      */
                                    /* alignment                          */

   Sankoff_relation_matrix rel_matrix; /* relation matrix                 */

   Sankoff_Option flags;            /* option flags                       */

} Sankoff_struct;

typedef Sankoff_struct *Sank_S;


/* global variables */
int   NInput_N;                     /* number of input sites              */
int   Input_N[20];                  /* array contain input sites          */
float Input_alpha;                  /* extra weight for alignments with   */
                                    /* input sites                        */
char command_string[MAX_STRING_LENGTH]; /* pusedo command line */
char *sargv[MAX_NUMBER_ARGS];           /* pusedo arguments */
int  sargc;                            /* number of pusedo arguments    */
#endif 




