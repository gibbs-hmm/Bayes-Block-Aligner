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

/**************************************************************************/
/* $Id: common.c,v 1.1 1997/02/14 20:28:36 junzhu Exp junzhu $            */
/*                                                                        */
/* Author:   Eric C. Rouchka, 1996.6                                      */
/*           Jun Zhu, 1996.8                                              */ 
/*                                                                        */
/* Description :  common subroutines                                      */
/**************************************************************************/


#include <string.h>
#include "common.h"

#define r_off   12
#define DIM(A) (sizeof(A)/sizeof((A)[0]))
/**----------------------------------------------------------------------**/
/**                        RANDOM NUMBER GENERATOR DECLARATIONS          **/
/**----------------------------------------------------------------------**/

static long     state[33] = {
  (long)0xd53f1852, (long)0xdfc78b83, (long)0x4f256096,  (long)0xe643df7,
  (long)0x82c359bf, (long)0xc7794dfa, (long)0xd5e9ffaa,  (long)0x2c8cb64a,
  (long)0x2f07b334, (long)0xad5a7eb5, (long)0x96dc0cde,  (long)0x6fc24589,
  (long)0xa5853646, (long)0xe71576e2, (long)0xdae30df,  (long)0xb09ce711,
  (long)0x5e56ef87, (long)0x4b4b0082, (long)0x6f4f340e,  (long)0xc5bb17e8,
  (long)0xd788d765, (long)0x67498087, (long)0x9d7aba26,  (long)0x261351d4,
  (long)0x411ee7ea, (long)0x393a263,  (long)0x2c5a5835,  (long)0xc115fcd8,
  (long)0x25e9132c, (long)0xd0c6e906, (long)0xc2bc5b2d,  (long)0x6c065c98,
  (long)0x6e37bd55 };

static long     *rJ = &state[r_off], *rK = &state[DIM(state)-1];



/**************************** SEQ2CH ************************************/
/* For convenient in calculation, we numerate amino acid to 1 to 20     */
/* corresponding to a to t.  Some special symbol we need to take care   */
/* are V,W,Y.                                                           */
/************************************************************************/

short SEQ2CH(char seq_type,short x)
{
  if (seq_type) { /* 1 for protein  0 for DNA */

  if(x==1) 
    {/* symbol for V  */
      return 21;
    }
  else if(x==9)
    { /* symbol for W */
      return 22;
    }
  else if(x==14)
    { /* symbol for Y */
      return 24;
    }
  else /* other symbols */
    return x;
  }
  else {
    if(x==1)
      { /* symbol for T */
       return 19;
      }
    else if(x==3)
      { /* symbol for G */
       return 6;
      }
    else if(x==4)
      { /* symbol for n */
        return 13;
      }
  }
 
}


void print_usage_bernoulli(void)
{
    printf("USAGE (site sampler) :: Gibbs -PBernoulli file lengths {flags}\n");
    printf("USAGE (motif sampler):: Gibbs -PBernoulli file lengths expect {flags}\n\n");
    printf("lengths = <int>[,<int>] : width of motif to be found\n");
    printf("expect  = <int>[,<int>] : expect number of motif elements\n");
    printf("\npossible flags:\n\n");
    printf("-W <pseudo_wt>       pseudosite weight (between 0 and 1)\n");
    printf("-w <pseduo_cnt_wt>   pseduocount weight\n");
    printf("-i <iteration_num>   number of iterations to try\n");
    printf("-S <num_Seeds>       number of seeds to try\n");
    printf("-s <seedval>         random number generator seed \n");
    printf("-p <plateau_per>     number of periods a maximum value");
    printf(" hasn't changed\n");
    printf("-P <prior_filename>      file of informative priors\n");
    printf("-o <out_filename>        file where results will be written\n");
    printf("-C <cutoff_value>        cutoff for near optimal sampler\n");
    printf("-c <mnum, beg, end>*     collapse alphabet between beg and end\n");
    printf("-R <mnum, beg, end>*     palandromic model between beg and end\n");
    printf("-a <mnum, beg, end, pal>*  Concentrate between beg and end\n");
    printf("-r                  turn off reverse complements with DNA\n");
    printf("-n                  Use nucleic acid alphabet\n");
    printf("-F                  Do not use fragmentation\n");
    printf("-x                  Do not remove low complexity regions\n");
    printf("-m                  Do not maximize after near optimal sampling\n");
    printf("-e                  Run Expectation/maximization algorithm\n");
    printf("-t                  Display sites used in near optimal sampling\n");
    printf("-I                  interactive mode\n");
    printf("-h                  help flag\n\n");
}

void print_usage_sankoff(void)
{
  printf("bayesaligner  -QUERY=queryfile -DATA=database\n");
  printf(" [-BASE_BLOCKSIZE=base_blocksize(def=15)]\n");
  printf(" [-MAX_NUM_BLOCKS=upper limit on number of blocks(def=20)]\n");
  printf(" [-NON_PROTEIN]\n");
  printf(" [-MATRIX_NUM=matrix_num]\n");
  printf(" [-DNA_MATRIX=matrix_start-matrix_end]\n");
  printf(" [-OUTFILE=outputfile]\n");
  printf(" [-BEST_ALIGN]\n");
  printf(" [-OUT_SEQ]\n");
  printf(" [-BACK_S_SIZE=back_sample_size]\n");
  printf(" [-BACK_S_OUTSEQ] \n");
  printf(" [-BACK_CUTOFF]\n");
  printf(" [-BACK_C_VALUE=back_cutoff_value]\n");
  printf(" [-MAX_SPEED]\n");
  printf(" [-EXACT_BLOCK]\n");
  printf(" [-NO_PROFILE]\n");
  printf(" [-DEBUG]\n");
}

void print_parameters_sankoff(FILE *fptr, Sank_S SK)
{
  int      i;
  char     **matrix_name=SK->rel_matrix.matrix_name;

  if(SK->flags.protein_sequence){
      fprintf(fptr, "#----parameters setting (protein sequence)---------\n");
  }
  else{
      fprintf(fptr, "#----parameters setting (DNA sequence)---------\n");
  }
  fprintf(fptr,"#-BASE_BLOCKSIZE=%d\n",SK->Sankoff_blocksize);
  fprintf(fptr,"#-SITE=");
  for(i=0;i<NInput_N;i++){
    fprintf(fptr,"%d ",Input_N[i]);
  }
  fprintf(fptr,"\n");  
  fprintf(fptr,"#-SITE_WEIGHT=%f\n",Input_alpha);  
  fprintf(fptr,"#-MATRIX=");
  for(i=0;i<SK->rel_matrix.Nmatrix;i++){
      fprintf(fptr,"%s, ",matrix_name[i]);
  }
  fprintf(fptr,"\n");
  if(SK->flags.back_sampling){
    fprintf(fptr,"#-BACK_S_BLOCKS=%d\n",SK->backsampling_blocks);
    fprintf(fptr,"#-BACK_S_SIZE=%d\n",SK->backsampling_number);
  }
  if (SK->flags.back_sampling_proportional) {
    fprintf(fptr,"#-BACK_S_PROPORTIONAL=TRUE\n");
  } 
  if(SK->flags.min_memory)
    fprintf(fptr,"#-MIN_MEMORY=TRUE\n");

  fprintf(fptr,"\n");
}



/*void print_EnterArgs(void)
{
   printf("\n\n\n\tARGUMENTS SHOULD BE ENTERED IN THE FOLLOWING USAGE :: \n\n");
   printf("(site sampler)  :: file lengths {flags}\n");
   printf("(motif sampler) :: file lengths expect {flags}\n");
   print_flags();
   printf("\n\nEnter in the arguments followed by <return> : \n");
}*/

void p_error(char *msg)
{
   printf("ERROR :: %s\nEXITING\n", msg);
   exit(0);
}

void p_usage_error(char *msg)
{
   printf("ERROR :: %s\n\n", msg);
   print_usage_bernoulli();
   printf("EXITING\n\n");
   exit(0);
}


/************************** findMaxMotifLen **********************/
/*                                                               */
/*   Find maximum motif length                                   */
/*****************************************************************/

int findMaxMotifLen(IPtype IP)
{
   int i;
   int max=0;
   for(i = 0; i < IP->nNumMotifTypes; i++) {
      if(IP->nMotifLen[i] > max)
         max = IP->nMotifLen[i];
   }
   return max;
}

/************************** findMaxNumMotif  *********************/
/*                                                               */
/*   Find maximum number of motif in a single type motif         */
/*****************************************************************/

int findMaxNumMotif(IPtype IP)
{
   int i;
   int max=0;
   int value;

   for(i = 0; i < IP->nNumMotifTypes; i++) {
      value = IP->nNumMotifs[i][FORWARD];
      if(IP->RevComplement)
         value += IP->nNumMotifs[i][REVERSE];
      if(value > max)
         max = value;
   }
   return max;
}

/******************* complement *****************************/
/*                                                          */
/*     find the complementary nucleotide                    */
/************************************************************/

char complement(char ch)
{
   char retval;

   switch(ch) {
      case 'a' :
         retval = 'b'; break;
      case 'b' :
         retval = 'a'; break;
      case 'c' :
         retval = 'd'; break;
      case 'd' :
         retval = 'c'; break;
      case 'n' :
         retval = 'e'; /* represent X here */
	 break;
      case 'x' :
         retval = 'e'; /* represent X here */
	 break;
      case 'X':
         retval = 'e'; /* represent X here */
	 break;	
      default:
         printf("CH = %c", ch);
         p_error("INVALID NUCLEOTIDE ENCOUNTERED"); break;
   }
   return retval;
}


/*
        Additive random number generator

        Modelled after "Algorithm A" in
        Knuth, D. E. (1981). The art of computer programming, volume 2, page 27
.

        7/26/90 WRG
*/

void sRandom(long x)
{
        register long   i;

        state[0] = x;
        /* linear congruential initializer */
        for (i=1; i<DIM(state); ++i) {
                state[i] = 1103515245*state[i-1] + 12345;
        }

        rJ = &state[r_off];
        rK = &state[DIM(state)-1];

        for (i=0; i<10*DIM(state); ++i) (void) Random();
}


/*
        Random --  return value in the range 0 <= x <= 2**31 - 1
*/
long Random(void)
{
        register long   r;

        r = *rK;
        r += *rJ--;
        *rK-- = r;

        if (rK < state)
                rK = &state[DIM(state)-1];
        else
                if (rJ < state)
                        rJ = &state[DIM(state)-1];
        return (r>>1)&0x7fffffff; /* discard the least-random bit */
}

void rm_file(char *fname)
{
   char *tmp;

   NEW(tmp, 4 + strlen(fname), char);
   sprintf(tmp, "rm %s", fname);
   system(tmp);
   free(tmp);
}


