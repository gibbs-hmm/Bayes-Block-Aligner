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

/*******************************************************************/
/* *
/* $Id: Sankoff.c,v 1.42 1997/09/03 20:01:40 junzhu Exp junzhu $   */
/* */
/* Author: Jun Zhu 1996/7/15                                        */
/*  Modified: Steve Carmack 1998/4/2                               */
/* */
/* */
/* Sankoff method: for detail read David Sankoff's paper in       */
/* Proc.Nat.Acad.Sci.USA Vol. 69, No.1, pp. 4-6, Jan. 1972        */
/* This is not strict Sankoff method, it is Sankoff-like method.  */
/* It uses the Sankoff_blocks as the limit of matching blocks.    */
/* */
/*******************************************************************/

#include <limits.h>
#include <math.h>
#include "common.h"
#include "blosum.h"
#include "Sankoff.h"
#include "mem_mgmt.h"
#include "routines.h"

void 
Sankoff(Model B, Sank_S SK)
{
   int             i, j, k, n;	/* general register */
   int             Next_Seq;	/* indicate if we have next sequence */
   int             Nmatrix;
   int           **s_pptr, *s_ptr;	/* pointer for quick access matrix */
   short          *x_ptr;
   double          sump;	/* sum of posterior probability */
   double          tmp, max_k;
   char          **matrix_name;
   /*   SimTime         time; */
   double          posterior;
   enum 
   {
      QUERY, DATABASE
   };
   int tused;         /*time elapsed */
   time_t tstart,tstop;

   time(&tstart); 
   /*  BeginTime(&time); *



   /* prepare the relation matrix */
   alloc_Sank_rel_matrix(B, SK);
   fill_Sank_rel_matrix(B, SK);
   Nmatrix = SK->rel_matrix.Nmatrix;

   /* print out general parameters */
   print_parameters_sankoff(B->IP->Datafiles->out_fpt, SK);

   /* read in query sequence */
   Sankoff_get_string(B, SK, QUERY);
   fprintf(B->IP->Datafiles->out_fpt, "Query file:\n");
   fprintf(B->IP->Datafiles->out_fpt,
	   "%s   LEN=%d\n", SK->comment[0], B->Seq->SeqLen[0]);
   fprintf(B->IP->Datafiles->out_fpt, "----------------------------\n\n");

   /* read in sequences in database one by one */
   Next_Seq = TRUE;
   while (Next_Seq)
   {				/* go scan it */
      Next_Seq = Sankoff_get_string(B, SK, DATABASE);
      fprintf(B->IP->Datafiles->out_fpt,
	      "\n%s      LEN=%d\n", SK->comment[1], B->Seq->SeqLen[1]);

      /* allocate space for Sankoff structure                   */
      /* as the space complexity is high, we can allocate some  */
      /* space at a time, then release if they are not needed.  */

      /* maximum number of blocks, the upper limit is 20 */
      if (B->Seq->SeqLen[0] > B->Seq->SeqLen[1])
      {
	 SK->Sankoff_blocks = B->Seq->SeqLen[1] / SK->Sankoff_blocksize;
      } else
      {
	 SK->Sankoff_blocks = B->Seq->SeqLen[0] / SK->Sankoff_blocksize;
      }
      /* only calculate what you need */
      if (SK->flags.back_sampling)
          SK->Sankoff_blocks = SK->backsampling_blocks;

      if (SK->Sankoff_blocks == 0)
	 SK->Sankoff_blocks = 1;
      if (SK->Sankoff_blocks > SK->Sankoff_max_blocks);
	 SK->Sankoff_blocks = SK->Sankoff_max_blocks;
      
      alloc_Sank(B, SK);	/* allocate basic memory  */
      alloc_Sank_N(B, SK);	/* allocate memory N      */
      Sankoff_probability_N(B, SK);	/* counting alignment     */
      free_Sank_N(B, SK);	/* free memory N          */
      alloc_Sank_W(B, SK);	/* allocate memory W      */
      Sankoff_probability_W(B, SK);	/* sum probability        */

      matrix_name = SK->rel_matrix.matrix_name;

      /* if no backward step free memory W */
      /*
       * if(!SK->flags.back_sampling && !SK->flags.exact_posterior )
       * free_Sank_W(B,SK);
       */
      /* optimal alignemnt */
      if (SK->flags.best_alignment)
      {				/* calculate the best alignment */
	 alloc_Sank_S(B, SK);	/* allocate memory for S        */
	 Sankoff_Make_S(B->Seq, SK);	/* calculate the matrix S       */
	 Sankoff_Trace_S(B, SK);/* find out the optimal alignment */

	 if (!SK->flags.back_sampling || !SK->flags.back_sampling_cutoff)
	    free_Sank_S(B, SK);
      } else
      {				/* just print out maginal likelihood */
	 fprintf(B->IP->Datafiles->out_fpt, "\n");

	 if (SK->flags.debug)
	 {			/* print detail information for debug */
	    for (n = 0; n < Nmatrix; n++)
	    {
	       fprintf(B->IP->Datafiles->out_fpt,
		       "%s:\n", matrix_name[n]);

	       for (k = 1; k <= SK->Sankoff_blocks; k++)
	       {
		  fprintf(B->IP->Datafiles->out_fpt,
			  "k=%d, W[k]=%g+%d, N[k]=%g, W[k]/N[k]=%g+%d\n",
		     k, SK->W[n][k], (SK->factor[n][k] + 1) * 300, SK->N[k],
		      SK->W[n][k] / SK->N[k], (SK->factor[n][k] + 1) * 300);
		  /* posterior of k>0 */
		  posterior = 0.0;
		  for (i = 1; i <= k; i++)
		     posterior += SK->PKSI[n][i];
		  posterior = posterior / (double) k;
		  fprintf(B->IP->Datafiles->out_fpt,
			  "P(k>0|R)=%g\n\n",
			  posterior / (SK->PKSI[n][0] + posterior));
	       }
	    }
	 } else
	 {			/* debug off */
	    if (SK->rel_matrix.Nmatrix == 1)
	    {
	       posterior = 0.0;
	       for (i = 1; i <= SK->Sankoff_blocks; i++)
		  posterior += SK->PKSI[0][i];
	       posterior = posterior / (double) SK->Sankoff_blocks;
	       fprintf(B->IP->Datafiles->out_fpt,
		       "P(k>0|R)=%g\n\n",
		       posterior / (SK->PKSI[0][0] + posterior));
	    }
	 }

	 /* print posterior P(k|R)  */
	 max_k = 0;
	 for (i = 1; i <= SK->Sankoff_blocks; i++)
	 {
	    if (max_k < SK->PK[i])
	       max_k = SK->PK[i];
	 }
         /*Fix the divide by zero problem*/
         if ( SK->PK[0] < 1.0e-307)
             SK->PK[0] = 1.0e-306;
	 tmp = max_k * (0.5 / SK->PK[0]) * (double) SK->Sankoff_blocks;
	 fprintf(B->IP->Datafiles->out_fpt,
		 "\nNumber of Blocks, k\n");
	 fprintf(B->IP->Datafiles->out_fpt,
		 "Bayesian Evidence against the null: P(k=k*|R)=%g\n",
		 tmp / (tmp + 0.5));
	 sump = 0;
	 for (i = 0; i <= SK->Sankoff_blocks; i++)
	 {
	    sump += SK->PK[i];
	    fprintf(B->IP->Datafiles->out_fpt,
		    "P(%d|R)=%g, P(k>%d|R)=%g\n",
		    i, SK->PK[i], i, 1. - sump);
	 }

	 if (SK->rel_matrix.Nmatrix > 1)
	 {			/* print posterior  P(psi|R) */
	    fprintf(B->IP->Datafiles->out_fpt,
		 "\n\nPosterior probability of score matrices, P(psi|R)\n");
	    for (n = 0; n < SK->rel_matrix.Nmatrix; n++)
	    {
	       fprintf(B->IP->Datafiles->out_fpt,
		       "%s: P(%d|R)=%g\n",
		       matrix_name[n], n, SK->PSI[n]);
	    }
	 }
      }

      if (SK->flags.exact_posterior)
      {
	 /* exact posterior */
	 alloc_Sank_BW(B, SK);
	 Sankoff_probability_BW(B, SK);
	 free_Sank_BW(B, SK);
	 if (!SK->flags.back_sampling)
	    free_Sank_W(B, SK);
      }
      if (SK->flags.back_sampling)
      {
	if (SK->backsampling_number > 0)
          { 
            if(SK->flags.min_memory)
   	        Sankoff_Back_Sampling_Min(B, SK);
            else
                Sankoff_Back_Sampling(B, SK);
	  }
	 if (SK->flags.best_alignment && SK->flags.back_sampling_cutoff)
	 {
	    Sankoff_Back_Sampling_cutoff(B, SK);
	    free_Sank_S(B, SK);
	 }
	 free_Sank_W(B, SK);
      }
      free_Sank(B, SK);
   }
  

   free_Sank_rel_matrix(B, SK);
   /* free memory */
   free(SK);
   if (B->IP->Datafiles->out_fpt != stdout)
      fclose(B->IP->Datafiles->out_fpt);
   free(B->IP->nPossSites);
   free(B->Seq->SeqLen);
   free(B->Seq->Orig);
   free(B->Seq->R[0]);
   free(B->Seq->nvEndLocs[0]);
   free(B->IP);
   free(B);
   /* EndTime(&time); */
   time(&tstop);
   tused = tstop - tstart;
   /*  tused = difftime(tstop,tstart); will not work on the sgi */

    printf("Job completed. Total Time %d sec (%f min)\n", tused, (double)tused/60.0);
   
}


/*******************  Sankoff_Back_Sampling ************************/
/*                                                                 */
/*  Back Sampling the alignment                                    */
/*  The results are saved in SK->start_x[0][k], and                */
/*  SK->end_x[0][k]                                                */
/*******************************************************************/

void Sankoff_Back_Sampling (Model B, Sank_S SK)
{
  int i,j,k,n,N;
  int start_x,start_y;
  int Nmatrix;
  int Len0,Len1;               /* sequence lengths */
  double* profile_ptr;         /* for quick access */
  float score;
  int Nblocks;
  int strl;
  int *Nback_sampling;         /* number of sample drawed for each matrix */
  int sample_limit;            /* maximum number of blocks to backsample */

  double tmp,tmp1,w;
  double ***WC,***WD,***WR;    /* for quick access */
  double **wd_pptr,**wd_pptr0; /* for quick access */
  double **wc_pptr,**wc_pptr0;    
  double **wr_pptr,**wr_pptr0;
  double ***pp_related;
  double sample_factor;        /* scale factor for random number(proportional)*/
  char   **matrix_name;
  short *Orig0=B->Seq->Orig[0];
  short *Orig1=B->Seq->Orig[1];
  float *prior_N=SK->prior_N;
  enum {DOWN, CENTRAL,RIGHT} last_step;
  Stype Seq=B->Seq;

  Len0=B->Seq->SeqLen[0];
  Len1=B->Seq->SeqLen[1];
  pp_related=SK->rel_matrix.pp_related;
  matrix_name=SK->rel_matrix.matrix_name;
  Nmatrix=SK->rel_matrix.Nmatrix;
  /* NEW(Nback_sampling,Nmatrix,int);*/
  Nback_sampling= (int *) calloc(Nmatrix,sizeof(int));
  /* allocate space for profile */
  NEWP(SK->profile,Len0+1 ,double);
  for(i=0;i<=Len0;i++){
      NEW(SK->profile[i],Len1+1,double);
  }

  /* initialize the profile */
  for(i=0; i<=Len0; i++) {
      profile_ptr=SK->profile[i];
      for(j=0; j<=Len1; j++) {
	  profile_ptr[j]=0;
      }
  }

  /* set seed */
  sRandom((long)time(NULL));

  fprintf(B->IP->Datafiles->out_fpt,"Back sampling.\n\n");
  if(Nmatrix>1) {/* multiple matrices, sampling according to P(psi|R)*/
      for(N=0; N<Nmatrix; N++) {
	  Nback_sampling[N]=(int)(SK->backsampling_number*SK->PSI[N]);
      }
  }
  else {
      Nback_sampling[0]=SK->backsampling_number;
  }

  for(N=Nmatrix-1; N>=0; N--) {
      if(N!=Nmatrix-1 &&  Nback_sampling[N]!=0){
	  /* need to recalculate the WC, WR,WD */
	  for(k=1;k<=SK->Sankoff_blocks;k++){
	      wd_pptr=SK->WD[k];
	      wd_pptr0=SK->WD[k-1];
	      wr_pptr=SK->WR[k];
	      wr_pptr0=SK->WR[k-1];
	      wc_pptr=SK->WC[k];
	      wc_pptr0=SK->WC[k-1];

	      /* keep track of decreasing factor */
	      SK->factor[N][k]= SK->factor[N][k-1]; 

	      tmp=wd_pptr0[Len0][Len1]+wr_pptr0[Len0][Len1]+
		  wc_pptr0[Len0][Len1];
	      if(tmp>1.0e30) { /*decrease by a factor */ 
		  for(i=1;i<Len0;i++) {
		      for(j=1;j<=Len1;j++){
			  wd_pptr0[i][j]=wd_pptr0[i][j]*1.0e-300;
			  wr_pptr0[i][j]=wr_pptr0[i][j]*1.0e-300;
			  wc_pptr0[i][j]=wc_pptr0[i][j]*1.0e-300;
		      }
		  }
		  SK->factor[N][k]++;
	      }
     
	      for(i=1;i<=Len0;i++) {
		  for(j=1;j<=Len1;j++){
		      wd_pptr[i][j]=wc_pptr[i-1][j]+wd_pptr[i-1][j];
		      wr_pptr[i][j]=wc_pptr[i][j-1]+wd_pptr[i][j-1]+
			            wr_pptr[i][j-1];
		      wc_pptr[i][j]=(wd_pptr0[i-1][j-1]+wr_pptr0[i-1][j-1]+
				     wc_pptr[i-1][j-1])*prior_N[i]*
			pp_related[N][Orig0[i-1]][Orig1[j-1]];
		  }
	      }
	  }
      }

      WC=SK->WC;
      WR=SK->WR;
      WD=SK->WD;
       
      fprintf(B->IP->Datafiles->out_fpt,
	      "%d samples are from Matrix %s.\n",
	      Nback_sampling[N],matrix_name[N]);
      if(SK->backsampling_blocks > SK->Sankoff_blocks) 
         (SK->flags.back_sampling_proportional = TRUE); /* preserve old option*/

      for(i=0;i<Nback_sampling[N];i++){
	  if(SK->flags.back_sampling_proportional){
	      /* sampling with random number of blocks according MagL */
	      /* get a random number between 0-1 */
	      tmp=Random()/(double)INT32_MAX;
	      /*            printf("%f\n",tmp); */
            if ( SK->backsampling_blocks > SK->Sankoff_blocks)
                sample_limit = SK->Sankoff_blocks;
            else
                sample_limit = SK->backsampling_blocks;
            
	    /* calculate factor  to scale random number */
         if (sample_limit != SK->Sankoff_blocks)
           {
            sample_factor = SK->PKSI[N][1];
            k = 1;
	    while (k < sample_limit )
            {
             k++;
             sample_factor += SK->PKSI[N][k];
	    }

            tmp = tmp * sample_factor;   /* scale the random number */
           }
	      tmp1=SK->PKSI[N][1];
	      k=1;
	      while(k < sample_limit && tmp1<tmp){
		  k++;
		  tmp1+=SK->PKSI[N][k];
	      }
	      Nblocks=k;
	  }
	  else{/* sampling with fix number of blocks */
	      Nblocks=SK->backsampling_blocks;
	  }

	  /* initializing */
	  for(k=1;k<=Nblocks;k++)
	      SK->start_x[0][k]=0;

	  start_x=Len0;
	  start_y=Len1;
	  k=Nblocks;
	  score=0;
	  last_step=RIGHT;

	  while(start_x >0 && start_y >0 && k>0){
	      /* get a random number between 0-1 */
	      tmp=Random()/(double)INT32_MAX;
	   
	      switch(last_step){
	      case DOWN:
		  w=WC[k][start_x][start_y]+
		    WD[k][start_x][start_y];	      
		  if(WC[k][start_x][start_y]/w >tmp){ 
		      /* take the center path */
		      /*score=score+ff_related[N]
		                       [B->Seq->Orig[0][start_x-1]]
		                       [B->Seq->Orig[1][start_y-1]];*/
		      SK->end_x[0][k]=start_x;
		      SK->end_y[0][k]=start_y;
		      SK->start_x[0][k]=start_x;
		      SK->start_y[0][k]=start_y;
		      start_x--;
		      start_y--;
		      last_step=CENTRAL;
		  }
		  else {/* go up */
		      start_x--;
		      last_step=DOWN;		  
		  }
		  break;
	      
	      case CENTRAL:
		  w=WC[k][start_x][start_y]+
		    WR[k-1][start_x][start_y]+
		    WD[k-1][start_x][start_y];
		  if(WC[k][start_x][start_y]/w >tmp){ 
		      /* take the center path */
		      /*score=score+ff_related[N]
		                       [B->Seq->Orig[0][start_x-1]]
		                       [B->Seq->Orig[1][start_y-1]];*/
		      SK->start_x[0][k]=start_x;
		      SK->start_y[0][k]=start_y;
		      start_x--;
		      start_y--;
		      last_step=CENTRAL;
		  }
		  else if((WC[k][start_x][start_y]+
			   WR[k-1][start_x][start_y])/w >tmp){
		      /* take a path to go left */
		      start_y--;
		      k--;
		      last_step=RIGHT;
		  }
		  else{ /* go up */
		      start_x--;
		      k--;
		      last_step=DOWN;		  
		  }
		  break;

	      case RIGHT:
		  w=WC[k][start_x][start_y]+
		    WR[k][start_x][start_y]+
		    WD[k][start_x][start_y];	      
		  if(WC[k][start_x][start_y]/w >tmp){
		      /* take the center path */
		      /*score=score+ff_related[N]
		                       [B->Seq->Orig[0][start_x-1]]
		                       [B->Seq->Orig[1][start_y-1]];*/	      
		      SK->end_x[0][k]=start_x;
		      SK->end_y[0][k]=start_y;
		      SK->start_x[0][k]=start_x;
		      SK->start_y[0][k]=start_y;
		      start_x--;
		      start_y--;
		      last_step=CENTRAL;
		  }
		  else if((WC[k][start_x][start_y]+
			   WR[k][start_x][start_y])/w > tmp){
		      /* take a path to the left */
		      start_y--;
		      last_step=RIGHT;
		  }
		  else {/* go up */
		      start_x--;
		      last_step=DOWN;		  
		  }
		  break;
	      }
	  }
        
	  /* print the score */
	  /*fprintf(B->IP->Datafiles->out_fpt,"k= %d , score= %3f\n",
		  Nblocks,score);*/
	  /* set profile */
	  for(k=1;k<=Nblocks;k++){
	      if(SK->start_x[0][k]!=0){/* find a block */
		  for(n=SK->start_x[0][k],j=SK->start_y[0][k];
		      n<=SK->end_x[0][k];n++,j++){
		      SK->profile[n][j]++;
		  }	
	      }
	  }
	  /* out put the matching string */
	  if(SK->flags.back_sampling_sequence){
	      for(k=1;k<=Nblocks;k++){
		  if(SK->start_x[0][k]!=0){
		      fprintf(B->IP->Datafiles->out_fpt,"%3d ",
			      SK->start_x[0][k]); 
		      for(n=SK->start_x[0][k];n<=SK->end_x[0][k];n++){
			  fprintf(B->IP->Datafiles->out_fpt,"%c",
				  (char )(SEQ2CH( SK->flags.protein_sequence,Seq->Orig[0][n-1])+65)); 
		      }
		      fprintf(B->IP->Datafiles->out_fpt,"% 3d\n",
			      SK->end_x[0][k]);
		      fprintf(B->IP->Datafiles->out_fpt,"%3d ",
			      SK->start_y[0][k]); 
		      for(n=SK->start_y[0][k];n<=SK->end_y[0][k];n++){
			  fprintf(B->IP->Datafiles->out_fpt,"%c",
				  (char )(SEQ2CH( SK->flags.protein_sequence,Seq->Orig[1][n-1])+65)); 
		      }
		      fprintf(B->IP->Datafiles->out_fpt,"% 3d\n",
			      SK->end_y[0][k]);
		      fprintf(B->IP->Datafiles->out_fpt,"\n");
		  }	
		  else{
		      fprintf(B->IP->Datafiles->out_fpt,"No matching.\n\n");
		  }
	      }
	  }
      }
  }

  /* normalize the profile */
  for(i=1;i<=Len0;i++){
      profile_ptr=SK->profile[i];
      for(j=1;j<=Len1;j++){
	  profile_ptr[j]=profile_ptr[j]/SK->backsampling_number;
      }
  }

  /* find the Bayesian alignment */
  Sankoff_Bayesian_align(B,SK);

  /* free the profile */
  for(i=0;i<=Len0;i++){
      free(SK->profile[i]);
  }
  free(SK->profile);

}



/*******************  Sankoff_Back_Sampling_Min ************************/
/* Minimum memory                                                 */
/* Back Sampling the alignment                                    */
/* The results are saved in SK->start_x[0][k], and                */
/* SK->end_x[0][k]                                                */
/*******************************************************************/

void 
Sankoff_Back_Sampling_Min(Model B, Sank_S SK)
{
   int             i, j, k, n, N;
   int             Nmatrix;
   int             Len0, Len1;	/* sequence lengths */
   double         *profile_ptr;	/* for quick access */
   float           score;
   int             Nblocks;
   int             strl;
   int            *Nback_sampling;	/* number of sample drawed for each
					 * matrix */
   int             sample_limit;   /* proportional sample limit */
   double          tmp, tmp1, w;
   double       ***WC, ***WD, ***WR;	/* for quick access */
   double        **wd_pptr, **wd_pptr0;	/* for quick access */
   double        **wc_pptr, **wc_pptr0;
   double        **wr_pptr, **wr_pptr0;
   double       ***pp_related;
   double        sample_factor;        /* scaling factor for random number */
   char          **matrix_name;
   short          *Orig0 = B->Seq->Orig[0];
   short          *Orig1 = B->Seq->Orig[1];
   float          *prior_N = SK->prior_N;
   enum
   {
      DOWN, CENTRAL, RIGHT
   }               last_step;
   Stype           Seq = B->Seq;
   save_loop      *current;
   save_loop      *get_next(save_loop *);
   save_loop      *alloc_min_Sank(int Nblocks, int Len0, int Len1, int last_step,int id);
   save_loop      *remove_back_sample(save_loop *top,save_loop *back_sample);
   void            print_link_list(save_loop *);
   void            Sankoff_probability_W(Model B, Sank_S SK);
   void            Sankoff_prob_RC(Model B, Sank_S SK, int level, int nMat);
   void            prob_RC(Model B, Sank_S SK, int level, int nMat);
   void            print_level(Model B, Sank_S SK, int level, int reallevel);
   int             level_change, level, id;
   int             work_blocks;
   int             pro_sum,block_sum;
   if (!SK->flags.min_memory) {
      fprintf(B->IP->Datafiles->out_fpt,"Program stopped. Error in Sankoff_Back_Sampling_Min. This routine should not be called unless flag is set\n");
      exit(1);
   }
   Len0 = B->Seq->SeqLen[0];
   Len1 = B->Seq->SeqLen[1];
   pp_related = SK->rel_matrix.pp_related;
   matrix_name = SK->rel_matrix.matrix_name;
   Nmatrix = SK->rel_matrix.Nmatrix;
   NEW(Nback_sampling, Nmatrix, int);
   /* allocate space for profile */
   NEWP(SK->profile, Len0 + 1, double);
   for (i = 0; i <= Len0; i++)
   { 
    NEW(SK->profile[i], Len1 + 1, double);
   }

   /* initialize the profile */
   for (i = 0; i <= Len0; i++)
   {
      profile_ptr = SK->profile[i];
      for (j = 0; j <= Len1; j++)
      {
	 profile_ptr[j] = 0;
      }
   }
    pro_sum=0;
    block_sum =0;

   /* set seed */
   sRandom((long) time(NULL));
   work_blocks = 2;


   fprintf(B->IP->Datafiles->out_fpt, "Back sampling.\n\n");
   if (Nmatrix > 1)
   {				/* multiple matrices, sampling according to
				 * P(psi|R) */
      for (N = 0; N < Nmatrix; N++)
      {
	 Nback_sampling[N] = (int) (SK->backsampling_number * SK->PSI[N]);
      }
   } else
   {
      Nback_sampling[0] = SK->backsampling_number;
   }
   /* init the id for debugging the backsamples */
    id=0;
   /* outer loop drives backsampling */

   for (N = Nmatrix - 1; N >= 0; N--)
   {
      if (N != Nmatrix - 1 && Nback_sampling[N] != 0)
      {
	 /* need to recalculate the WC, WR,WD */
	     Sankoff_prob_RC(B, SK, SK->Sankoff_blocks+1, N); 
      }
      WC = SK->WC;
      WR = SK->WR;
      WD = SK->WD;

      fprintf(B->IP->Datafiles->out_fpt,
	      "%d samples are from Matrix %s.\n",
	      Nback_sampling[N], matrix_name[N]);
      SK->top = NULL;		/* set the current task to the top */
      /* sample proportionally */  
      if(SK->backsampling_blocks > SK->Sankoff_blocks) 
             (SK->flags.back_sampling_proportional = TRUE); /* preserve old option */   
    
      for (i = 0; i < Nback_sampling[N]; i++)
      {
	 if (SK->flags.back_sampling_proportional)
	 {
	    /* sampling with random number of blocks according MagL */
	    /* get a random number between 0-1 */
	    tmp = Random() / (double) INT32_MAX;
	    /*    printf("Level 1 random number is %f\n",tmp); */
            if ( SK->backsampling_blocks > SK->Sankoff_blocks)
                sample_limit = SK->Sankoff_blocks;
            else
                sample_limit = SK->backsampling_blocks;
            
	    /* calculate factor to scale random number */
            if (sample_limit != SK->Sankoff_blocks)
              { 
                sample_factor = SK->PKSI[N][1];
                k = 1;
	        while (k < sample_limit )
                 {
                  k++;
                  sample_factor += SK->PKSI[N][k];
	          }

                 tmp = tmp * sample_factor;   /* scale the random number */
	      }
	    tmp1 = SK->PKSI[N][1];
	    k = 1;
	    while (k < sample_limit && tmp1 < tmp)
	    {
	       k++;
	       tmp1 += SK->PKSI[N][k];
	    }
	    Nblocks = k;
	 } else
	 {			/* sampling with fix number of blocks */
	    Nblocks = SK->backsampling_blocks;
	 }
	 /* allocate space to save partial results */
	 /* build task list */


	 if (SK->top == NULL)
	 {
	    SK->top = alloc_min_Sank(Nblocks, Len0, Len1, RIGHT,id++);
	    current = SK->top;
	 } else
	 {
	    current->next = alloc_min_Sank(Nblocks, Len0, Len1, RIGHT,id++);
	    current = current->next;
	 }

	 /*
	  * print_link_list(SK->top); *
	  */ 
	 } 
      /*   printf("Now processing %d backsamples\n",id); */
          /* init level */
      /*	  print_link_list(SK->top); */
	
         level = -10; /* force a recalculation of the probability matrix */
	 while (current = get_next(SK->top))
	 {
	    /* printf("current level=%d\n",level); */
	    if (level != current->k)
	    {
	      /* printf("need to recalculate probability matrix current level=%d new level=%d \n",level, current->k);*/


	         Sankoff_prob_RC(B, SK, current->k, N);
		 /* Sankoff_probability_W(B,SK); */
		 /*      print_level(B, SK, current->k % work_blocks, current->k);*/
	       WC = SK->WC;
	       WR = SK->WR;
	       WD = SK->WD;
	       level = current->k;
       	    }
	    level_change = FALSE;
	    while ((current->startb_x > 0 && current->startb_y > 0 && current->k > 0) &&  level_change != TRUE)
	    {
	       /* Get a random number between 0-1 */
	       tmp = Random() / (double) INT32_MAX;
	       /* printf("random number is %f\n",tmp); */
	       switch (current->last_step)
	       {
	       case DOWN:
		  w = WC[current->k % work_blocks][current->startb_x][current->startb_y] +
		     WD[current->k % work_blocks][current->startb_x][current->startb_y];
		  if (WC[current->k % work_blocks][current->startb_x][current->startb_y] / w > tmp)
		  {
		     /* take the center path */
		     /*
		      * score=score+ff_related[N]
		      * [B->Seq->Orig[0][current->startb_x-1]]
		      * [B->Seq->Orig[1][current->startb_y-1]];
		      */
		     current->end_x[0][current->k] = current->startb_x;
		     current->end_y[0][current->k] = current->startb_y;
		     current->start_x[0][current->k] = current->startb_x;
		     current->start_y[0][current->k] = current->startb_y;
		     current->startb_x--;
		     current->startb_y--;
		     current->last_step = CENTRAL;
		  } else
		  {		/* go up */
		     current->startb_x--;
		     current->last_step = DOWN;
		  }
		  break;

	       case CENTRAL:
		  w = WC[current->k % work_blocks][current->startb_x][current->startb_y] +
		     WR[(current->k - 1) % work_blocks][current->startb_x][current->startb_y] +
		     WD[(current->k - 1) % work_blocks][current->startb_x][current->startb_y];
		  if (WC[current->k % work_blocks][current->startb_x][current->startb_y] / w > tmp)
		  {
		     /* take the center path */
		     /*
		      * score=score+ff_related[N]
		      * [B->Seq->Orig[0][current->startb_x-1]]
		      * [B->Seq->Orig[1][current->startb_y-1]];
		      */
		     current->start_x[0][current->k] = current->startb_x;
		     current->start_y[0][current->k] = current->startb_y;
		     current->startb_x--;
		     current->startb_y--;
		     current->last_step = CENTRAL;
		  } else if ((WC[current->k % work_blocks][current->startb_x][current->startb_y] +
			      WR[(current->k - 1) % work_blocks][current->startb_x][current->startb_y]) / w > tmp)
		  {
		     /* take a path to go left */
		     current->startb_y--;
		     level = current->k;
		     current->k--;
		     level_change = TRUE;
		     current->last_step = RIGHT;
		  } else
		  {		/* go up */
		     current->startb_x--;
		     level = current->k;
		     current->k--;
		     level_change = TRUE;
		     current->last_step = DOWN;
		  }
		  break;

	       case RIGHT:
		  w = WC[current->k % work_blocks][current->startb_x][current->startb_y] +
		     WR[current->k % work_blocks][current->startb_x][current->startb_y] +
		     WD[current->k % work_blocks][current->startb_x][current->startb_y];
		  if (WC[current->k % work_blocks][current->startb_x][current->startb_y] / w > tmp)
		  {
		     /* take the center path */
		     /*
		      * score=score+ff_related[N]
		      * [B->Seq->Orig[0][current->startb_x-1]]
		      * [B->Seq->Orig[1][current->startb_y-1]];
		      */
		     current->end_x[0][current->k] = current->startb_x;
		     current->end_y[0][current->k] = current->startb_y;
		     current->start_x[0][current->k] = current->startb_x;
		     current->start_y[0][current->k] = current->startb_y;
		     current->startb_x--;
		     current->startb_y--;
		     current->last_step = CENTRAL;
		  } else if ((WC[current->k % work_blocks][current->startb_x][current->startb_y] +
			      WR[current->k % work_blocks][current->startb_x][current->startb_y]) / w > tmp)
		  {
		     /* take a path to the left */
		     current->startb_y--;
		     current->last_step = RIGHT;
		  } else
		  {		/* go up */
		     current->startb_x--;
		     current->last_step = DOWN;
		  }
		  break;
	       }
	    }
	    if (!(current->startb_x > 0) || !(current->startb_y > 0) || !(current->k > 0))
	    {
	       current->done = TRUE;
	       /* printf("completed backsample id=%d\n",current->id); */

	       /* print the score */
	       /*
	        * fprintf(B->IP->Datafiles->out_fpt,"k= %d , score= %3f\n",
	        * Nblocks,score);
	        */
	       /* set profile */
	       for (k = 1; k <= current->Nblocks; k++)
	       {
		  if (current->start_x[0][k] != 0)
		  {		/* find a block */
                      block_sum++;
		     for (n = current->start_x[0][k], j = current->start_y[0][k];
			  n <= current->end_x[0][k]; n++, j++)
		     {
			SK->profile[n][j]++;
                        pro_sum++;
		     }
		  }
	       }


	       /* print_link_list(SK->top); */
	       /* out put the matching string */
	       if (SK->flags.back_sampling_sequence)
	       {
		  /* printf("outputing sequences\n"); */
		  for (k = 1; k <= current->Nblocks; k++)
		  {
		     if (current->start_x[0][k] != 0)
		     {
			fprintf(B->IP->Datafiles->out_fpt, "%3d ",
				current->start_x[0][k]);
			for (n = current->start_x[0][k]; n <= current->end_x[0][k]; n++)
			{
			   fprintf(B->IP->Datafiles->out_fpt, "%c",
				 (char) (SEQ2CH( SK->flags.protein_sequence,Seq->Orig[0][n - 1]) + 65));
			}
			fprintf(B->IP->Datafiles->out_fpt, "% 3d\n",
				current->end_x[0][k]);
			fprintf(B->IP->Datafiles->out_fpt, "%3d ",
				current->start_y[0][k]);
			for (n = current->start_y[0][k]; n <= current->end_y[0][k]; n++)
			{
			   fprintf(B->IP->Datafiles->out_fpt, "%c",
				 (char) (SEQ2CH( SK->flags.protein_sequence,Seq->Orig[1][n - 1]) + 65));
			}
			fprintf(B->IP->Datafiles->out_fpt, "% 3d\n",
				current->end_y[0][k]);
			fprintf(B->IP->Datafiles->out_fpt, "\n");
		     } else
		     {
			fprintf(B->IP->Datafiles->out_fpt, "No matching.\n\n");
		     }
		  }
	       }
	       /* is done so remove it */
                 SK->top = remove_back_sample(SK->top,current);
	    } else
	    {
	       /* printf("level change to %d \n",current->k); */
	    }
	 }
      }
   /* printf("backsampling number is %d\n",SK->backsampling_number);
   printf("profile sum =%d    block sum=%d\n",pro_sum,block_sum); */
   /* normalize the profile */
   for (i = 1; i <= Len0; i++)
   {
      profile_ptr = SK->profile[i];
      for (j = 1; j <= Len1; j++)
      {
	 profile_ptr[j] = profile_ptr[j] / SK->backsampling_number;
	 /*    if (profile_ptr[j] > .000001)
            printf("i=%d,j=%d,p=%f\n",i,j, profile_ptr[j]); */
      }
   }

   /* find the Bayesian alignment */
   Sankoff_Bayesian_align(B, SK);

   /* free the profile */
   for (i = 0; i <= Len0; i++)
   {
      free(SK->profile[i]);
   }
   free(SK->profile);

}
/*******************remove_back_sample**************/
/* remove a back_sample from the list              */
/**************************************************/
save_loop * remove_back_sample(save_loop *top,save_loop *back_sample)
{
 save_loop *tmp_save,*tmp_top, *previous;
 char found;
 if(top==NULL || back_sample == NULL)
   return(NULL);
 found = FALSE;
 tmp_top = top;
while(tmp_top)
  {
     
    if (tmp_top == back_sample) {
        found= TRUE;
        break;
    }
   previous = tmp_top;
   tmp_top = tmp_top->next;
  }

  if (!found) {
       printf("error! In remove_back_sample cound not find backsample to remove\n");
       exit (1);
  }
  if (tmp_top == top) {
    /* should return storage */
   return(top->next);
  }
  else {
    previous->next = tmp_top->next;
    /* should return storage */
   return (top);
  }
}    
     
  
           
             

/******************get_next*************************/
/* */
/* get the next backsample at the maximum level    */
/***************************************************/
save_loop      *
get_next(save_loop * top)
{
   save_loop      *tmp_max;
   if (top == NULL)
     return(NULL);
   tmp_max = top;
   while (top != NULL)
   {
      if (tmp_max->k <= top->k)
	 tmp_max = top;
      if (top->k <= 0) {
        printf("Error in get_next k=%d \n",top->k);
        exit(1);
      }
      top = top->next;
   }
    return(tmp_max); /* it will return null if it is an empty list */
}

/*****************print_link_list*****************/
/* print the contents of the link list        */
/*************************************************/
void 
print_link_list(save_loop * top)
{
   while (top)
   {
      printf("last_step=%d, Nblocks=%d    K=%d    done=%d  startb_x=%d startb_y=%d, id=%d\n",
	     top->last_step, top->Nblocks, top->k, top->done, top->startb_x, top->startb_y,top->id);
      top = top->next;
   }
}


/************ alloc_min_Sank************************/
/*                                              */
/* Allocate the basic memory need by Sankoff   */
/* procedure or minimum memory mode.           */
/***********************************************/

save_loop      *
alloc_min_Sank(int Nblocks, int Len0, int Len1, int last_step, int id)
{
   int             i, j, k, n;
   save_loop      *tmp_save;

   NEW(tmp_save, 1, save_loop);
   /* allocate space for terminals of matching blocks */
   NEWP(tmp_save->start_x, Nblocks + 1, short);
   NEWP(tmp_save->start_y, Nblocks + 1, short);
   NEWP(tmp_save->end_x, Nblocks + 1, short);
   NEWP(tmp_save->end_y, Nblocks + 1, short);
   for (k = 0; k <= Nblocks; k++)
   {
      NEW(tmp_save->start_x[k], Nblocks + 1, short);
      NEW(tmp_save->start_y[k], Nblocks + 1, short);
      NEW(tmp_save->end_x[k], Nblocks + 1, short);
      NEW(tmp_save->end_y[k], Nblocks + 1, short);
   }
   /* initialize the terminal matching blocks */
   for (i = 0; i <= tmp_save->Nblocks; i++)
   {
      for (j = 0; j <= Nblocks; j++)
      {
	 tmp_save->start_x[i][j] = 0;
	 tmp_save->start_y[i][j] = 0;
	 tmp_save->end_x[i][j] = 0;
	 tmp_save->end_y[i][j] = 0;
      }
   }
   tmp_save->last_step = last_step;
   tmp_save->startb_x = Len0;
   tmp_save->startb_y = Len1;
   tmp_save->Nblocks = Nblocks;
   tmp_save->k = Nblocks;
   tmp_save->id = id;
   tmp_save->next = NULL;

   return (tmp_save);
}


/**************  Sankoff_prob_RC         ***************************/
/* */
/* This function calculates the sum of all probability,           */
/* then use it to calculate the posterior probabilty.             */
/* it is used for the minimum memory version                       */
/* */
/*******************************************************************/

void 
Sankoff_prob_RC(Model B, Sank_S SK, int level, int n_matrix)
{
   register int    i, i1, j, j1, k, n, N;
   int             Nmatrix;	/* number of matrices used  */
   double        **wd_pptr, **wr_pptr;	/* pointer for quick access */
   double        **wd_pptr0, **wr_pptr0;	/* pointer for quick access */
   double        **wc_pptr, **wc_pptr0;
   double         *wc_ptr, *wc_ptr1, *wd_ptr;
   double         *wd_ptr1, *wr_ptr, *wd_ptr01, *wr_ptr01;
   double          tmp;
   double       ***pp_related;
   short          *Orig0 = B->Seq->Orig[0];
   short          *Orig1 = B->Seq->Orig[1];
   int             work_blocks;	/* number of working matching blocks */
   int             Len0, Len1;
   int             maxfactor;	/* max factor of decrease            */
   long double     maxdecrease;	/* max of decrease                   */
   long double   **W;
   long double     sum, sumk;
   float          *prior_N = SK->prior_N;
   void            print_level(Model B, Sank_S SK, int level, int reallevel);
   Len0 = B->Seq->SeqLen[0];
   Len1 = B->Seq->SeqLen[1];
   pp_related = SK->rel_matrix.pp_related;

   
 

      /* only two working array are allocated */
      /* hard coded */
      if (SK->flags.min_memory)
	 work_blocks = 2;
      else
	 work_blocks = SK->Sankoff_blocks + 1;
      N = n_matrix;

      if (work_blocks != SK->Sankoff_blocks + 1)
      {
	 /* need to reset the boundary condition */
	 /* for k=0                              */

	 /* initialize the WC[0][i][j]=0 */
	 /* WD[0][i][j]=0, WR[0][i][j]=1 */
	 for (i = 1; i <= Len0; i++)
	 {
	    wc_ptr = SK->WC[0][i];
	    wd_ptr = SK->WD[0][i];
	    wr_ptr = SK->WR[0][i];
	    for (j = 1; j <= Len1; j++)
	    {
	       wc_ptr[j] = 0.0;
	       wd_ptr[j] = 0.0;
	       wr_ptr[j] = 1.0e-300;
	    }
	 }
	 /* initialize WC[0][i][0]=0 */
	 /* initialize WD[0][i][0]=1 */
	 /* initialize WR[0][i][0]=0 */
	 /* for i!=0                   */
	 for (i = 1; i <= Len0; i++)
	 {
	    SK->WC[0][i][0] = 0.0;
	    SK->WD[0][i][0] = 1.0e-300;
	    SK->WR[0][i][0] = 0.0;
	 }

	 /* initialize WC[0][0][j]=0 */
	 /* initialize WD[0][0][j]=0 */
	 /* initialize WR[0][0][j]=1 */
	 /* for any j                */
	 wc_ptr = SK->WC[0][0];
	 wd_ptr = SK->WD[0][0];
	 wr_ptr = SK->WR[0][0];
	 for (j = 0; j <= Len1; j++)
	 {
	    wc_ptr[j] = 0.0;
	    wd_ptr[j] = 0.0;
	    wr_ptr[j] = 1.0e-300;
	 }
      }
      /* need to recalculate the WC, WR,WD */
      for (k = 1; k <= level; k++)
      {
	 wd_pptr = SK->WD[k % work_blocks];
	 wd_pptr0 = SK->WD[(k - 1) % work_blocks];
	 wr_pptr = SK->WR[k % work_blocks];
	 wr_pptr0 = SK->WR[(k - 1) % work_blocks];
	 wc_pptr = SK->WC[k % work_blocks];
	 wc_pptr0 = SK->WC[(k - 1) % work_blocks];

	 /* keep track of decreasing factor */
	 SK->factor[N][k] = SK->factor[N][k - 1];

	 if (k == work_blocks)
	 {			/* reset the boundary condition */
	    for (i = 0; i <= Len0; i++)
	    {
	       wd_pptr[i][0] = 0.0;
	    }
	    for (j = 0; j <= Len1; j++)
	    {
	       wr_pptr[0][j] = 0.0;
	    }
	 }
	 tmp = wd_pptr0[Len0][Len1] + wr_pptr0[Len0][Len1] +
	    wc_pptr0[Len0][Len1];
	 if (tmp > 1.0e30)
	 {			/* decrease by a factor */
	    for (i = 1; i < Len0; i++)
	    {
	       for (j = 1; j <= Len1; j++)
	       {
		  wd_pptr0[i][j] = wd_pptr0[i][j] * 1.0e-300;
		  wr_pptr0[i][j] = wr_pptr0[i][j] * 1.0e-300;
		  wc_pptr0[i][j] = wc_pptr0[i][j] * 1.0e-300;
	       }
	    }
	    SK->factor[N][k]++;
	 }
	 for (i = 1; i <= Len0; i++)
	 {
	    for (j = 1; j <= Len1; j++)
	    {
	       wd_pptr[i][j] = wc_pptr[i - 1][j] + wd_pptr[i - 1][j];
	       wr_pptr[i][j] = wc_pptr[i][j - 1] + wd_pptr[i][j - 1] +
		  wr_pptr[i][j - 1];
	       wc_pptr[i][j] = (wd_pptr0[i - 1][j - 1] + wr_pptr0[i - 1][j - 1] +
				wc_pptr[i - 1][j - 1]) *
		  pp_related[N][Orig0[i - 1]][Orig1[j - 1]];
	    }
	 }
      }

      /* print_level(B,SK,level%work_blocks,level); */
      /* print_level(B,SK,(level-1)%work_blocks,level-1); */

   
}
/*********************print_level *********/
/* routine to print data structure for    */
/* minimum memory version                 */
/******************************************/ 
void 
print_level(Model B, Sank_S SK, int level, int reallevel)
{
   int             i, j, k, Len0, Len1,L0,L1,E0,E1;
   Len0 = B->Seq->SeqLen[0];
   L0=Len0;
   Len1 = B->Seq->SeqLen[1];
   L1= Len1;
   E0 = Len0-10;
   E1 = Len1-10;
   if (Len0 > 10)
      Len0 = 10;
   if (Len1 > 10)
      Len1 = 10;
   printf("level=%d  reallevel=%d\n", level, reallevel);
   for (i = 0; i <= Len0; i++)
   {
      for (j = 0; j <= Len1; j++)
      {
	 printf("i=%d,  j=%d  WC=%e  WD=%e  WR=%e\n", i, j,
	     SK->WC[level][i][j], SK->WD[level][i][j], SK->WR[level][i][j]);
      }
   }
  for (i = E0; i <= L0; i++)
   {
      for (j =E1; j <= L1; j++)
      {
	 printf("i=%d,  j=%d  WC=%e  WD=%e  WR=%e\n", i, j,
	     SK->WC[level][i][j], SK->WD[level][i][j], SK->WR[level][i][j]);
      }
   }
}

/*******************  Sankoff_Back_Sampling_cutoff *****************/
/* */
/* Back Sampling the alignment                                    */
/*******************************************************************/

void 
Sankoff_Back_Sampling_cutoff(Model B, Sank_S SK)
{

   int             i, k, n;
   int             Nblocks;
   double          sumS, Ln2 = log(2);
   enum Last_step
   {
      DOWN, CENTRAL, RIGHT
   };
   struct Sank_BS_Cutoff_node
   {
      enum Last_step  last_step;
      float           score;	/* score from the end to this node */
      int             k;	/* number of block from the start  */
      /* to this blocks                  */
      int             i, j;	/* current position                */
      int            *start_x, *start_y;	/* starts of matching
						 * segments     */
      int            *end_x, *end_y;	/* ends of matching segments       */
      struct Sank_BS_Cutoff_node *last;	/* address for the last one        */
   };
   typedef struct Sank_BS_Cutoff_node Sank_BS_Cutoff_node;
   Sank_BS_Cutoff_node *current_node, *node_tmp;
   struct
   {
      Sank_BS_Cutoff_node *header;
      int             num;
   }               node_list;

   short          *Orig0 = B->Seq->Orig[0];	/* for quick access                */
   short          *Orig1 = B->Seq->Orig[1];	/* for quick access                */
   float        ***ff_related = SK->rel_matrix.ff_related;

   fprintf(B->IP->Datafiles->out_fpt,
       "Back sampling above cutoff=%f.\n\n", SK->backsampling_cutoff_value);

   sumS = 0;
   Nblocks = SK->backsampling_blocks;

   /* set up the first  node */
   NEW(node_tmp, 1, Sank_BS_Cutoff_node);
   NEW(node_tmp->start_x, Nblocks + 1, int);
   NEW(node_tmp->start_y, Nblocks + 1, int);
   NEW(node_tmp->end_x, Nblocks + 1, int);
   NEW(node_tmp->end_y, Nblocks + 1, int);
   node_tmp->last_step = RIGHT;
   node_tmp->score = 0;
   node_tmp->k = Nblocks;
   node_tmp->i = B->Seq->SeqLen[0];
   node_tmp->j = B->Seq->SeqLen[1];
   node_tmp->last = NULL;

   node_list.header = node_tmp;
   node_list.num = 1;

   while (node_list.num > 0)
   {
      /* pop up a node */
      current_node = node_list.header;
      node_list.header = node_list.header->last;
      node_list.num--;

      if (current_node->score +
	  SK->SW[0][current_node->k][current_node->i][current_node->j] >=
	  SK->backsampling_cutoff_value)
      {				/* explore this node */
	 if (current_node->i > 0 && current_node->j > 0 &&
	     current_node->k > 0)
	 {			/* explore this node */
	    switch (current_node->last_step)
	    {
	    case DOWN:		/* two choices */
	       /* path to go up */
	       NEW(node_tmp, 1, Sank_BS_Cutoff_node);
	       NEW(node_tmp->start_x, Nblocks + 1, int);
	       NEW(node_tmp->start_y, Nblocks + 1, int);
	       NEW(node_tmp->end_x, Nblocks + 1, int);
	       NEW(node_tmp->end_y, Nblocks + 1, int);
	       memcpy(node_tmp->start_x, current_node->start_x,
		      sizeof(int) * (Nblocks + 1));
	       memcpy(node_tmp->start_y, current_node->start_y,
		      sizeof(int) * (Nblocks + 1));
	       memcpy(node_tmp->end_x, current_node->end_x,
		      sizeof(int) * (Nblocks + 1));
	       memcpy(node_tmp->end_y, current_node->end_y,
		      sizeof(int) * (Nblocks + 1));
	       node_tmp->last_step = DOWN;
	       node_tmp->score = current_node->score;
	       node_tmp->k = current_node->k;
	       node_tmp->i = current_node->i - 1;
	       node_tmp->j = current_node->j;
	       node_tmp->last = node_list.header;
	       node_list.header = node_tmp;
	       node_list.num++;

	       /* path central */
	       NEW(node_tmp, 1, Sank_BS_Cutoff_node);
	       NEW(node_tmp->start_x, Nblocks + 1, int);
	       NEW(node_tmp->start_y, Nblocks + 1, int);
	       NEW(node_tmp->end_x, Nblocks + 1, int);
	       NEW(node_tmp->end_y, Nblocks + 1, int);
	       memcpy(node_tmp->start_x, current_node->start_x,
		      sizeof(int) * (Nblocks + 1));
	       memcpy(node_tmp->start_y, current_node->start_y,
		      sizeof(int) * (Nblocks + 1));
	       memcpy(node_tmp->end_x, current_node->end_x,
		      sizeof(int) * (Nblocks + 1));
	       memcpy(node_tmp->end_y, current_node->end_y,
		      sizeof(int) * (Nblocks + 1));
	       node_tmp->last_step = CENTRAL;
	       node_tmp->score = current_node->score +
		  ff_related[0][Orig0[current_node->i - 1]]
		  [Orig1[current_node->j - 1]];
	       node_tmp->k = current_node->k;
	       node_tmp->i = current_node->i - 1;
	       node_tmp->j = current_node->j - 1;
	       node_tmp->end_x[current_node->k] = current_node->i;
	       node_tmp->end_y[current_node->k] = current_node->j;
	       node_tmp->start_x[current_node->k] = current_node->i;
	       node_tmp->start_y[current_node->k] = current_node->j;
	       node_tmp->last = node_list.header;
	       node_list.header = node_tmp;
	       node_list.num++;
	       break;

	    case CENTRAL:	/* three choice */
	       /* take the center path */
	       NEW(node_tmp, 1, Sank_BS_Cutoff_node);
	       NEW(node_tmp->start_x, Nblocks + 1, int);
	       NEW(node_tmp->start_y, Nblocks + 1, int);
	       NEW(node_tmp->end_x, Nblocks + 1, int);
	       NEW(node_tmp->end_y, Nblocks + 1, int);
	       memcpy(node_tmp->start_x, current_node->start_x,
		      sizeof(int) * (Nblocks + 1));
	       memcpy(node_tmp->start_y, current_node->start_y,
		      sizeof(int) * (Nblocks + 1));
	       memcpy(node_tmp->end_x, current_node->end_x,
		      sizeof(int) * (Nblocks + 1));
	       memcpy(node_tmp->end_y, current_node->end_y,
		      sizeof(int) * (Nblocks + 1));
	       node_tmp->last_step = CENTRAL;
	       node_tmp->score = current_node->score +
		  ff_related[0][Orig0[current_node->i - 1]]
		  [Orig1[current_node->j - 1]];
	       node_tmp->k = current_node->k;
	       node_tmp->i = current_node->i - 1;
	       node_tmp->j = current_node->j - 1;
	       node_tmp->start_x[current_node->k] = current_node->i;
	       node_tmp->start_y[current_node->k] = current_node->j;
	       node_tmp->last = node_list.header;
	       node_list.header = node_tmp;
	       node_list.num++;

	       /* take a path to go left */
	       NEW(node_tmp, 1, Sank_BS_Cutoff_node);
	       NEW(node_tmp->start_x, Nblocks + 1, int);
	       NEW(node_tmp->start_y, Nblocks + 1, int);
	       NEW(node_tmp->end_x, Nblocks + 1, int);
	       NEW(node_tmp->end_y, Nblocks + 1, int);
	       memcpy(node_tmp->start_x, current_node->start_x,
		      sizeof(int) * (Nblocks + 1));
	       memcpy(node_tmp->start_y, current_node->start_y,
		      sizeof(int) * (Nblocks + 1));
	       memcpy(node_tmp->end_x, current_node->end_x,
		      sizeof(int) * (Nblocks + 1));
	       memcpy(node_tmp->end_y, current_node->end_y,
		      sizeof(int) * (Nblocks + 1));
	       node_tmp->last_step = RIGHT;
	       node_tmp->score = current_node->score;
	       node_tmp->k = current_node->k - 1;
	       node_tmp->i = current_node->i;
	       node_tmp->j = current_node->j - 1;
	       node_tmp->last = node_list.header;
	       node_list.header = node_tmp;
	       node_list.num++;

	       /* go up */
	       NEW(node_tmp, 1, Sank_BS_Cutoff_node);
	       NEW(node_tmp->start_x, Nblocks + 1, int);
	       NEW(node_tmp->start_y, Nblocks + 1, int);
	       NEW(node_tmp->end_x, Nblocks + 1, int);
	       NEW(node_tmp->end_y, Nblocks + 1, int);
	       memcpy(node_tmp->start_x, current_node->start_x,
		      sizeof(int) * (Nblocks + 1));
	       memcpy(node_tmp->start_y, current_node->start_y,
		      sizeof(int) * (Nblocks + 1));
	       memcpy(node_tmp->end_x, current_node->end_x,
		      sizeof(int) * (Nblocks + 1));
	       memcpy(node_tmp->end_y, current_node->end_y,
		      sizeof(int) * (Nblocks + 1));
	       node_tmp->last_step = DOWN;
	       node_tmp->score = current_node->score;
	       node_tmp->k = current_node->k - 1;
	       node_tmp->i = current_node->i - 1;
	       node_tmp->j = current_node->j;
	       node_tmp->last = node_list.header;
	       node_list.header = node_tmp;
	       node_list.num++;
	       break;

	    case RIGHT:	/* three choice */
	       /* take the center path */
	       NEW(node_tmp, 1, Sank_BS_Cutoff_node);
	       NEW(node_tmp->start_x, Nblocks + 1, int);
	       NEW(node_tmp->start_y, Nblocks + 1, int);
	       NEW(node_tmp->end_x, Nblocks + 1, int);
	       NEW(node_tmp->end_y, Nblocks + 1, int);
	       memcpy(node_tmp->start_x, current_node->start_x,
		      sizeof(int) * (Nblocks + 1));
	       memcpy(node_tmp->start_y, current_node->start_y,
		      sizeof(int) * (Nblocks + 1));
	       memcpy(node_tmp->end_x, current_node->end_x,
		      sizeof(int) * (Nblocks + 1));
	       memcpy(node_tmp->end_y, current_node->end_y,
		      sizeof(int) * (Nblocks + 1));
	       node_tmp->last_step = CENTRAL;
	       node_tmp->score = current_node->score +
		  ff_related[0][Orig0[current_node->i - 1]]
		  [Orig1[current_node->j - 1]];
	       node_tmp->k = current_node->k;
	       node_tmp->i = current_node->i - 1;
	       node_tmp->j = current_node->j - 1;
	       node_tmp->end_x[current_node->k] = current_node->i;
	       node_tmp->end_y[current_node->k] = current_node->j;
	       node_tmp->start_x[current_node->k] = current_node->i;
	       node_tmp->start_y[current_node->k] = current_node->j;
	       node_tmp->last = node_list.header;
	       node_list.header = node_tmp;
	       node_list.num++;

	       /* take a path to the left */
	       NEW(node_tmp, 1, Sank_BS_Cutoff_node);
	       NEW(node_tmp->start_x, Nblocks + 1, int);
	       NEW(node_tmp->start_y, Nblocks + 1, int);
	       NEW(node_tmp->end_x, Nblocks + 1, int);
	       NEW(node_tmp->end_y, Nblocks + 1, int);
	       memcpy(node_tmp->start_x, current_node->start_x,
		      sizeof(int) * (Nblocks + 1));
	       memcpy(node_tmp->start_y, current_node->start_y,
		      sizeof(int) * (Nblocks + 1));
	       memcpy(node_tmp->end_x, current_node->end_x,
		      sizeof(int) * (Nblocks + 1));
	       memcpy(node_tmp->end_y, current_node->end_y,
		      sizeof(int) * (Nblocks + 1));
	       node_tmp->last_step = RIGHT;
	       node_tmp->score = current_node->score;
	       node_tmp->k = current_node->k;
	       node_tmp->i = current_node->i;
	       node_tmp->j = current_node->j - 1;
	       node_tmp->last = node_list.header;
	       node_list.header = node_tmp;
	       node_list.num++;

	       /* go up */
	       NEW(node_tmp, 1, Sank_BS_Cutoff_node);
	       NEW(node_tmp->start_x, Nblocks + 1, int);
	       NEW(node_tmp->start_y, Nblocks + 1, int);
	       NEW(node_tmp->end_x, Nblocks + 1, int);
	       NEW(node_tmp->end_y, Nblocks + 1, int);
	       memcpy(node_tmp->start_x, current_node->start_x,
		      sizeof(int) * (Nblocks + 1));
	       memcpy(node_tmp->start_y, current_node->start_y,
		      sizeof(int) * (Nblocks + 1));
	       memcpy(node_tmp->end_x, current_node->end_x,
		      sizeof(int) * (Nblocks + 1));
	       memcpy(node_tmp->end_y, current_node->end_y,
		      sizeof(int) * (Nblocks + 1));
	       node_tmp->last_step = DOWN;
	       node_tmp->score = current_node->score;
	       node_tmp->k = current_node->k;
	       node_tmp->i = current_node->i - 1;
	       node_tmp->j = current_node->j;
	       node_tmp->last = node_list.header;
	       node_list.header = node_tmp;
	       node_list.num++;
	       break;
	    }			/* end of switch block */
	 } else
	 {			/* not explore this node */

	    /* make sure do not report same alignment multiple time */
	    if (current_node->score >= SK->backsampling_cutoff_value &&
		(current_node->k == 0 && current_node->last_step == RIGHT))
	    {			/* report the alignment */
	       fprintf(B->IP->Datafiles->out_fpt, "score=% 3f\n",
		       current_node->score);
	       sumS = sumS + exp(current_node->score / 2. * Ln2);
	       /* out put the matching string */
	       if (SK->flags.back_sampling_sequence)
	       {
		  for (k = 1; k <= Nblocks; k++)
		  {
		     if (current_node->start_x[k] != 0)
		     {
			fprintf(B->IP->Datafiles->out_fpt, "%3d ",
				current_node->start_x[k]);
			for (n = current_node->start_x[k];
			     n <= current_node->end_x[k]; n++)
			{
			   fprintf(B->IP->Datafiles->out_fpt, "%c",
				   (char) (SEQ2CH( SK->flags.protein_sequence,Orig0[n - 1]) + 65));
			}
			fprintf(B->IP->Datafiles->out_fpt, "% 3d\n",
				current_node->end_x[k]);
			fprintf(B->IP->Datafiles->out_fpt, "%3d ",
				current_node->start_y[k]);
			for (n = current_node->start_y[k];
			     n <= current_node->end_y[k]; n++)
			{
			   fprintf(B->IP->Datafiles->out_fpt, "%c",
				   (char) (SEQ2CH( SK->flags.protein_sequence,Orig1[n - 1]) + 65));
			}
			fprintf(B->IP->Datafiles->out_fpt, "% 3d\n",
				current_node->end_y[k]);
			fprintf(B->IP->Datafiles->out_fpt, "\n");
		     }
		  }		/* end of for block */
	       }		/* end of
				 * if(SK->flags.back_sampling_sequence) */
	    }			/* end of if(current_node->score
				 * >=SK->backsampling_cutoff) */
	 }			/* end of else block */

      }
      /* release the space for current node */
      free(current_node->end_x);
      free(current_node->end_y);
      free(current_node->start_x);
      free(current_node->start_y);
      free(current_node);
   }				/* end of while node */

   fprintf(B->IP->Datafiles->out_fpt, "SumS(S>=%4f,k=%d)/N[%d]=%g  ",
	   SK->backsampling_cutoff_value,
	   SK->backsampling_blocks,
	   SK->backsampling_blocks,
	   sumS / SK->N[SK->backsampling_blocks]);

   fprintf(B->IP->Datafiles->out_fpt, "SumS(S>=%4f,k=%d)/W[%d]=%g\n",
	   SK->backsampling_cutoff_value,
	   SK->backsampling_blocks,
	   SK->backsampling_blocks,
	   sumS / SK->W[0][SK->backsampling_blocks]);

}


/******************* Sankoff_Make_S  ********************************/
/* */
/* fill in score matrix with specific blocks                    */
/********************************************************************/

void 
Sankoff_Make_S(Stype Seq, Sank_S SK)
{
   register int    i, j, k, n, m;	/* general registers   */
   int             L0, L1;	/* length of sequences */
   short          *Orig0, *Orig1;	/* for quick access    */
   float         **sv_pptr, **sv_pptr1;
   float         **sw_pptr, **sw_pptr1;
   float        ***ff_related = SK->rel_matrix.ff_related;
   int             Nmatrix;

   Nmatrix = SK->rel_matrix.Nmatrix;
   L0 = Seq->SeqLen[0];
   L1 = Seq->SeqLen[1];
   Orig0 = Seq->Orig[0];
   Orig1 = Seq->Orig[1];

   for (n = 0; n < Nmatrix; n++)
   {
      for (k = 1; k <= SK->Sankoff_blocks; k++)
      {
	 sv_pptr = SK->SV[n][k];
	 sv_pptr1 = SK->SV[n][k - 1];
	 sw_pptr = SK->SW[n][k];
	 sw_pptr1 = SK->SW[n][k - 1];
	 for (i = 1; i <= L0; i++)
	 {
	    for (j = 1; j <= L1; j++)
	    {
	       sv_pptr[i][j] = ff_related[n]
		  [Orig0[i - 1]][Orig1[j - 1]] +
		  max3(sv_pptr[i - 1][j - 1], sw_pptr1[i - 1][j - 1], -99);
	       sw_pptr[i][j] = max3(sw_pptr[i - 1][j], sv_pptr[i][j],
				    sw_pptr[i][j - 1]);
	    }
	 }
      }
   }
}

/*******************  Sankoff_Trace_S  ********************************/
/* */
/* trace back score matrix with specific blocks                  */
/**********************************************************************/

void 
Sankoff_Trace_S(Model B, Sank_S SK)
{
   int             i, j, k, K, m, n, N, t;
   Stype           Seq = B->Seq;
   float         **sv_pptr, **sw_pptr, **sw_pptr1;
   float           tmp, tmp_max;
   int             Nmatrix;
   int             block_start;
   double          Ln2 = log(2);
   double          posterior, map, margin;
   float        ***ff_related = SK->rel_matrix.ff_related;
   char          **matrix_name = SK->rel_matrix.matrix_name;

   Nmatrix = SK->rel_matrix.Nmatrix;
   for (N = 0; N < Nmatrix; N++)
   {
      /* output matrix name */
      fprintf(B->IP->Datafiles->out_fpt,
	      "%s\n\n", matrix_name[N]);

      for (K = 1; K <= SK->Sankoff_blocks; K++)
      {
	 k = 0;
	 t = 1;
	 sw_pptr1 = SK->SW[N][K - k - 1];
	 sw_pptr = SK->SW[N][K - k];
	 sv_pptr = SK->SV[N][K - k];

	 i = Seq->SeqLen[0];
	 j = Seq->SeqLen[1];
	 block_start = TRUE;

	 tmp_max = sw_pptr[i][j];
	 fprintf(B->IP->Datafiles->out_fpt,
		 "k=%d, Max score=%g\n", K, tmp_max);

	 while (i >= 1 && j >= 1 && k <= K - 1)
	 {
	    if (block_start)
	    {			/* could be start of a block */
	       if (ff_related[N]
		   [Seq->Orig[0][i - 1]][Seq->Orig[1][j - 1]] >= 0)
	       {		/* it is possible start of a block */
		  if (sw_pptr[i][j] > sw_pptr[i - 1][j] &&
		      sw_pptr[i][j] > sw_pptr[i][j - 1])
		     /* because the error over floating addition */
		     /* we have to be careful of comparison     */
		  {
		     /* start of an alignment block */
		     SK->end_x[K][t] = i;
		     SK->end_y[K][t] = j;
		     block_start = FALSE;
		     if (sw_pptr1[i - 1][j - 1] > sv_pptr[i - 1][j - 1])
		     {
			/* also end of an alignment block */
			SK->start_x[K][t] = i;
			SK->start_y[K][t] = j;
			block_start = TRUE;
			t++;
			/* check the next matrix */
			k++;
			sw_pptr1 = SK->SW[N][K - k - 1];
			sw_pptr = SK->SW[N][K - k];
			sv_pptr = SK->SV[N][K - k];
		     }
		     i--;
		     j--;
		  } else
		  {		/* non-negative, but could not start of */
		     /* an alignment                         */
		     if (sw_pptr[i - 1][j] > sw_pptr[i][j - 1])
			i--;
		     else
			j--;
		  }
	       } else
	       {		/* negative, could not be a start of */
		  /* an alignment                      */
		  if (sw_pptr[i - 1][j] > sw_pptr[i][j - 1])
		     i--;
		  else
		     j--;
	       }
	    } else
	    {			/* extend or end a block */
	       if (ff_related[N]
		   [Seq->Orig[0][i - 1]][Seq->Orig[1][j - 1]] >= 0)
	       {		/* could be the end of a block */
		  if (sw_pptr1[i - 1][j - 1] > sv_pptr[i - 1][j - 1])
		  {
		     /* end of an alignment block */
		     SK->start_x[K][t] = i;
		     SK->start_y[K][t] = j;
		     block_start = TRUE;
		     t++;
		     /* check the next matrix */
		     k++;
		     sw_pptr1 = SK->SW[N][K - k - 1];
		     sw_pptr = SK->SW[N][K - k];
		     sv_pptr = SK->SV[N][K - k];
		  }
		  i--;
		  j--;
	       } else
	       {		/* must be extension of a block */
		  i--;
		  j--;
	       }
	    }
	 }
	 /* output the matching string */
	 if (SK->flags.output_sequence)
	 {
	    for (k = K; k >= 1; k--)
	    {
	       if (SK->start_x[K][k] != 0)
	       {
		  fprintf(B->IP->Datafiles->out_fpt, "%3d ",
			  SK->start_x[K][k]);
		  for (n = SK->start_x[K][k]; n <= SK->end_x[K][k]; n++)
		  {
		     fprintf(B->IP->Datafiles->out_fpt, "%c",
			     (char) (SEQ2CH( SK->flags.protein_sequence,Seq->Orig[0][n - 1]) + 65));
		  }
		  fprintf(B->IP->Datafiles->out_fpt, "% 3d\n",
			  SK->end_x[K][k]);
		  fprintf(B->IP->Datafiles->out_fpt, "%3d ",
			  SK->start_y[K][k]);
		  for (n = SK->start_y[K][k]; n <= SK->end_y[K][k]; n++)
		  {
		     fprintf(B->IP->Datafiles->out_fpt, "%c",
			     (char) (SEQ2CH( SK->flags.protein_sequence,Seq->Orig[1][n - 1]) + 65));
		  }
		  fprintf(B->IP->Datafiles->out_fpt, "% 3d\n",
			  SK->end_y[K][k]);
		  fprintf(B->IP->Datafiles->out_fpt, "\n");
	       }
	    }
	 }
	 /* posterior probability */
	 map = exp(SK->SW[N][K][Seq->SeqLen[0]][Seq->SeqLen[1]] / 2. * Ln2) /
	    SK->N[K];

	 margin = SK->W[N][K] / SK->N[K];
	 posterior = map / margin;
	 fprintf(B->IP->Datafiles->out_fpt,
		 "k=%d, W[k]=%g, N[k]=%g\n"
		 "S[k]/N[k]=%g, W[k]/N[k]=%g, S[k]/W[k]=%g\n",
		 K, SK->W[N][K], SK->N[K], map, margin, posterior);

	 /* posterior of k>0 */
	 posterior = 0;
	 for (i = 1; i <= K; i++)
	    posterior += SK->W[N][i] / SK->N[i];
	 posterior = posterior / (float) K;
	 fprintf(B->IP->Datafiles->out_fpt,
		 "P(k>0|R)=%g\n", posterior / (1. + posterior));

	 /* score difference and expect score difference */
	 fprintf(B->IP->Datafiles->out_fpt,
		 "score difference=%3f, expected score difference=%f\n\n",
		 SK->SW[N][K][Seq->SeqLen[0]][Seq->SeqLen[1]] -
		 SK->SW[N][K - 1][Seq->SeqLen[0]][Seq->SeqLen[1]],
		 SC(B, SK, K));
      }
   }
}

/**************  Sankoff_probability_N   ***************************/
/* */
/* This function calculates total number of alignment.            */
/* */
/*******************************************************************/

void 
Sankoff_probability_N(Model B, Sank_S SK)
{
   register int    i, j, k, n;
   double        **nd_pptr, **nr_pptr;	/* pointer for quick access */
   double        **nd_pptr0, **nr_pptr0;	/* pointer for quick access */
   double        **nc_pptr;
   int             Len0, Len1;
   short          *Orig0 = B->Seq->Orig[0];
   short          *Orig1 = B->Seq->Orig[1];
   float          *prior_N = SK->prior_N;

   Len0 = B->Seq->SeqLen[0];
   Len1 = B->Seq->SeqLen[1];


   /* calculating N */
   /* value for no matching blocks */
   SK->N[0] = 1;

   for (k = 1; k <= SK->Sankoff_blocks; k++)
   {
      nd_pptr = SK->ND[k % 2];
      nd_pptr0 = SK->ND[(k - 1) % 2];
      nr_pptr = SK->NR[k % 2];
      nr_pptr0 = SK->NR[(k - 1) % 2];
      nc_pptr = SK->NC[k % 2];

      if (k == 2)
      {				/* need to reset the boundary condition */
	 for (i = 0; i <= Len0; i++)
	 {
	    nd_pptr[i][0] = 0;
	 }
	 for (j = 0; j <= Len1; j++)
	 {
	    nr_pptr[0][j] = 0;
	 }
      }
      for (i = 1; i <= Len0; i++)
      {
	 /*
	  * printf ("%4.2f ",prior_N[i]); if(i%10==0) printf ("\n");
	  */
	 for (j = 1; j <= Len1; j++)
	 {
	    nd_pptr[i][j] = nc_pptr[i - 1][j] + nd_pptr[i - 1][j];
	    nr_pptr[i][j] = nc_pptr[i][j - 1] + nd_pptr[i][j - 1] +
	       nr_pptr[i][j - 1];
	    nc_pptr[i][j] = (nd_pptr0[i - 1][j - 1] + nr_pptr0[i - 1][j - 1] +
			     nc_pptr[i - 1][j - 1]) * prior_N[i];
	 }
      }
      SK->N[k] = nd_pptr[Len0][Len1] +
	 nr_pptr[Len0][Len1] +
	 nc_pptr[Len0][Len1];
   }
}

/**************  Sankoff_probability_W   ***************************/
/* */
/* This function calculates the sum of all probability,           */
/* then use it to calculate the posterior probabilty.             */
/* */
/*******************************************************************/

void 
Sankoff_probability_W(Model B, Sank_S SK)
{
   register int    i, i1, j, j1, k, n;
   int             Nmatrix;	/* number of matrices used  */
   double        **wd_pptr, **wr_pptr;	/* pointer for quick access */
   double        **wd_pptr0, **wr_pptr0;	/* pointer for quick access */
   double        **wc_pptr, **wc_pptr0;
   double         *wc_ptr, *wc_ptr1, *wd_ptr;
   double         *wd_ptr1, *wr_ptr, *wd_ptr01, *wr_ptr01;
   double          tmp;
   double       ***pp_related;
   short          *Orig0 = B->Seq->Orig[0];
   short          *Orig1 = B->Seq->Orig[1];
   int             work_blocks;	/* number of working matching blocks */
   int             Len0, Len1;
   int             maxfactor;	/* max factor of decrease            */
   long double     maxdecrease;	/* max of decrease                   */
   long double   **W;
   long double     sum, sumk;
   float          *prior_N = SK->prior_N;

   Len0 = B->Seq->SeqLen[0];
   Len1 = B->Seq->SeqLen[1];
   pp_related = SK->rel_matrix.pp_related;

   Nmatrix = SK->rel_matrix.Nmatrix;

   if (SK->flags.back_sampling || SK->flags.exact_posterior)
   {
      /* allocated full space */
      work_blocks = SK->Sankoff_blocks + 1;
   } else
   {				/* only two working array are allocated */
      work_blocks = 2;
   }
   /*  set up minimum working space */
   if (SK->flags.min_memory)
      work_blocks = 2;
   else
      work_blocks = SK->Sankoff_blocks + 1;

   /* value for no matching blocks */
   for (n = 0; n < Nmatrix; n++)
      SK->W[n][0] = 1.0e-300;

   for (n = 0; n < Nmatrix; n++)
   {
      if (n != 0 && work_blocks != SK->Sankoff_blocks + 1)
      {
	 /* need to reset the boundary condition */
	 /* for k=0                              */

	 /* initialize the WC[0][i][j]=0 */
	 /* WD[0][i][j]=0, WR[0][i][j]=1 */
	 for (i = 1; i <= Len0; i++)
	 {
	    wc_ptr = SK->WC[0][i];
	    wd_ptr = SK->WD[0][i];
	    wr_ptr = SK->WR[0][i];
	    for (j = 1; j <= Len1; j++)
	    {
	       wc_ptr[j] = 0.0;
	       wd_ptr[j] = 0.0;
	       wr_ptr[j] = 1.0e-300;
	    }
	 }
	 /* initialize WC[0][i][0]=0 */
	 /* initialize WD[0][i][0]=1 */
	 /* initialize WR[0][i][0]=0 */
	 /* for i!=0                   */
	 for (i = 1; i <= Len0; i++)
	 {
	    SK->WC[0][i][0] = 0.0;
	    SK->WD[0][i][0] = 1.0e-300;
	    SK->WR[0][i][0] = 0.0;
	 }

	 /* initialize WC[0][0][j]=0 */
	 /* initialize WD[0][0][j]=0 */
	 /* initialize WR[0][0][j]=1 */
	 /* for any j                */
	 wc_ptr = SK->WC[0][0];
	 wd_ptr = SK->WD[0][0];
	 wr_ptr = SK->WR[0][0];
	 for (j = 0; j <= Len1; j++)
	 {
	    wc_ptr[j] = 0.0;
	    wd_ptr[j] = 0.0;
	    wr_ptr[j] = 1.0e-300;
	 }
      }
      /* calculate the sum of probability */
      for (k = 1; k <= SK->Sankoff_blocks; k++)
      {
	 wd_pptr = SK->WD[k % work_blocks];
	 wd_pptr0 = SK->WD[(k - 1) % work_blocks];
	 wr_pptr = SK->WR[k % work_blocks];
	 wr_pptr0 = SK->WR[(k - 1) % work_blocks];
	 wc_pptr = SK->WC[k % work_blocks];
	 wc_pptr0 = SK->WC[(k - 1) % work_blocks];

	 /* keep track of decreasing factor */
	 SK->factor[n][k] = SK->factor[n][k - 1];

	 if (k == work_blocks)
	 {			/* reset the boundary condition */
	    for (i = 0; i <= Len0; i++)
	    {
	       wd_pptr[i][0] = 0.0;
	    }
	    for (j = 0; j <= Len1; j++)
	    {
	       wr_pptr[0][j] = 0.0;
	    }
	 }
	 tmp = wd_pptr0[Len0][Len1] + wr_pptr0[Len0][Len1] +
	    wc_pptr0[Len0][Len1];
	 if (tmp > 1.0e30)
	 {			/* decrease by a factor */
	    for (i = 0; i < Len0; i++)
	    {
	       for (j = 0; j <= Len1; j++)
	       {
		  wd_pptr0[i][j] = wd_pptr0[i][j] * 1.0e-300;
		  wr_pptr0[i][j] = wr_pptr0[i][j] * 1.0e-300;
		  wc_pptr0[i][j] = wc_pptr0[i][j] * 1.0e-300;
	       }
	    }
	    SK->factor[n][k]++;
	 }

	 for (i = 1, i1 = 0; i <= Len0; i++, i1++)
	 {
	    wc_ptr = wc_pptr[i];
	    wc_ptr1 = wc_pptr[i1];
	    wd_ptr = wd_pptr[i];
	    wd_ptr1 = wd_pptr[i1];
	    wr_ptr = wr_pptr[i];
	    wd_ptr01 = wd_pptr0[i1];
	    wr_ptr01 = wr_pptr0[i1];

	    for (j = 1, j1 = 0; j <= Len1; j++, j1++)
	    {
	       /* wd_pptr[i][j]=wc_pptr[i-1][j]+wd_pptr[i-1][j]; */
	       wd_ptr[j] = wc_ptr1[j] + wd_ptr1[j];
	       /*
	        * wr_pptr[i][j]=wc_pptr[i][j-1]+wd_pptr[i][j-1]+
	        * wr_pptr[i][j-1];
	        */
	       wr_ptr[j] = wc_ptr[j1] + wd_ptr[j1] +
		  wr_ptr[j1];
	       /*
	        * wc_pptr[i][j]=(wd_pptr0[i-1][j-1]+wr_pptr0[i-1][j-1]+
	        * wc_pptr[i-1][j-1])*
	        * p_related[matrix_num][Orig0[i-1]][Orig1[j-1]]* prior_N[i];
	        */
	       wc_ptr[j] = (wd_ptr01[j1] + wr_ptr01[j1] +
			    wc_ptr1[j1]) *
		  pp_related[n][Orig0[i1]][Orig1[j1]];

	       /******************************debug*********************/
	       /*    printf("k=%d,i=%d,j=%d,wd=%g,wr=%g,wc=%g,matrix=%g  %c %c \n ",
                  k-1,i,j,wd_ptr[j],wr_ptr[j],wc_ptr[j],
                   pp_related[n][Orig0[i1]][Orig1[j1]],
                   (char) (SEQ2CH( SK->flags.protein_sequence,
                   Orig0[i1])+ 65),
                   (char) (SEQ2CH( SK->flags.protein_sequence,
                   Orig1[j1])+65)  );*/
             
	       /********************************************************/
	    }
	 }

	 SK->W[n][k] = wd_pptr[Len0][Len1] +
	    wr_pptr[Len0][Len1] +
	    wc_pptr[Len0][Len1];
      }
   }

   /* find out the max of decrease factor */
   maxfactor = 0;
   for (n = 0; n < Nmatrix; n++)
   {
      if (SK->factor[n][SK->Sankoff_blocks] > maxfactor)
      {
	 maxfactor = SK->factor[n][SK->Sankoff_blocks];
      }
   }

   /* allocate temporary spaces */
   /* SK->W still hold the exact sum, W will hold rescaled sum */
   NEWP(W, Nmatrix, long double);
   for (n = 0; n < Nmatrix; n++)
   {
      NEW(W[n], SK->Sankoff_blocks + 1, long double);
   }

   /* let the W[k] decrease the same amount of factor */
   for (n = 0; n < Nmatrix; n++)
   {
      for (k = 0; k <= SK->Sankoff_blocks; k++)
      {
	 W[n][k] = SK->W[n][k];
	 for (i = 0; i < (maxfactor - SK->factor[n][k]); i++)
	 {
	    W[n][k] = W[n][k] * 1.0e-300;
	 }
      }
   }

   /* constant decrease factor */
   maxdecrease = 1;
   for (i = 0; i <= maxfactor; i++)
   {
      maxdecrease = maxdecrease * 1.0e-300;
   }

   /* calculate PK,PSI,PKSI */
   for (n = 0; n < Nmatrix; n++)
   {
      sum = 0;
      for (k = 1; k <= SK->Sankoff_blocks; k++)
      {
	 sum += (W[n][k] / SK->N[k]);
      }
      sum = sum * 0.5 / (SK->Sankoff_blocks) + 0.5 * maxdecrease;
      SK->PM = sum;
      SK->PKSI[n][0] = 0.5 * maxdecrease / sum;
      for (k = 1; k <= SK->Sankoff_blocks; k++)
      {
	 SK->PKSI[n][k] = (W[n][k] / SK->N[k]) * 0.5 /
	    (SK->Sankoff_blocks * sum);
      }
   }

   if (Nmatrix > 1)
   {				/* multiple matrices */
      sum = 0;
      for (n = 0; n < Nmatrix; n++)
      {
	 for (k = 1; k <= SK->Sankoff_blocks; k++)
	 {
	    sum += (W[n][k] / SK->N[k]);
	 }
      }
      sum = sum * 0.5 / (Nmatrix * SK->Sankoff_blocks) + 0.5 * maxdecrease;
      SK->PM = sum;
      /* calculate PK */
      SK->PK[0] = 0.5 * maxdecrease / sum;
      for (k = 1; k <= SK->Sankoff_blocks; k++)
      {
	 sumk = 0;
	 for (n = 0; n < Nmatrix; n++)
	 {
	    sumk += W[n][k] / SK->N[k];
	 }
	 sumk = sumk * 0.5 / (Nmatrix * SK->Sankoff_blocks);
	 SK->PK[k] = sumk / (sum);
      }

      /* calculate PSI */
      sum = 0;
      for (n = 0; n < Nmatrix; n++)
      {
	 for (k = 1; k <= SK->Sankoff_blocks; k++)
	 {
	    sum += (W[n][k] / SK->N[k]);
	 }
      }
      for (n = 0; n < Nmatrix; n++)
      {
	 sumk = 0;
	 for (k = 1; k <= SK->Sankoff_blocks; k++)
	 {
	    sumk += (W[n][k] / SK->N[k]);
	 }
	 SK->PSI[n] = sumk / sum;
      }
   } else
   {				/* single matrix */
      for (k = 0; k <= SK->Sankoff_blocks; k++)
      {
	 SK->PK[k] = SK->PKSI[0][k];
      }
   }

   /* release space for W */
   for (n = 0; n < Nmatrix; n++)
   {
      free(W[n]);
   }
   free(W);
    /* print_prob(SK) */ ;
}

/**************  Sankoff_probability_BW   **************************/
/* */
/* This function calculates the sum of all probability backward,  */
/* then use it to calculate the posterior probabilty.             */
/* */
/*******************************************************************/

void 
Sankoff_probability_BW(Model B, Sank_S SK)
{
   register int    i, j, k, k1, n;
   int             Nmatrix;	/* number of matrices used  */
   double        **bwd_pptr, **bwr_pptr;	/* pointer for quick access */
   double        **bwd_pptr0, **bwr_pptr0;	/* pointer for quick access */
   double        **bwc_pptr, **bwc_pptr0;
   double         *bwc_ptr, *bwc_ptr0, *bwd_ptr, *bwd_ptr0, *bwr_ptr, *bwr_ptr0;
   double        **wd_pptr, **wr_pptr;	/* pointer for quick access */
   double        **wd_pptr0, **wr_pptr0;	/* pointer for quick access */
   double        **wc_pptr, **wc_pptr0;
   double         *wc_ptr, *wc_ptr0, *wd_ptr, *wd_ptr0, *wr_ptr, *wr_ptr0;
   double       ***pp_related;
   double          sum, sumk;
   double          tmp;
   short          *Orig0 = B->Seq->Orig[0];
   short          *Orig1 = B->Seq->Orig[1];
   int             work_blocks;	/* number of working matching blocks */
   int             Len0, Len1;
   int             maxdecrease;	/* maximum decrease factor          */
   float          *prior_N = SK->prior_N;
   double         *profile_ptr;
   double          P_M;

   Len0 = B->Seq->SeqLen[0];
   Len1 = B->Seq->SeqLen[1];
   pp_related = SK->rel_matrix.pp_related;
   Nmatrix = SK->rel_matrix.Nmatrix;

   /* allocate space for profile */
   NEWP(SK->profile, Len0 + 1, double);
   for (i = 0; i <= Len0; i++)
   {
      NEW(SK->profile[i], Len1 + 1, double);
   }
   /* initialize the profile */
   for (i = 0; i <= Len0; i++)
   {
      profile_ptr = SK->profile[i];
      for (j = 0; j <= Len1; j++)
      {
	 profile_ptr[j] = 0;
      }
   }

   maxdecrease = 0;
   for (n = 0; n < Nmatrix; n++)
   {
      if (maxdecrease < SK->factor[n][SK->Sankoff_blocks])
      {
	 maxdecrease = SK->factor[n][SK->Sankoff_blocks];
      }
   }

   /* calculate posterior probability */
   for (n = Nmatrix - 1; n >= 0; n--)
   {
      if (n != Nmatrix - 1)
      {				/* need to recalculate the WC, WR,WD */
	 /* As the memories located for WC,WR and WD is */
	 /* (SK->Sankoff_Blocks+1)*(Len0+1)*(Len1+1), so */
	 /* we do not need to re-initialize them here.  */
	 /* Also, we do not need to keep track decrease */
	 /* factors because they were already filled.   */

	 /* calculate the sum of probability */
	 for (k = 1; k <= SK->Sankoff_blocks; k++)
	 {
	    wd_pptr = SK->WD[k];
	    wd_pptr0 = SK->WD[k - 1];
	    wr_pptr = SK->WR[k];
	    wr_pptr0 = SK->WR[k - 1];
	    wc_pptr = SK->WC[k];
	    wc_pptr0 = SK->WC[k - 1];

	    /* keep the sum in range */
	    tmp = wd_pptr0[Len0][Len1] + wr_pptr0[Len0][Len1] +
	       wc_pptr0[Len0][Len1];
	    if (tmp > 1.0e30)
	    {			/* decrease by a factor */
	       for (i = 1; i < Len0; i++)
	       {
		  for (j = 1; j <= Len1; j++)
		  {
		     wd_pptr0[i][j] = wd_pptr0[i][j] * 1.0e-300;
		     wr_pptr0[i][j] = wr_pptr0[i][j] * 1.0e-300;
		     wc_pptr0[i][j] = wc_pptr0[i][j] * 1.0e-300;
		  }
	       }
	    }
	    for (i = 1; i <= Len0; i++)
	    {
	       for (j = 1; j <= Len1; j++)
	       {
		  wd_pptr[i][j] = wc_pptr[i - 1][j] + wd_pptr[i - 1][j];
		  wr_pptr[i][j] = wc_pptr[i][j - 1] + wd_pptr[i][j - 1] +
		     wr_pptr[i][j - 1];
		  wc_pptr[i][j] = (wd_pptr0[i - 1][j - 1] +
				   wr_pptr0[i - 1][j - 1] +
				   wc_pptr[i - 1][j - 1]) * prior_N[i] *
		     pp_related[n][Orig0[i - 1]][Orig1[j - 1]];
	       }
	    }
	 }
      }
      /* calculate the sum of probability backward */
      for (k = 1; k <= SK->Sankoff_blocks; k++)
      {
	 bwd_pptr = SK->BWD[k];
	 bwd_pptr0 = SK->BWD[k - 1];
	 bwr_pptr = SK->BWR[k];
	 bwr_pptr0 = SK->BWR[k - 1];
	 bwc_pptr = SK->BWC[k];
	 bwc_pptr0 = SK->BWC[k - 1];

	 /* keep the sum in range */
	 tmp = bwd_pptr0[1][1] + bwr_pptr0[1][1] +
	    bwc_pptr0[1][1];
	 if (tmp > 1.0e30)
	 {			/* decrease by a factor */
	    for (i = 1; i < Len0 + 1; i++)
	    {
	       for (j = 1; j <= Len1 + 1; j++)
	       {
		  bwd_pptr0[i][j] = bwd_pptr0[i][j] * 1.0e-300;
		  bwr_pptr0[i][j] = bwr_pptr0[i][j] * 1.0e-300;
		  bwc_pptr0[i][j] = bwc_pptr0[i][j] * 1.0e-300;
	       }
	    }
	 }
	 for (i = Len0; i >= 1; i--)
	 {
	    for (j = Len1; j >= 1; j--)
	    {
	       bwd_pptr[i][j] = bwc_pptr[i + 1][j] + bwd_pptr[i + 1][j] +
		  bwr_pptr[i + 1][j];
	       bwr_pptr[i][j] = bwc_pptr[i][j + 1] + bwr_pptr[i][j + 1];
	       bwc_pptr[i][j] = (bwd_pptr0[i + 1][j + 1] +
				 bwr_pptr0[i + 1][j + 1] +
				 bwc_pptr[i + 1][j + 1]) * prior_N[i] *
		  pp_related[n][Orig0[i - 1]][Orig1[j - 1]];
	    }
	 }
      }

      /*
       * printf("\n\n"); for(k=1;k<=SK->Sankoff_blocks;k++){
       * printf("wc[%d]=%g, wr[%d]=%g, wd[%d]=%g, w[%d]=%g\n",
       * k,SK->WC[k][Len0][Len1], k,SK->WR[k][Len0][Len1],
       * k,SK->WD[k][Len0][Len1],
       * k,SK->WD[k][Len0][Len1]+SK->WR[k][Len0][Len1]+SK->WC[k][Len0][Len1]);
       * printf("bwc[%d]=%g, bwr[%d]=%g, bwd[%d]=%g,bw[%d]=%g\n",
       * k,SK->BWC[k][1][1], k,SK->BWR[k][1][1], k,SK->BWD[k][1][1],
       * k,SK->BWD[k][1][1]+SK->BWR[k][1][1]+SK->BWC[k][1][1]); }
       */

      /* Need to keep every matrix in the same scale */
      /* and consider profile is floating point vaviable */
      /* for(k=SK->Sankoff_blocks; k ?????????k */

      /* add profile */
      for (k = SK->Sankoff_blocks; k >= 1; k--)
      {
	 for (i = 1; i <= Len0; i++)
	 {
	    wc_ptr = SK->WC[k][i];
	    profile_ptr = SK->profile[i];
	    for (j = 1; j <= Len1; j++)
	    {
	       for (k1 = 0; k1 <= SK->Sankoff_blocks - k; k1++)
	       {
		  tmp = wc_ptr[j] * (SK->BWC[k1 + 1][i + 1][j + 1] +
				     SK->BWR[k1][i][j + 1] +
				     SK->BWD[k1][i + 1][j]) / SK->N[k + k1];

		  profile_ptr[j] += wc_ptr[j] *
		     (SK->BWC[k1 + 1][i + 1][j + 1] +
		      SK->BWR[k1][i][j + 1] +
		      SK->BWD[k1][i + 1][j]) / SK->N[k + k1];
	       }
	    }
	 }
      }
   }

   /* normalize the profile */
   P_M = 0.5 / (SK->PM * Nmatrix * SK->Sankoff_blocks);
   for (i = 1; i <= Len0; i++)
   {
      profile_ptr = SK->profile[i];
      for (j = 1; j <= Len1; j++)
      {
	 profile_ptr[j] = profile_ptr[j] * P_M;
      }
   }

   /* find the Bayesian alignment */
   Sankoff_Bayesian_align(B, SK);

   /* free the profile */
   for (i = 0; i <= Len0; i++)
   {
      free(SK->profile[i]);
   }
   free(SK->profile);

}


/******************************** SC ************************/
/* */
/* calculate the expected score from statistics            */
/************************************************************/

double 
SC(Model B, Sank_S SK, int K)
{
   int             i, j;
   int             N = 0;	/* sum(unmatched X*unmatch Y)  */
   int             start_x, start_y;	/* start of unmatching blocks  */
   int             nx, ny;	/* length of unmatching blocks */
   double          Kp, score;

   start_x = 1;
   start_y = 1;
   for (i = K - 1; i >= 1; i--)
   {
      if (SK->start_x[K][i] != 0)
      {
	 nx = SK->start_x[K][i] - start_x;
	 ny = SK->start_y[K][i] - start_y;
	 N = N + nx * ny;
	 start_x = SK->end_x[K][i] + 1;
	 start_y = SK->end_y[K][i] + 1;
      }
   }
   /* the unmatching blocks at the tail */
   if (start_x <= B->Seq->SeqLen[0] && start_y <= B->Seq->SeqLen[1])
      N = N + (B->Seq->SeqLen[0] - start_x) * (B->Seq->SeqLen[1] - start_y);

   /* calculate the expected score with 0.1 probability */
   Kp = 7.227e7;
   score = log(Kp * N);
   return score;
}

/* find the maximum number */
float 
max3(float a, float b, float c)
{
   float           d;
   d = a;
   if (b > d)
      d = b;
   if (c > d)
      d = c;
   return d;
}
