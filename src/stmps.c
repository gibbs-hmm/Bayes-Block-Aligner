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


/*******************  Sankoff_Back_Sampling ************************/
/* */
/* Back Sampling the alignment                                    */
/* The results are saved in SK->start_x[0][k], and                */
/* SK->end_x[0][k]                                                */
/*******************************************************************/
#include "Sankoff.h"
void Sankoff_Back_Sampling(Model B, Sank_S SK)
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
   double          tmp, tmp1, w;
   double       ***WC, ***WD, ***WR;	/* for quick access */
   double        **wd_pptr, **wd_pptr0;	/* for quick access */
   double        **wc_pptr, **wc_pptr0;
   double        **wr_pptr, **wr_pptr0;
   double       ***pp_related;
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
   void            print_link_list(save_loop *);
   void            Sankoff_probability_W(Model B, Sank_S SK);
   void            Sankoff_prob_RC(Model B, Sank_S SK, int level, int nMat);
   void            prob_RC(Model B, Sank_S SK, int level, int nMat);
   void            print_level(Model B, Sank_S SK, int level, int reallevel);
   int             level_change, level, id;
   int             work_blocks;

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

   /* set seed */
   sRandom((long) time(NULL));

	 /* hard coded */
	 /* work_blocks=SK->Sankoff_blocks+1; */
	 if (SK->flags.min_memory)
	    work_blocks = 2;
	 else
	    work_blocks = SK->Sankoff_blocks + 1;


   /* sRandom(01010101L); */
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
	 if (SK->flags.min_memory)
	 {
	     Sankoff_prob_RC(B, SK, SK->Sankoff_blocks, N);
	   /*   Sankoff_probability_W(B,SK); */
	 } else
	 {
	    for (k = 1; k <= SK->Sankoff_blocks; k++)
	    {
	       print_level(B, SK, k, k);
	       wd_pptr = SK->WD[k];
	       wd_pptr0 = SK->WD[k - 1];
	       wr_pptr = SK->WR[k];
	       wr_pptr0 = SK->WR[k - 1];
	       wc_pptr = SK->WC[k];
	       wc_pptr0 = SK->WC[k - 1];

	       /* keep track of decreasing factor */
	       SK->factor[N][k] = SK->factor[N][k - 1];

	       tmp = wd_pptr0[Len0][Len1] + wr_pptr0[Len0][Len1] +
		  wc_pptr0[Len0][Len1];
	       if (tmp > 1.0e30)
	       {		/* decrease by a factor */
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
	 }
      }
      WC = SK->WC;
      WR = SK->WR;
      WD = SK->WD;

      fprintf(B->IP->Datafiles->out_fpt,
	      "%d samples are from Matrix %s.\n",
	      Nback_sampling[N], matrix_name[N]);
      SK->top = NULL;		/* set the current task to the top */
      for (i = 0; i < Nback_sampling[N]; i++)
      {
	 if (SK->backsampling_blocks > SK->Sankoff_blocks)
	 {
	    /* sampling with random number of blocks according MagL */
	    /* get a random number between 0-1 */
	    tmp = Random() / (double) LONG_MAX;
	    tmp1 = SK->PKSI[N][1];
	    k = 1;
	    while (k < SK->Sankoff_blocks && tmp1 < tmp)
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
           printf("Now processing %d backsamples\n",id);
          /* init level */
	  print_link_list(SK->top);

	 current = get_next(SK->top);
	 level = current->k;
	 Sankoff_prob_RC(B, SK, current->k, N);
	 print_level(B, SK, current->k % work_blocks, current->k);
	 while (current = get_next(SK->top))
	 {
	    /* printf("current level=%d\n",level); */
	    if (level != current->k)
	    {
	       printf("need to recalculate probability matrix\n");


	         Sankoff_prob_RC(B, SK, current->k, N);
		 /* Sankoff_probability_W(B,SK); */
	       print_level(B, SK, current->k % work_blocks, current->k);
	       WC = SK->WC;
	       WR = SK->WR;
	       WD = SK->WD;
	       level = current->k;
	    }
	    level_change = FALSE;
	    while (current->startb_x > 0 && current->startb_y > 0 && current->k > 0 && level_change != TRUE)
	    {
	       /* Get a random number between 0-1 */
	       tmp = Random() / (double) LONG_MAX;

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
                printf("completed backsample id=%d\n",current->id);

	       /* print the score */
	       /*
	        * fprintf(B->IP->Datafiles->out_fpt,"k= %d , score= %3f\n",
	        * Nblocks,score);
	        */
	       /* set profile */
	       for (k = 1; k <= current->Nblocks; k++)
	       {
		  if (current->start_x[0][current->k] != 0)
		  {		/* find a block */
		     for (n = current->start_x[0][k], j = current->start_y[0][k];
			  n <= current->end_x[0][k]; n++, j++)
		     {
			SK->profile[n][j]++;
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
				 (char) (SEQ2CH(Seq->Orig[0][n - 1]) + 65));
			}
			fprintf(B->IP->Datafiles->out_fpt, "% 3d\n",
				current->end_x[0][k]);
			fprintf(B->IP->Datafiles->out_fpt, "%3d ",
				current->start_y[0][k]);
			for (n = current->start_y[0][k]; n <= current->end_y[0][k]; n++)
			{
			   fprintf(B->IP->Datafiles->out_fpt, "%c",
				 (char) (SEQ2CH(Seq->Orig[1][n - 1]) + 65));
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
	    } else
	    {
	       /* printf("level change to %d \n",current->k); */
	    }
	 }
      }
   printf("backsampling number is %d\n",SK->backsampling_number);
   /* normalize the profile */
   for (i = 1; i <= Len0; i++)
   {
      profile_ptr = SK->profile[i];
      for (j = 1; j <= Len1; j++)
      {
	 profile_ptr[j] = profile_ptr[j] / SK->backsampling_number;
       if (profile_ptr[j] > .000001)
            printf("i=%d,j=%d,p=%f\n",i,j, profile_ptr[j]);
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

/******************get_next*************************/
/* */
/* get the next backsample at the maximum level    */
/***************************************************/
save_loop      *
get_next(save_loop * top)
{
   save_loop      *tmp_max;
   tmp_max = top;
   while (top != NULL)
   {
      if (tmp_max->k <= top->k && top->startb_x > 0 && top->startb_y > 0 &&
	  (top->done != TRUE))
	 tmp_max = top;
      top = top->next;
   }
   if (tmp_max->k == 0 || tmp_max->done == TRUE)
   {
      return (NULL);
   } else
   {
      return (tmp_max);
   }
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
/* */
/* Allocate the basic memory need by Sankoff   */
/* procedure.                                  */
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
