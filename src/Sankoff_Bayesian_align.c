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
/*                                                                 */
/* $Id: Sankoff_Bayesian_align.c,v 1.6 1997/11/21 22:08:54 junzhu Exp junzhu $   */
/*                                                                 */
/* Author: Jun Zhu                                                 */
/* Date:   July 1, 1996                                            */
/*                                                                 */
/*******************************************************************/

#include "string.h"
#include "common.h"
#include "Sankoff.h"

/********** Sankoff_Bayesian_align(Model B,Sank_S SK) **************/
/*                                                                 */
/*  Find the alignment according to joint probability              */
/*******************************************************************/

void Sankoff_Bayesian_align(Model B, Sank_S SK)
{

  int i,j,k,w;
  int i1,i2,j1,j2,pos;
  int *rank,*max_i,*max_j;
  double *max,tmp;
  int exist;
  double sum,sparse_cutoff;
  double  *marginal_prob; /*marginal probability of the query aligns with data sequence */
  int Nblocks,Pblocks;
  double **tmp_profile;                 /* tmp copy of profile */
  char xlabel[20],ylabel[20];
  char filename[80],extension[10];
  char filename_tmp0[80],filename_tmp1[80],filename_tmp2[80];
  char tmpch;
  FILE *align_fptr, *profile_fptr, *marginal_fptr, *bayes_info_fptr;
  int X_start,X_end,Y_start,Y_end;
  int *x_start,*x_end,*y_start,*y_end;
  int strl;                             /* temporary string length */         

   NEW(marginal_prob,B->Seq->SeqLen[0] + 1, double); /* allocate space for marginal probability */

  /* copy SK->profile to tmp_profile */
  NEWP(tmp_profile,B->Seq->SeqLen[0]+1 ,double);
  for(i=0;i<=B->Seq->SeqLen[0];i++){
      NEW(tmp_profile[i],B->Seq->SeqLen[1]+1,double);
      memcpy(tmp_profile[i],SK->profile[i],
	     sizeof(double)*(B->Seq->SeqLen[1]+1));
     marginal_prob[i]=0;  /* initialize */
  }

  Nblocks=SK->Sankoff_blocks;
  NEW(rank,Nblocks+1,int);
  NEW(max,Nblocks+1,double);
  NEW(max_i,Nblocks+1,int);
  NEW(max_j,Nblocks+1,int);

  /* allocate space for ends of aligned segments */
  NEW(x_start,Nblocks+1,int);
  NEW(x_end,Nblocks+1,int);
  NEW(y_start,Nblocks+1,int);
  NEW(y_end,Nblocks+1,int);

  /* construct an output file a file  use first 15 characters */
  j=0;
  strl=strlen(SK->comment[0]);
  if(strl>15) strl=15;

  for(i=0;i<strl;i++){
      if(SK->comment[0][i]!=' '&& SK->comment[0][i]!='>' &&
	 SK->comment[0][i]!='|'){
	  filename[j]=SK->comment[0][i];
          xlabel[j]=SK->comment[0][i];  /* get a label for the profile plot */
	  j++;
      }
  }
  filename[j] ='-'; 
  xlabel[j]='\0';
  j++;
   
  strl=strlen(SK->comment[1]);
  if(strl>15) strl=15;
  k=0;
  for(i=0;i<strl;i++) {
      if(SK->comment[1][i]!=' '&& SK->comment[1][i]!='>' &&
	 SK->comment[1][i]!='|')
      { 
	  filename[j]=SK->comment[1][i];
          ylabel[k]=SK->comment[1][i]; /* get y axis label for profile plot */
          k++;
	  j++;
      }
  }
  ylabel[k]='\0' ;
  filename[j] ='\0'; 
  i=0;
  exist=1;
  while(exist==1){
      strcpy(filename_tmp0,filename);
      strcpy(filename_tmp1,filename);
      strcpy(filename_tmp2,filename);
      strcat(filename_tmp0,".profile");
      strcat(filename_tmp1,".alignment");
      strcat(filename_tmp2,".marginal");
      sprintf(extension,".%d",i);
      strcat(filename_tmp0,extension);  
      strcat(filename_tmp1,extension);
      strcat(filename_tmp2,extension);  
      align_fptr =NULL;
      profile_fptr=NULL;
      marginal_fptr=NULL;
      /* detect whether the same file exist before */
      if((align_fptr   = fopen(filename_tmp0,"r"))==NULL && 
	 (profile_fptr = fopen(filename_tmp1,"r"))==NULL){
	  exist=0;
      }
      else{/* file exist, increase extension */
	  if(align_fptr!=NULL)   fclose(align_fptr);
	  if(profile_fptr!=NULL) fclose(profile_fptr);
      }
      i++;
  }
  /* write names of profile,marginal,alignment files to bayes.info file */
  if ((bayes_info_fptr = fopen("bayes.info","w")) != NULL) {
      fprintf(bayes_info_fptr,"%d\n",B->Seq->SeqLen[0]);
      fprintf(bayes_info_fptr,"%d\n",B->Seq->SeqLen[1]);
      fprintf(bayes_info_fptr,"%s\n", filename_tmp1);
      fprintf(bayes_info_fptr,"%s\n", filename_tmp0);
      fprintf(bayes_info_fptr,"%s\n", xlabel);
      fprintf(bayes_info_fptr,"%s\n", ylabel);
      fprintf(bayes_info_fptr,"%s\n", filename_tmp2);
      fclose(bayes_info_fptr);
  }
  else {
       fprintf(stderr,"unable to write to bayes.info\n");
  }

  /* profile */
  profile_fptr=fopen(filename_tmp0,"w");
  print_parameters_sankoff(profile_fptr,SK);
  fprintf(profile_fptr,"#%s\n",SK->comment[0]);
  fprintf(profile_fptr,"#%s\n",SK->comment[1]);
  if (SK->flags.profile == FALSE)
    fprintf(profile_fptr,"***********no profile data written because of NO_PROFILE OPTION *********\n");
   if (SK->flags.full_profile)
      sparse_cutoff=-1.0; /* -1 assures us that we will print out all values */
   else
      sparse_cutoff= SPARSE_CUTOFF_LIMIT;
      
  for(i=1;i<=B->Seq->SeqLen[0];i++){
      for(j=1;j<=B->Seq->SeqLen[1];j++){
	if ((SK->flags.profile == TRUE ) && ( SK->profile[i][j] > sparse_cutoff) || 
             (i==1 && j==1) || (i==B->Seq->SeqLen[0] && j==B->Seq->SeqLen[1]))
                    fprintf(profile_fptr,"%d %d %g\n",i,j,SK->profile[i][j]);
      
          marginal_prob[i]+=SK->profile[i][j];   /*add up all of the times i aligns with some j*/ 
      }
  }
        if (SK->flags.profile == TRUE)
            fprintf(profile_fptr,"\n\n");
  fclose(profile_fptr);

  /* print out the marginal probability that the query sequence aligns with some fashion with the data string */

  if (!(marginal_fptr=fopen(filename_tmp2,"w"))) {
      printf("unable to open %s\n",filename_tmp2);
      }
  else {     
       fprintf(marginal_fptr,"Title=\"%s\"\n",filename);
       fprintf(marginal_fptr,"Xaxis=\"Query string index\"\n");
       fprintf(marginal_fptr,"Yaxis=\"Probability of Alignment\"\n");
       fprintf(marginal_fptr,"Line=\"Probability\"\n");
       fprintf(marginal_fptr,"#%s\n",SK->comment[0]);
       fprintf(marginal_fptr,"#%s\n",SK->comment[1]);
       for(i=1;i<=B->Seq->SeqLen[0];i++){
         fprintf(marginal_fptr,"%d %g\n",i,marginal_prob[i]);
       }
       fclose(marginal_fptr);
 }
  /* do it again under another name, kluggy but it keeps me from having to pass the name back to the GUI interface */
 if (!(marginal_fptr=fopen("marginal.txt","w"))) {
      printf("unable to open marginal.txt \n");
      }
  else {
       fprintf(marginal_fptr,"Title=\"%s\"\n",filename);
       fprintf(marginal_fptr,"Xaxis=\"Query string index\"\n");
       fprintf(marginal_fptr,"Yaxis=\"Prob. of Alignment\"\n");
       fprintf(marginal_fptr,"Line=\"Probability\"\n");
       for(i=1;i<=B->Seq->SeqLen[0];i++){
         fprintf(marginal_fptr,"%d %g\n",i,marginal_prob[i]);
       }
       fclose(marginal_fptr);
 }

 free(marginal_prob);

  /* alignment */
  align_fptr=fopen(filename_tmp1,"w");
  print_parameters_sankoff(align_fptr,SK);
  fprintf(align_fptr,"Query %s\n",SK->comment[0]);
  fprintf(align_fptr,"Data %s\n\n",SK->comment[1]);

  /* find out the alignment */
  for(k=0;k<=Nblocks;k++){
      rank[k]=k;
      max[k]=0;
      /* find out the highest peak at the moment */
      for(i=1;i<=B->Seq->SeqLen[0];i++){
	  for(j=1;j<=B->Seq->SeqLen[1];j++){
	      if(SK->profile[i][j]>max[k]){
		  max[k]=SK->profile[i][j];
		  max_i[k]=i;
		  max_j[k]=j;
	      }
	  }
      }
      if(max[k]==0){
	  Nblocks=k;
	  break;
      }
     /*fprintf(align_fptr,"k=%d, max=%d, i=%d, j=%d\n",k,max[k],max_i,max_j);*/
      /* mark the alignment */
      /* go down stream */
      i1=max_i[k];
      j1=max_j[k];
      while(i1>=1 && j1>=1 && SK->profile[i1][j1]>=0.25*max[k] ){
	  /*  0.67 is set to detect two overlapped peaks */
	  SK->profile[i1][j1]=0;
	  i1--;
	  j1--;
      }
      /*fprintf(align_fptr," x1=%d,  y1=%d\n",i1+1,j1+1);*/
      x_start[k]=i1+1;
      y_start[k]=j1+1;
      /* go up stream */
      i2=max_i[k]+1;
      j2=max_j[k]+1;
      while(i2<=B->Seq->SeqLen[0] && j2<=B->Seq->SeqLen[1] 
	    && SK->profile[i2][j2]>=0.25*max[k] ){
	  SK->profile[i2][j2]=0;
	  i2++;
	  j2++;
      }
      /*fprintf(align_fptr," x2=%d,  y2=%d\n\n",i2-1,j2-1);*/
      x_end[k]=i2-1;
      y_end[k]=j2-1;
      /* mark (i1,j1) -(i2,j2) to zero */
      for(i=0;i<=x_end[k];i++)
	  for(j=y_start[k];j<=B->Seq->SeqLen[1];j++)
	      SK->profile[i][j]=0;
      for(i=x_start[k];i<=B->Seq->SeqLen[0];i++)
	  for(j=0;j<=y_end[k];j++)
	      SK->profile[i][j]=0;
      /* if the current block is the extension of an existing block */
      /* cencatenate them together  */
      for(i=0; i<k; i++){
	  if(x_end[i]+1==x_start[k] && y_end[i]+1==y_start[k]){
	      x_end[i]=x_end[k];
	      y_end[i]=y_end[k];
	      k--;
	  }
	  else if(x_start[i]-1==x_end[k] && y_start[i]-1==y_end[k]){
	      x_start[i]=x_start[k];
	      y_start[i]=y_start[k];
	      k--;
	  }
      }
  } 

  /* find out how many probable aligned segments */
  Pblocks=0;
  sum=SK->PK[0];
  while(1.0-sum >0.8){
      Pblocks++;
      sum+=SK->PK[Pblocks];
  }

  /* put the first Pblocks in order */
  for(k=0;k<Pblocks-1;k++){
      for(i=k+1;i<Pblocks;i++){
          if(x_start[i]<x_start[k]){
	      pos=x_start[i];
	      x_start[i]=x_start[k];
	      x_start[k]=pos;

	      pos=y_start[i];
	      y_start[i]=y_start[k];
	      y_start[k]=pos;

	      pos=x_end[i];
	      x_end[i]=x_end[k];
	      x_end[k]=pos;

	      pos=y_end[i];
	      y_end[i]=y_end[k];
	      y_end[k]=pos;

	      tmp=max[i];
	      max[i]=max[k];
	      max[k]=tmp;	      

	      pos=max_i[i];
	      max_i[i]=max_i[k];
	      max_i[k]=pos;

	      pos=max_j[i];
	      max_j[i]=max_j[k];
	      max_j[k]=pos;

	      tmp=rank[i];
	      rank[i]=rank[k];
	      rank[k]=tmp;
	  }
      }
  }

  /* print the alignment with weight graph  */
  /* assume the terminal is 80 columns wide */
  fprintf(align_fptr,
"\n-----------------Probable aligned segments-----------------\n\n");
  for(k=0;k<Pblocks;k++){
      fprintf(align_fptr,"Rank=%d, peak=%6.3f at (%d,%d)\n",rank[k]+1,
	      max[k],max_i[k],max_j[k]);    

      X_start=x_start[k];
      Y_start=y_start[k];
      X_end=(X_start+80<x_end[k])?X_start+80:x_end[k];
      Y_end=(Y_start+80<y_end[k])?Y_start+80:y_end[k];
      while(X_start <x_end[k]) {
	  fprintf(align_fptr,"%4d ",X_start);
	  for(i=X_start;i<=X_end;i++) {
	      if(SK->flags.protein_sequence) { /* protein */
	          fprintf(align_fptr,"%c",
		    (char )(SEQ2CH(SK->flags.protein_sequence,B->Seq->Orig[0][i-1])+65));
	      }
	      else  { /* DNA */
		  switch (B->Seq->Orig[0][i-1]) {
		      case 0: tmpch='A'; break;
		      case 1: tmpch='T'; break;
		      case 2: tmpch='C'; break;
		      case 3: tmpch='G'; break;
		  }
		  fprintf(align_fptr,"%c",tmpch);
	      }
	  }	       
	  fprintf(align_fptr," %4d\n",X_end);
	  /* print weights */
	  for(w=(int)(max[k]*10/2)*2;w>=0;w-=2){
	      fprintf(align_fptr,"     ");
	      for(i=0;i<=X_end-X_start;i++){
		  if(tmp_profile[i+X_start][i+Y_start]>=
		     (w+1.5)/10)
		  { /* print : */
		      fprintf(align_fptr,":");
		  }
		  else if(tmp_profile[i+X_start][i+Y_start]>=
		      (w+0.5)/10)
		  { /* print . */
		      fprintf(align_fptr,".");
		  }
		  else /* print \s */
		      fprintf(align_fptr," ");		
	      }
	      fprintf(align_fptr,"\n"); 
	  }
	  fprintf(align_fptr,"%4d ",Y_start);
	  for(j=Y_start;j<=Y_end;j++) {
	      if(SK->flags.protein_sequence) { /* protein */
	          fprintf(align_fptr,"%c",
		    (char )(SEQ2CH(SK->flags.protein_sequence,B->Seq->Orig[1][j-1])+65));
	      }
	      else  { /* DNA */
		  switch (B->Seq->Orig[1][j-1]) {
		      case 0: tmpch='A'; break;
		      case 1: tmpch='T'; break;
		      case 2: tmpch='C'; break;
		      case 3: tmpch='G'; break;
		  }
		  fprintf(align_fptr,"%c",tmpch);
	      }
	  }
	  fprintf(align_fptr," %4d\n\n",Y_end);

	  /* update X,Y start and end */
	  X_start=X_end+1;
	  Y_start=Y_end+1;
	  X_end=(X_start+80<x_end[k])?X_start+80:x_end[k];
	  Y_end=(Y_start+80<y_end[k])?Y_start+80:y_end[k];
      }
  }

  fprintf(align_fptr,
	  "\n------------The rest of (%d-%d) seqments----------\n\n",
	  Nblocks,Pblocks);
  for(k=Pblocks;k<Nblocks;k++){
      fprintf(align_fptr,"Rank=%3d, peak=%6.3f at (%d,%d)\n",
	      k+1,max[k],max_i[k],max_j[k]);
      fprintf(align_fptr,"%3d ",x_start[k]);
      for(i=x_start[k];i<=x_end[k];i++) {
	  if(SK->flags.protein_sequence) { /* protein */
	      fprintf(align_fptr,"%c",
		      (char )(SEQ2CH(SK->flags.protein_sequence,B->Seq->Orig[0][i-1])+65));
	  }
	  else  { /* DNA */
	      switch (B->Seq->Orig[0][i-1]) {
	          case 0: tmpch='A'; break;
	          case 1: tmpch='T'; break;
		  case 2: tmpch='C'; break;
		  case 3: tmpch='G'; break;
	      }
	      fprintf(align_fptr,"%c",tmpch);
	  }
      }
      fprintf(align_fptr," %3d\n",x_end[k]);
      fprintf(align_fptr,"%3d ",y_start[k]);
      for(j=y_start[k];j<=y_end[k];j++) {
	  if(SK->flags.protein_sequence) { /* protein */
	      fprintf(align_fptr,"%c",
		      (char )(SEQ2CH(SK->flags.protein_sequence,B->Seq->Orig[1][j-1])+65));
	  }
	  else  { /* DNA */
	      switch (B->Seq->Orig[1][j-1]) {
	          case 0: tmpch='A'; break;
	          case 1: tmpch='T'; break;
	          case 2: tmpch='C'; break;
	          case 3: tmpch='G'; break;
	      }
	      fprintf(align_fptr,"%c",tmpch);
	  }
      }
      fprintf(align_fptr," %3d\n\n",y_end[k]);
  }
  fclose(align_fptr);

  free(rank);
  free(max);
  free(max_i);
  free(max_j);

  /* free space for ends of aligned segments */
  free(x_start);
  free(x_end);
  free(y_start);
  free(y_end);

}










