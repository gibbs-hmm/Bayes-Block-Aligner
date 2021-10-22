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

/***************************************************************************/
/* $Id: mem_mgmt.c,v 1.9 1997/03/21 19:30:45 junzhu Exp junzhu $           */
/*                                                                         */
/* Author:         Eric C. Rouchka, 1996.6                                 */
/*                 Jun Zhu, 1996.8                                         */
/*                                                                         */
/***************************************************************************/

#include <math.h>
#include "common.h"

/******************* alloc_Model ***************************************/
/*                                                                     */
/* DESCRIPTION : Allocates the memory for the basic Model and returns  */
/*               a structure                                           */
/***********************************************************************/

Model alloc_Model(void)

{
   Model B;

   NEW(B, 1, Modelstruct);
   NEW(B->IP, 1, InputParams);
   NEW(B->C, 1, Counts);
   NEW(B->First, 1, Counts);
   NEW(B->Seq, 1, Stringstruct);
   NEWP(B->Seq->R, 1, char);
   NEWP(B->Seq->nvEndLocs, 1, int);
   NEWP(B->Seq->ProcessedSTR, 1, char);
   return B;
}
 
/********************  alloc_Counts   *******************************/

void alloc_Counts(Model B)

{
   int    i,j,t,maxlen;
   IPtype IP;           /* for quick access data */
   Ctype  First, C;     /* for quick access data */
 
   IP=B->IP;
   First=B->First;
   C=B->C;

   /* maximum motif length */
   maxlen = findMaxMotifLen(IP);           

   /* C->dSumPseudo[MotifType][BG/MOTIF][nAlphaLen] */
   NEWPP(C->dSumPseudo,  IP->nNumMotifTypes, double);               
   NEWPP(First->dSumPseudo, IP->nNumMotifTypes, double);

   /* C->dvPseudoCounts[MotifType][PositionInMotif][nAlphaLen]  */
   NEWPP(C->dPseudoCounts,IP->nNumMotifTypes , double);
   NEWPP(First->dPseudoCounts,IP->nNumMotifTypes, double);

   for(t=0;t< IP->nNumMotifTypes;t++){
     NEWP(First->dSumPseudo[t],2,double);
     NEWP(C->dSumPseudo[t],2,double);

     NEW(C->dSumPseudo[t][BG], IP->nAlphaLen, double);  
     NEW(C->dSumPseudo[t][MOTIF], IP->nAlphaLen, double);
     NEW(First->dSumPseudo[t][BG], IP->nAlphaLen, double); 
     NEW(First->dSumPseudo[t][MOTIF], IP->nAlphaLen, double);
     
     NEWP(C->dPseudoCounts[t],maxlen+1, double);
     NEWP(First->dPseudoCounts[t],maxlen+1, double);    
     for(j = 0; j <= maxlen; j++) {
         NEW(C->dPseudoCounts[t][j], IP->nAlphaLen, double);
         NEW(First->dPseudoCounts[t][j],IP->nAlphaLen, double);
      }
   }
   First->nTotBack = 0;
   NEWP(First->dmodel_sites, IP->nNumMotifTypes, double);
   NEWP(C->dmodel_sites, IP->nNumMotifTypes, double);
   NEWP(First->dmodel_pseudo, IP->nNumMotifTypes, double);
   NEWP(C->dmodel_pseudo, IP->nNumMotifTypes, double);
   NEW(First->dtot_sites, IP->nNumMotifTypes, double);
   NEW(C->dtot_sites, IP->nNumMotifTypes, double);
   NEW(First->dtot_pseudo, IP->nNumMotifTypes, double);
   NEW(C->dtot_pseudo, IP->nNumMotifTypes, double);
   NEW(First->dTot, IP->nNumMotifTypes, double);
   NEW(C->dTot, IP->nNumMotifTypes, double);
   NEW(First->dbg_pseudo, IP->nNumMotifTypes, double);
   NEW(C->dbg_pseudo, IP->nNumMotifTypes, double);

   /* First->fCounts[MotifType][PositionInMotif][nAlphabet] */
   NEWPP(First->fCounts,IP->nNumMotifTypes, float);
   NEWPP(C->fCounts, IP->nNumMotifTypes, float);
   for(t=0;t<IP->nNumMotifTypes;t++){
     NEWP(First->fCounts[t],maxlen+1, float);
     NEWP(C->fCounts[t], maxlen+1, float);
     for(i = 0; i <= maxlen; i++) { 
      /* we also allocate space for symbol represent any A.A's 
         or nucleotides                                          */
      NEW(First->fCounts[t][i], IP->nAlphaLen+1, float);
      NEW(C->fCounts[t][i], IP->nAlphaLen+1, float);
     }
   }
   for(t = 0; t < IP->nNumMotifTypes; t++) {
       NEW(First->dmodel_sites[t], IP->RevComplement + 1, double);
       NEW(C->dmodel_sites[t], IP->RevComplement + 1, double);
       NEW(First->dmodel_pseudo[t], IP->RevComplement + 1, double);
       NEW(C->dmodel_pseudo[t], IP->RevComplement + 1, double);
   }
}


/************ alloc_min_Sank************************/
/*                                             */
/* Allocate the basic memory need by Sankoff   */
/* procedure.                                  */
/***********************************************/




/************ alloc_Sank************************/
/*                                             */
/* Allocate the basic memory need by Sankoff   */
/* procedure.                                  */
/***********************************************/

void alloc_Sank(Model B, Sank_S SK)
{
     int        i,j,k,n;
     float      *s_ptr,*sum_ptr;
     short      *x_ptr;
     double     *p_ptr;

     /* allocate space for terminals of matching blocks */
     NEWP(SK->start_x,SK->Sankoff_blocks+1,short);
     NEWP(SK->start_y,SK->Sankoff_blocks+1,short);
     NEWP(SK->end_x,SK->Sankoff_blocks+1,short);
     NEWP(SK->end_y,SK->Sankoff_blocks+1,short);
     for(k=0;k<=SK->Sankoff_blocks;k++){
         NEW(SK->start_x[k],SK->Sankoff_blocks+1,short);
	 NEW(SK->start_y[k],SK->Sankoff_blocks+1,short);
	 NEW(SK->end_x[k],SK->Sankoff_blocks+1,short);
	 NEW(SK->end_y[k],SK->Sankoff_blocks+1,short);
     }
     /* initialize the terminal matching blocks */
     for(i=0;i<=SK->Sankoff_blocks;i++) {
         for(j=0;j<=SK->Sankoff_blocks;j++){
	     SK->start_x[i][j]=0;
	     SK->start_y[i][j]=0;
	     SK->end_x[i][j]=0;
	     SK->end_y[i][j]=0;
	 }
     }
  
     /* allocate space for PK, PSI,PKSI, W and N */
     NEW(SK->PK,SK->Sankoff_blocks+1,double); 
     NEW(SK->PSI,SK->rel_matrix.Nmatrix,double);       
     NEWP(SK->PKSI,SK->rel_matrix.Nmatrix,double); 
     NEWP(SK->W,SK->rel_matrix.Nmatrix,double);
     NEWP(SK->factor,SK->rel_matrix.Nmatrix,int);       
     NEW(SK->N,SK->Sankoff_blocks+1,double);
     for (i=0;i<=SK->rel_matrix.Nmatrix-1;i++){
         NEW(SK->PKSI[i],SK->Sankoff_blocks+1,double);	 
	 NEW(SK->W[i],SK->Sankoff_blocks+1,double);
	 NEW(SK->factor[i],SK->Sankoff_blocks+1,int);
	 for(j=0;j<SK->Sankoff_blocks+1;j++){
	     SK->factor[i][j]=0;
	 }
     }

     NEW(SK->prior_N,B->Seq->SeqLen[0]+1,float);
     for (i=1;i<=B->Seq->SeqLen[0];i++){
	SK->prior_N[i]=1.0;
     }     
    
     for(i=0;i<NInput_N;i++){
         SK->prior_N[Input_N[i]]+=Input_alpha;
     }

}
/************ alloc_Sank_S ************************/
/*                                                */
/*  Allocate space for S                          */
/**************************************************/

void alloc_Sank_S(Model B, Sank_S SK)
{
     register int i,j,k,n;
     int   Nmatrix;
     int   L0,L1;  /* sequence length */
     float **sv_pptr,**sw_pptr;

     Nmatrix=SK->rel_matrix.Nmatrix;
     L0=B->Seq->SeqLen[0];
     L1=B->Seq->SeqLen[1];

     /* allocate score space  */
     NEWP3(SK->SV,Nmatrix,float); 
     NEWP3(SK->SW,Nmatrix,float);    
     for(n=0; n<Nmatrix; n++) {
         NEWPP(SK->SV[n],SK->Sankoff_blocks+1,float);
         NEWPP(SK->SW[n],SK->Sankoff_blocks+1,float);
	 for(k=0; k<=SK->Sankoff_blocks; k++){
	     NEWP(SK->SV[n][k],L0+1,float);
	     NEWP(SK->SW[n][k],L0+1,float);
	     for(i=0; i<=L0; i++){
	         NEW(SK->SV[n][k][i],L1+1,float);
	         NEW(SK->SW[n][k][i],L1+1,float);
	     }
	 }
     }

     /* initialize the matrix SV and SW  */
     /* initialize SV[n][0][i][j]=0,SW[n][0][i][j]=0 */
     for(n=0; n<Nmatrix; n++) {
         sv_pptr=SK->SV[n][0];
         sw_pptr=SK->SW[n][0];
         for(i=0; i<=L0; i++) {
	     for(j=0; j<=L1; j++) {
	         sv_pptr[i][j]=0;
		 sw_pptr[i][j]=0;
	     }
	 }
     }

     /* initialize SV[n][k][0][j]=0,SW[n][k][0][j]=0  */
     /* initialize SV[n][k][i][0]=0,SW[n][k][i][0]=0  */
     for(n=0; n<Nmatrix; n++) {
         for(k=0; k<=SK->Sankoff_blocks; k++) {
	     sv_pptr=SK->SV[n][k];
	     sw_pptr=SK->SW[n][k];
	     for(i=0; i<=L0; i++) {
	         sv_pptr[i][0]=0;
		 sw_pptr[i][0]=0;
	     }
	     for(j=0; j<=L1; j++) {
	         sv_pptr[0][j]=0; 
		 sw_pptr[0][j]=0;
	     }
	 }
     }
}
     
/************ alloc_Sank_N  ************************/
/*                                                 */
/*  Allocate space for N                           */
/***************************************************/
     
void alloc_Sank_N(Model B, Sank_S SK)
{
     register int i,j,k,n;
     int      Len0,Len1;               /* for quick access */
     double   *nd_ptr,*nr_ptr,*nc_ptr; /* for quick access */

     Len0=B->Seq->SeqLen[0];
     Len1=B->Seq->SeqLen[1];

     /* allocate space for ND, NR and NC */
     NEWPP(SK->ND,2,double);
     NEWPP(SK->NR,2,double);
     NEWPP(SK->NC,2,double);
     for(i=0;i<=1;i++){
         NEWP(SK->ND[i],Len0+1,double);
         NEWP(SK->NR[i],Len0+1,double);
         NEWP(SK->NC[i],Len0+1,double);
         for(j=0;j<=Len0;j++){
	     NEW(SK->ND[i][j],Len1+1,double);
	     NEW(SK->NR[i][j],Len1+1,double);
	     NEW(SK->NC[i][j],Len1+1,double);
	 }
     }

     /* initialize the NC[0][i][j]=0              */
     /* ND[0][i][j]=0, NR[0][i][j]=1  for i,j !=0 */
     for(i=1;i<=Len0;i++){
	 nc_ptr=SK->NC[0][i];
	 nd_ptr=SK->ND[0][i];
	 nr_ptr=SK->NR[0][i];

	 for(j=1;j<=Len1;j++){
	     nc_ptr[j]=0.0;
	     nd_ptr[j]=0.0;
	     nr_ptr[j]=1.0;
	 }
     }

     /* initialize the NC[0][i][0]=0, NC[0][0][j]=0 */
     /* initialize the ND[0][i][0]=1, ND[0][0][j]=0 */
     /* initialize the NR[0][i][0]=0, NR[0][0][j]=1 */
     /* for i,j !=0                                 */

     for(i=1;i<=Len0;i++){
         SK->NC[0][i][0]=0;
	 SK->ND[0][i][0]=1;
	 SK->NR[0][i][0]=0;
     }
     
     nc_ptr=SK->NC[0][0];
     nd_ptr=SK->ND[0][0];
     nr_ptr=SK->NR[0][0];

     for(j=1;j<=Len1;j++){
         nc_ptr[j]=0;
	 nd_ptr[j]=0;
	 nr_ptr[j]=1;
     }

     /* initialize the NC[0][0][0]=0 */
     /* initialize the ND[0][0][0]=0 */
     /* initialize the NR[0][0][0]=1 */
     SK->NC[0][0][0]=0;
     SK->ND[0][0][0]=0;
     SK->NR[0][0][0]=1;

     /* initialize the NC[k][i][j]=0 */
     /* initialize the ND[k][i][j]=0 */
     /* initialize the NR[k][i][j]=0 */
     /* for k!=0, and i=0 or j=0     */
     for(i=0; i<=Len0; i++) {
         SK->NC[1][i][0]=0;
	 SK->ND[1][i][0]=0;
	 SK->NR[1][i][0]=0;
     }   
     nc_ptr=SK->NC[1][0];
     nd_ptr=SK->ND[1][0];
     nr_ptr=SK->NR[1][0];        
     for(j=0; j<=Len1; j++) {
         nc_ptr[j]=0;
	 nd_ptr[j]=0;
	 nc_ptr[j]=0;
     }
}

/************ alloc_Sank_W  ************************/
/*                                                 */
/*  Allocate space for W                           */
/***************************************************/
     
void alloc_Sank_W(Model B, Sank_S SK)
{
     register int i,j,k,n;
     int      Len0,Len1;               /* for quick access */
     double   *wd_ptr,*wr_ptr,*wc_ptr; /* for quick access */
     int      work_blocks;            
   
     Len0=B->Seq->SeqLen[0];
     Len1=B->Seq->SeqLen[1];


     if(SK->flags.back_sampling || SK->flags.exact_posterior){
         /* allocate whole space */
         work_blocks=SK->Sankoff_blocks+1;
     }
     else{/* only two array for working space  */
         work_blocks=2;
     }
     /* minimum memory mode even for backsampling */
     if(SK->flags.min_memory && !SK->flags.exact_posterior)
      work_blocks=2; 

     /* allocate space for WD,WR and WC */
     NEWPP(SK->WD,work_blocks,double);
     NEWPP(SK->WR,work_blocks,double);
     NEWPP(SK->WC,work_blocks,double);

     for(i=0;i<work_blocks;i++){
         NEWP(SK->WD[i],Len0+1,double);
	 NEWP(SK->WR[i],Len0+1,double);
	 NEWP(SK->WC[i],Len0+1,double);
	 for(j=0;j<=Len0;j++){
	     NEW(SK->WD[i][j],Len1+1,double);
	     NEW(SK->WR[i][j],Len1+1,double);
	     NEW(SK->WC[i][j],Len1+1,double);
	 }
     }

     /* initialize the WC[0][i][j]=0  */
     /* WD[0][i][j]=0, WR[0][i][j]=1  */
     /* for i,j !=0                   */

     for(i=1;i<=Len0;i++){
         wc_ptr=SK->WC[0][i];
	 wd_ptr=SK->WD[0][i];
	 wr_ptr=SK->WR[0][i];

	 for(j=1;j<=Len1;j++){
	     wc_ptr[j]=0.0;
	     wd_ptr[j]=0.0;
	     wr_ptr[j]=1.0e-300;	
	 }
     }

     /* initialize WC[0][i][0]=0 */
     /* initialize WD[0][i][0]=1 */
     /* initialize WR[0][i][0]=0 */
     /* for i!=0                     */
     for(i=1;i<=Len0;i++){
         SK->WC[0][i][0]=0.0;
	 SK->WD[0][i][0]=1.0e-300;
	 SK->WR[0][i][0]=0.0;
     }

     /* initialize WC[0][0][j]=0 */
     /* initialize WD[0][0][j]=0 */
     /* initialize WR[0][0][j]=1 */
     /* for any j                */
     wc_ptr=SK->WC[0][0];
     wd_ptr=SK->WD[0][0];
     wr_ptr=SK->WR[0][0];
     for(j=0;j<=Len1;j++){
         wc_ptr[j]=0.0;
	 wd_ptr[j]=0.0;
	 wr_ptr[j]=1.0e-300;
     }

     /* initialize WC[k][i][j]=0 */
     /* initialize WD[k][i][j]=0 */
     /* initialize WR[k][i][j]=0 */
     /* for k!=0 and i=0 or j=0  */
     for(k=1;k<work_blocks;k++){
         for(i=0;i<=Len0;i++){
	     SK->WC[k][i][0]=0;
	     SK->WD[k][i][0]=0;
	     SK->WR[k][i][0]=0;
	 }
	 wc_ptr=SK->WC[k][0];
	 wd_ptr=SK->WD[k][0];
	 wr_ptr=SK->WR[k][0];
	 for(j=0;j<=Len1;j++){
	     wc_ptr[j]=0.0;
	     wd_ptr[j]=0.0;
	     wr_ptr[j]=0.0;
	 }
     }
}

/************ alloc_Sank_BW  ************************/
/*                                                  */
/*  Allocate space for BW                           */
/****************************************************/
     
void alloc_Sank_BW(Model B, Sank_S SK)
{
     register int i,j,k,n;
     int      Len0,Len1;                  /* for quick access */
     double   *bwd_ptr,*bwr_ptr,*bwc_ptr; /* for quick access */
             
     Len0=B->Seq->SeqLen[0];
     Len1=B->Seq->SeqLen[1];

     /* allocate space for WD,WR and WC */
     NEWPP(SK->BWD,SK->Sankoff_blocks+1,double);
     NEWPP(SK->BWR,SK->Sankoff_blocks+1,double);
     NEWPP(SK->BWC,SK->Sankoff_blocks+1,double);

     for(i=0;i<=SK->Sankoff_blocks;i++) {
         NEWP(SK->BWD[i],Len0+2,double);
	 NEWP(SK->BWR[i],Len0+2,double);
	 NEWP(SK->BWC[i],Len0+2,double);
	 for(j=0;j<=Len0+1;j++) {
	     NEW(SK->BWD[i][j],Len1+2,double);
	     NEW(SK->BWR[i][j],Len1+2,double);
	     NEW(SK->BWC[i][j],Len1+2,double);
	 }
     }

     /* initialize the BWC[0][i][j]=0 */
     /* BWD[0][i][j]=0, BWR[0][i][j]=1 */
     for(i=1;i<=Len0+1;i++){
         bwc_ptr=SK->BWC[0][i];
	 bwd_ptr=SK->BWD[0][i];
	 bwr_ptr=SK->BWR[0][i];

	 for(j=1;j<=Len1+1;j++){
	     bwc_ptr[j]=0.0;
	     bwd_ptr[j]=0.0;
	     bwr_ptr[j]=1.0e-300;	
	 }
     }

     /* initialize the BWC[0][i][Len1+1]=0, BWC[0][Len0+1][j]=0 */
     /* initialize the BWD[0][i][Len1+1]=1, BWD[0][Len0+1][j]=0 */
     /* initialize the BWR[0][i][Len1+1]=0, BWR[0][Len0+1][j]=1 */
     for(i=1;i<=Len0+1;i++){
         SK->BWC[0][i][Len1+1]=0;
	 SK->BWD[0][i][Len1+1]=1.0e-300;
	 SK->BWR[0][i][Len1+1]=0;
     }
     bwc_ptr=SK->BWC[0][Len0+1];
     bwd_ptr=SK->BWD[0][Len0+1];
     bwr_ptr=SK->BWR[0][Len0+1];

     for(j=0;j<=Len1+1;j++){
         bwc_ptr[j]=0;
	 bwd_ptr[j]=0;
	 bwr_ptr[j]=1.0e-300;
     }

     for(k=1; k<=SK->Sankoff_blocks; k++) {
         for(i=1; i<=Len0+1; i++) {
	     bwc_ptr=SK->BWC[k][i];
	     bwd_ptr=SK->BWD[k][i];
	     bwr_ptr=SK->BWR[k][i];
	     for(j=0; j<=Len1+1; j++) {
	         bwc_ptr[j]=0.0;
		 bwd_ptr[j]=0.0;
		 bwr_ptr[j]=0.0;
	     }
	 }
     }	     
}

/************ alloc_Sank_rel_matrix  ****************/
/*                                                  */
/*  Allocate space for relation matrix              */
/****************************************************/
     
void alloc_Sank_rel_matrix(Model B, Sank_S SK)
{
    register   int i,k;
    int        AlphaLen,Nmatrix;
    double     ***pp_related;
    float      ***ff_related;
    char       **matrix_name;

    Nmatrix=SK->rel_matrix.Nmatrix;
    AlphaLen=B->IP->nAlphaLen;

    NEWPP(pp_related,Nmatrix,double);
    NEWPP(ff_related,Nmatrix,float);
    NEWP(matrix_name,Nmatrix,char);
    for(k=0; k<Nmatrix; k++) {
        NEWP(pp_related[k],AlphaLen,double);
        NEWP(ff_related[k],AlphaLen,float);
	for(i=0; i<AlphaLen; i++) {
	    NEW(pp_related[k][i],AlphaLen,double);
	    NEW(ff_related[k][i],AlphaLen,float);
	}
	NEW(matrix_name[k],30,char);
    }
    SK->rel_matrix.pp_related=pp_related;
    SK->rel_matrix.ff_related=ff_related;
    SK->rel_matrix.matrix_name=matrix_name;
}


/*****************  init_maxdata ********************/

void init_maxdata(MaxResults *tmp)
{
   tmp->F = NULL;
   tmp->nNumMotifs = NULL;
   tmp->nMotifLoc = NULL;
   tmp->RevComp = NULL;
   tmp->dvMotifProb = NULL;
   tmp->frequency = NULL;
}

/***************  FreeData(Model B) **************************/
/*                                                           */
/*       free space in a model                               */
/*************************************************************/

void FreeData(Model B)
{
   int i, j, t, maxlen;

   maxlen = findMaxMotifLen(B->IP)+1;

   FREEPP(B->C->dSumPseudo,  B->IP->nNumMotifTypes,2);
   FREEPP(B->First->dSumPseudo,  B->IP->nNumMotifTypes,2);
   FREEPP(B->C->dPseudoCounts, B->IP->nNumMotifTypes,maxlen);
   FREEPP(B->First->dPseudoCounts, B->IP->nNumMotifTypes,maxlen);
   FREEP(B->First->dmodel_sites, B->IP->nNumMotifTypes);
   FREEP(B->C->dmodel_sites, B->IP->nNumMotifTypes);
   FREEP(B->First->dmodel_pseudo, B->IP->nNumMotifTypes);
   FREEP(B->C->dmodel_pseudo, B->IP->nNumMotifTypes);
   free(B->First->dtot_sites);
   free(B->C->dtot_sites);
   free(B->First->dtot_pseudo);
   free(B->C->dtot_pseudo);
   free(B->First->dTot);
   free(B->C->dTot);
   free(B->First->dbg_pseudo);
   free(B->C->dbg_pseudo);
   FREEPP(B->First->fCounts, B->IP->nNumMotifTypes,maxlen);
   FREEPP(B->C->fCounts, B->IP->nNumMotifTypes,maxlen);
   if(B->IP->Datafiles->output_filename != NULL)
      free(B->IP->Datafiles->output_filename);
   if(B->IP->Datafiles->prior_filename != NULL)
      free(B->IP->Datafiles->prior_filename);
   free(B->IP->Datafiles);
   free(B->IP->dposterior_prob);
   free(B->IP->nMotifLen);
   free(B->IP->DOF);
   free(B->IP->nPossSites);
   for(t = 0; t < B->IP->nNumMotifTypes; t++) {
      free(B->IP->nNumMotifs[t]);
      if(B->IP->AltModel->Collapsed != NULL)
         free(B->IP->AltModel->Collapsed[t]);
      if(B->IP->AltModel->Palandromic != NULL)
         free(B->IP->AltModel->Palandromic[t]);
      if(B->IP->AltModel->Concen != NULL)
         free(B->IP->AltModel->Concen->Concentrated[t]);
      if(B->F != NULL)  {
         free(B->F->nColMask[t]);
         free(B->F->nOldColMask[t]);
      } 
   }
   if(B->IP->AltModel->Collapsed != NULL)
      free(B->IP->AltModel->Collapsed);
   if(B->IP->AltModel->Palandromic != NULL)
      free(B->IP->AltModel->Palandromic);
   if(B->IP->AltModel->Concen != NULL) {
      free(B->IP->AltModel->Concen->Concentrated);
      free(B->IP->AltModel->Concen->NumConcen);
      free(B->IP->AltModel->Concen);
   }
   if(B->F != NULL) {
      free(B->F->nColMask);
      free(B->F->nOldColMask);
      free(B->F->nMaxLen);
      free(B->F->FragWidth);
      free(B->F->shift);
      free(B->F); 
   }
   free(*B->Seq->R);
   free(*B->Seq->ProcessedSTR);
   free(B->IP->nNumMotifs);
   free(B->IP->Datafiles->desc_filename);
   free(B->Seq->ProcessedSTR);
   free(B->Seq->nvEndLocs);
   free(B->Seq->R);
   /* free individual sequences */
   for(i=0;i<=B->IP->nNumSequences-1;i++)
     free(B->Seq->Orig[i]);
   free(B->Seq->Orig);
   free(B->Seq->SeqLen);
   free(B->First);
   free(B->C); 
   free(B->Seq);
   free(B->IP);
   free(B);
   B = NULL;
}



/************************** free_maxdata  ************************/

void free_maxdata(MaxResults *tmp, IPtype IP)
 
{
   int maxlen=0, i, j;

   if(tmp->F != NULL) {
      free(tmp->F->nMaxLen);
      free(tmp->F->FragWidth);
      free(tmp->F->shift);
      if(tmp->F->nColMask != NULL)
         FREEP(tmp->F->nColMask, IP->nNumMotifTypes);
      if(tmp->F->nOldColMask != NULL)
         FREEP(tmp->F->nOldColMask, IP->nNumMotifTypes);
      free(tmp->F);
   }
   if(tmp->nNumMotifs != NULL) {
      for(i = 0; i < IP->nNumMotifTypes; i++) {
         if(tmp->nNumMotifs[i] > maxlen)
            maxlen = tmp->nNumMotifs[i];
      } 
      free(tmp->nNumMotifs);
   }
   if(tmp->nMotifLoc != NULL)
      FREEP(tmp->nMotifLoc, maxlen);
   if(tmp->RevComp != NULL)
      FREEP(tmp->RevComp, maxlen);
   if(tmp->dvMotifProb != NULL)
      FREEP(tmp->dvMotifProb, maxlen);
   if(tmp->frequency != NULL)
      FREEP(tmp->frequency, IP->nSeqLen);
}

/****************** free_Sank***********************/

void free_Sank(Model B,Sank_S SK)
{
     int i,k;

     for(k=0;k<=SK->Sankoff_blocks;k++){
         free(SK->start_x[k]);
	 free(SK->start_y[k]);
	 free(SK->end_x[k]);
	 free(SK->end_y[k]);
     }
     free(SK->start_x);
     free(SK->start_y);
     free(SK->end_x);
     free(SK->end_y);

     if(SK->rel_matrix.Nmatrix>1){ /* use multiple matrices */
         for (i=0;i<=SK->rel_matrix.Nmatrix-1;i++){
	     free(SK->PKSI[i]);	 
	     free(SK->W[i]);
	 }
     }
     else {/* for a single matrix */
         free(SK->PKSI[0]);	 
	 free(SK->W[0]);
     }
     free(SK->PSI);       
     free(SK->PKSI); 
     free(SK->PK);
     free(SK->W);
     free(SK->N);
     free(SK->prior_N);

}

/****************** free_Sank_S ***********************/

void free_Sank_S(Model B,Sank_S SK)
{
  
     register int i,j,k,n;

     /* free score space */
     for(n=0;n<SK->rel_matrix.Nmatrix;n++){
         for(k=0;k<=SK->Sankoff_blocks;k++){
	     for(i=0;i<=B->Seq->SeqLen[0];i++){
	         free(SK->SV[n][k][i]);
	         free(SK->SW[n][k][i]);
	     }
	     free(SK->SV[n][k]);
	     free(SK->SW[n][k]);
	 }
	 free(SK->SV[n]);
	 free(SK->SW[n]);
     }
     free(SK->SV);
     free(SK->SW);

}

/****************** free_Sank_N ***********************/

void free_Sank_N(Model B,Sank_S SK)
{
  
     register int i,k;

     /* free space for ND, NR and NC*/
     for(k=0;k<=1;k++){
	 for(i=0;i<=B->Seq->SeqLen[0];i++){
	     free(SK->ND[k][i]);
	     free(SK->NR[k][i]);
	     free(SK->NC[k][i]);
	 }
	 free(SK->ND[k]);
	 free(SK->NR[k]);
	 free(SK->NC[k]);
     }

     free(SK->ND);
     free(SK->NR);
     free(SK->NC);


}


/****************** free_Sank_W ***********************/

void free_Sank_W(Model B,Sank_S SK)
{
  
     register int i,k,n;
     int      work_blocks;

     if(SK->flags.back_sampling || SK->flags.exact_posterior){
         /* allocated full space */
         work_blocks=SK->Sankoff_blocks+1;
     }
     else{/* only two array of working space allocated */
         work_blocks=2;
     }
  
      if(SK->flags.min_memory && !SK->flags.exact_posterior)
         work_blocks=2; 
   
     /* free space for WD, WR, WC */     
     for(k=0;k<work_blocks;k++){
         for(i=0;i<=B->Seq->SeqLen[0];i++){
	     free(SK->WD[k][i]);
	     free(SK->WR[k][i]);
	     free(SK->WC[k][i]);
	 }
	 free(SK->WD[k]);
	 free(SK->WR[k]);
	 free(SK->WC[k]);
     }

     free(SK->WD);
     free(SK->WC);
     free(SK->WR);
}

/************ free_Sank_BW  *************************/
/*                                                  */
/*  Free space for BW                               */
/****************************************************/
     
void free_Sank_BW(Model B, Sank_S SK)
{
     register int i,j,k,n;
     int      Len0;                       /* for quick access */
             
     Len0=B->Seq->SeqLen[0];

     for(k=0;k<=SK->Sankoff_blocks;k++){
	 for(i=0;i<=Len0+1;i++){
	     free(SK->BWD[k][i]);
	     free(SK->BWR[k][i]);
	     free(SK->BWC[k][i]);
	 }
         free(SK->BWD[k]);
	 free(SK->BWR[k]);
	 free(SK->BWC[k]);
     }
     free(SK->BWD);
     free(SK->BWR);
     free(SK->BWC);     
}

/************ free_Sank_rel_matrix  *****************/
/*                                                  */
/*  Free space for relation matrix                  */
/****************************************************/
     
void free_Sank_rel_matrix(Model B, Sank_S SK)
{
    register int   i,k;
    int            AlphaLen,Nmatrix;
    double         ***pp_related;
    float          ***ff_related;
    char           **matrix_name;

    pp_related=SK->rel_matrix.pp_related;
    ff_related=SK->rel_matrix.ff_related;
    matrix_name=SK->rel_matrix.matrix_name;
    Nmatrix=SK->rel_matrix.Nmatrix;
    AlphaLen=B->IP->nAlphaLen;

    for(k=0; k<Nmatrix; k++) {
	for(i=0; i<AlphaLen; i++) {
	    free(pp_related[k][i]);
	    free(ff_related[k][i]);
	}
        free(pp_related[k]);
	free(matrix_name[k]);
    }
    free(pp_related);
    free(ff_related);
    free(matrix_name);
}    










