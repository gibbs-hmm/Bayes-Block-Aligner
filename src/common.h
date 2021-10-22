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
/* $Id: common.h,v 1.4 1996/10/25 19:32:57 junzhu Exp junzhu $            */
/*                                                                        */
/* Author:       Eric C. Rouchka,1996,7.                                  */
/*               Jun Zhu, 1996,8.                                         */
/**************************************************************************/

#ifndef COMMON_H
#define COMMON_H
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <time.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>
#include "Sankoff.h"


/*========================================================================*/
/*=================================CONSTANTS==============================*/
/*========================================================================*/
#ifndef INT32_MAX
#define INT32_MAX  2147483647
#endif
#define DEBUG 1
#define MAX_WIDTH_MULT 5
#define PLATEAU_PERIODS 20
#define NUM_SEEDS 10
#define MAX_ITERATIONS 500
#define PSEUDO_CNT_WT 0.1
#define PSEUDO_SITE_WT 0.8
#define NEIGHBOR_RESIDUES 5
#define NEAROPT_CUT 0.5
#define MIN_PSEUDOCNT 0.001
#define FALSE 0
#define TRUE 1


/*========================================================================*/
/*==================================MACROS================================*/
/*========================================================================*/

#define MEW(x,n,t)      (( (x=(t*) malloc(((n)*sizeof(t))))==NULL) ? \
                         (t*) (intptr_t) (fprintf(stderr,"Out of Memory."),exit(1),0):x)
#define NEW(x,n,t)      (( (x=(t*) calloc(n,sizeof(t)))==NULL) ? \
                         (t*) (intptr_t) (fprintf(stderr,"Out of Memory."),exit(1),0):x)
#define NEWP(x,n,t)     (( (x=(t**) calloc(n,sizeof(t*)))==NULL) ? \
                        (t**) (intptr_t) (fprintf(stderr,"Out of Memory."),exit(1),0):x)

#define NEWPP(x,n,t)    (( (x=(t***) calloc(n,sizeof(t**)))==NULL) ? \
                        (t***) (intptr_t) (fprintf(stderr,"Out of Memory."),exit(1),0):x)

#define NEWP3(x,n,t)    (( (x=(t****) calloc(n,sizeof(t***)))==NULL) ? \
                        (t****) (intptr_t) (fprintf(stderr,"Out of Memory."),exit(1),0):x)

/* Note for FREEP need local variable i defined  */
/* and for FREEPP need variables i and j defined */

#define FREEP(x, n) 	{for(i=0;i < (n);i++) \
                           if((x)[i] != NULL) free((x)[i]);\
                        if((x) != NULL) free((x));(x) = NULL;}

#define FREEPP(x,ii,jj)	{for(i=0;i <(ii);i++){\
			   for(j=0;j<(jj);j++) \
                             free((x)[i][j]);\
                           free((x)[i]);}\
			(x)=NULL;}



#define NUMMOTIFS(x)  ((x)[FORWARD] + (x)[REVERSE])

#define CH2INT(n,B)	(int)(*(B)->Seq->R)[(n)] - 97

#define CH2INTCOMP(n,B)	(int)(complement((*(B)->Seq->R)[(n)])) - 97

#define GETSTRING(S) i = 0; Finished = FALSE;\
                     while(!Finished) {\
                        scanf("%c", &(S)[i]);\
                        if((S)[i] == '\n') {\
                           Finished = TRUE;\
                           (S)[i] = '\0';\
                        }\
                        i++;\
                     }


#define max(a, b)	((a) > (b) ? (a) : (b))
#define min(a, b)	((a) < (b) ? (a) : (b))


/*================== define enumberate constant =======================*/

enum flags     {DELETE, ADD};
enum where     {BG, MOTIF};
enum indicator {POSSIBLE, ENDS, MARKED};
enum direction {FORWARD, REVERSE};
enum col_frag  {COL_OFF, COL_ON, COL_BLOCKED};
enum cl_params {cl_a, cl_c, cl_e, cl_h, cl_i, cl_m, cl_o, cl_p, cl_r,cl_s,
		cl_t,cl_w, cl_x, cl_C, cl_F, cl_I, cl_P, cl_R, cl_S, cl_W};

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/* cl_a : Concentrated Region           cl_t : Sequence type         */
/* cl_c : Collapsed Alphabet            cl_w : Pseudocount weight    */
/* cl_e : Use Expectation/Maximization  cl_x : Don't Xnu sequence    */
/* cl_h : Help flag                     cl_C : Near optimal cutoff   */
/* cl_i : Number of iterations          cl_F : Don't fragment        */
/* cl_m : Don't use map maximization    cl_I : Interactive mode      */
/* cl_o : Output file                   cl_P : Informed priors file  */
/* cl_p : Plateau periods               cl_R : palandromic sequence  */
/* cl_r : Reverse complement            cl_S : Numer of seeds        */
/* cl_s : Seed Value                    cl_W : Pseudosite weight     */
/*                                                                   */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

/*=======================================================================*/
/*================================DATA TYPES=============================*/
/*=======================================================================*/
typedef struct M1{
            int Mtype;
            int pos;
            int seq_num;
            int left, right;
            char *sMotif;
            short RevComp;
            struct M1 *next;
   } MlistEl;

typedef struct {
            int nNumMotifs;
            int nMotifLen;
            MlistEl *Motifs;
   } Mlist_struct;

typedef Mlist_struct **Mlist;

typedef struct {
            int nPossStartPos;
            int nInMotif;
            int nMotifStartPos;
            short RevComp;               /* Is Motif Reverse complement?*/
   } PoSition;

typedef struct {
            double dResidueInBGProb;
            double dResidueInMotifProb;
            double dinBGProb;
            double dinMotifProb;
            double dBGProbDenom;
            double dMotifProbFact;
            double ***dvInMotifProb;          /* For Each individual residue */
            double **dvInBGProb;
            unsigned short update;
   } ProbStruct;

typedef struct {
            float  ***fCounts;              /* Observed count      */
                                            /* fCounts[Motiftype]  */
                                            /*        [PositioninMotif] */
                                            /*        [Alphabet]    */
            double ***dPseudoCounts;        /* Pseudo Counts        */
                                            /* dPseudoCounts[Motiftype]  */
                                            /*        [PositioninMotif] */
                                            /*        [Alphabet]    */
            double ***dSumPseudo;           /* sum pseudo counts    */
                                            /* dSumPseudo[Motiftype]*/
					    /*        [Background/Motif]  */
                                            /*        [Alphabet]    */
                                            /* arrange memory in this way */
                                            /* is not very efficient, but */
                                            /* it is clear when used with */
                                            /* above arrays               */

            double dSumBGPseudo;
            int    nTotBack;                /* Total Observed BG Residues    */
            int    nTotMotif;               /* Total Observed Motif Residues */
            double dTotSumMotifPseudo;

            double **dmodel_sites;          /* See protein science paper */
            double **dmodel_pseudo;         /* Appendix and p 1627 (3)   */
            double *dtot_sites;             /* for the description       */
            double *dtot_pseudo;
            double *dbg_pseudo;
            double *dTot;

   } Counts;

typedef Counts *Ctype;

typedef struct {
            FILE *fpt;                          /* sequence file pointer    */
            FILE *prior_fpt;                    /* beta priors file pointer */
            FILE *out_fpt;                      /* output file pointer      */
            FILE *desc_file;

            char *filename;
            char *prior_filename;
            char *output_filename;
            char *desc_filename;
   } files;

typedef struct {
            int* SeqLen;            /* length of each sequence  */
            short **Orig;           /* sequence in number       */
            char **R;               /* all strings concatenated */
            char **ProcessedSTR;    /* sequence after Xnu       */
            int  **nvEndLocs;       /* end location of each sequence
				       in  **R                  */

   } Stringstruct;

typedef Stringstruct *Stype;

typedef struct {
   short **Concentrated;
   short *NumConcen;
} Concenstruct;

typedef Concenstruct *ConType;

typedef struct {
   int **Collapsed;        /* Array indicating positions that are collapsed   */
   int **Palandromic;      /* Array indicating positions that are palandromic */
   ConType Concen;
} AltModelStruct;

typedef AltModelStruct *AltModelType;


typedef struct  {
            short is_defined[30];
            short site_samp;
            int *nMotifLen;              /* Motif Residue Lengths       */
            int nSeqLen;                 /* Sequence Residue Length     */
            int nAlphaLen;               /* # of characters in alphabet */
            int **nNumMotifs;            /* number of initial motifs    */
            int nNumSequences;           /* number of different seq.    */
            files *Datafiles;            /* datafiles structure         */
            double dPseudoCntWt;
            double dPseudoSiteWt;
            int nMaxIterations;
            long lSeedVal;
            int nPlateauPeriods;
            int nSeeds;
            double *dposterior_prob;
            int *nPossSites;
            double dnull_map;
            int col_shift;
            int nNumMotifTypes;
            double dCutoff;
            short RevComplement;
            int nNumProcessed;
            int *DOF;                   /* degree of freedom for each
					   motif type                    */
            AltModelType AltModel;
            double Logl;
   } InputParams;

typedef InputParams *IPtype;

typedef struct {
            int *nMaxLen;              /* 5 * width of motif         */
            int **nColMask;            /* Column Mask for each motif */
            int **nOldColMask;
            int *FragWidth;
            double ***nvFragCnts;
            int *shift;                /* how far first column has shifted */
   } FragStruct;

typedef FragStruct *Ftype;

typedef struct {
            IPtype      IP;
            Ctype       C;
            Ctype       First;
            ProbStruct  *P;
            files       *Datafiles;
            Ftype       F;
            Stype       Seq;
   } Modelstruct;

typedef Modelstruct * Model;

typedef struct {
            double frequency;              /* sampling frequency of motif  */
            double rev_freq;               /* Reverse Complement frequency */
            double fwd_freq;               /* Forward frequency */
   } Frequency;

typedef struct {
            double dProbability;
            int *nNumMotifs;
            int **nMotifLoc;
            short **RevComp;
            long nseed;
            int nIterationNum;
            double **dvMotifProb;
            Frequency **frequency;           /* frequency structure */
            Ftype F;
   } MaxResults;

typedef struct {
            long begin;
            long end;
            long total;
   } SimTime;

typedef struct {
            double abs_logl;
            double pct_logl;
            double abs_tot;
            double pct_tot;
            double last_logl;
            double last_tot;
            double new_tot;
            double logl;
            double likelihood;
            double cost;
   } EMStats;

typedef struct MC {
   char *desc;
   char *cmp_desc;
   double logl;
   double cmp_logl;
   double dof;
   double lrt;
   double pval;
   struct MC *next;
}  MCstruct;

/*-------------------------------------------------------------------------*/
/*-------------------------------FUNCTION PROTOTYPES-----------------------*/
/*-------------------------------------------------------------------------*/
short SEQ2CH(char,short);
void  print_usage_bernoulli (void);
void  print_usage_sankoff(void);
/*void  print_parameters_sankoff(FILE* fptr, Sank_S SK); */
void  p_error (char *msg);
void  p_usage_error (char *msg);
int   findMaxNumMotif (IPtype IP);
int   findMaxMotifLen (IPtype IP);
char  complement (char ch);
void  sRandom (long x);
long  Random (void);
void  rm_file(char *fname);
void  free_maxdata(MaxResults *tmp, IPtype IP);
void  FreeData(Model);
void  print_EnterArgs(void);
void  BeginTime(SimTime *S);
void  EndTime(SimTime *S);
void  FreeData(Model B);
void print_parameters_sankoff(FILE *fptr, Sank_S SK);

#endif
