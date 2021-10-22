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

/**********************************************************************/
/*                                                                    */
/* $Id: Gibbs_NON_GUI.c,v 1.3 1997/03/18 03:23:12 junzhu Exp junzhu $ */
/*                                                                    */
/* Author: Jun Zhu                                                    */
/* Date:   July 1, 1996                                               */
/* Description: non_graphic user interface for Gibbs,                 */
/*              working for batch job                                 */
/**********************************************************************/

#include "routines.h"
#include "mem_mgmt.h"


int Gibbs_NON_GUI(int argc, char **argv)
{
  int    i,j,n;
  int    matrix_start,matrix_end;  /* start and end matrix # for DNA */
  int    matrix_step;              /* for DNA */
  Model  M;
  Sank_S SK;
  char   *cval;
  char   option[30];
  char   *queryfile,*datafile,*outputfile;
  FILE   *fptr=NULL;
  int    q_file=0;                          /* flag for query file */
  int    d_file=0;
  char   tmp_filename[] = "/tmp/tmpXXXXXX";
  
  /* allocate space for GUI_argv, we only allocate array of pointer,    */
  /* we do not need to locate space for each array pointed by pointers. */
  /* The pointers will accept values from other array, such as          */
  /* GUI_argv[i]=argv[i+1]                                              */
  /* this is for plug-in interface of Bernoulli program                 */
  NEWP(GUI_argv,30,char);


  /* Sankoff (Bayesian) method */
      M = alloc_Model();
      NEW(M->IP->Datafiles,1,files);
      M->IP->Datafiles->out_fpt=stdout;
      NEW(M->IP->Datafiles->desc_filename, 256, char);
      
      /* M->IP->Datafiles->desc_filename = tempnam("/tmp", NULL); */
      int fd = mkstemp(tmp_filename);
      close(fd);
      M->IP->Datafiles->desc_filename = tmp_filename;
      
      M->IP->nAlphaLen=20;  /* protein sequence */

      SK=(Sank_S)malloc(sizeof(Sankoff_struct));
      /* default values */
      SK->flags.best_alignment=FALSE;
      SK->flags.output_sequence=FALSE;
      SK->Sankoff_blocksize=BASE_BLOCKSIZE;
      SK->flags.back_sampling=TRUE;

      SK->backsampling_number=1000;
      SK->flags.back_sampling_sequence=FALSE;
      SK->flags.back_sampling_cutoff=FALSE;
      SK->backsampling_cutoff_value=1e10;
      SK->flags.debug=FALSE;
      SK->flags.exact_posterior=FALSE;
      SK->flags.protein_sequence=FALSE;
      SK->flags.min_memory=TRUE;
      SK->flags.back_sampling_proportional = TRUE;
      SK->Sankoff_max_blocks = DEFAULT_NUM_BLOCKS;
      SK->backsampling_blocks=SK->Sankoff_max_blocks;
      SK->flags.profile=TRUE;
      SK->flags.full_profile = FALSE;

      SK->rel_matrix.matrix_num=5;  /* use blosum62 matrix */
      NInput_N=0;                   /* no special sites */
      Input_alpha=0;                /* no special site weight */
      SK->flags.protein_sequence=1; /* protein sequence */


      /* scan command line options */
      for(n=1;n<argc;n++) {
	  i=0;
	  while((argv[n][i]!='=' ||argv[n][i]==(char) 0)&& i < 29) {
	      option[i]=argv[n][i];
	      i++;
	 }
	 /* add terminator the end of option string */
	 option[i]=(char) 0;
	 i++;
         if(strcmp(option, "-PBayes_align")==0) {
	   /* ignore this obsolete option */
         }
	 else if(strcmp(option,"-QUERY")==0) { /* input data file */
	     queryfile=&argv[n][i];
	     if((fptr=fopen(queryfile,"r"))==NULL) {
	         printf("Check qury file name!!!\n");
		 return 0;
	     }
	     SK->queryfile_ptr=fptr;
	     q_file=1;
	 }
	 else if(strcmp(option,"-DATA")==0) {/*  data file */
	     datafile=&argv[n][i];
	     if((fptr=fopen(datafile,"r"))==NULL){
	         printf("Check data file name!!!\n");
		 return 0;
	     }
	     M->IP->Datafiles->fpt=fptr;
	     M->IP->Datafiles->filename=
	       (char *)malloc((strlen(datafile)+1)*sizeof(char));
	     strcpy(M->IP->Datafiles->filename,datafile);
	     d_file=1;
	 }
	 else if(strcmp(option,"-OUTFILE")==0) { /* output file */
	     outputfile=&argv[n][i];
	     M->IP->Datafiles->out_fpt=fopen(outputfile,"w");
	     if(M->IP->Datafiles->out_fpt==NULL)
	         M->IP->Datafiles->out_fpt=stdout;
	 }
	 else if(strcmp(option,"-BASE_BLOCKSIZE")==0) {
	     /* Maximum Number of Blocks */
	     cval=&argv[n][i];
	     sscanf(cval,"%d", &SK->Sankoff_blocksize);
	 }
	 else if(strcmp(option,"-MAX_NUM_BLOCKS")==0) {
	     /* Maximum Number of Blocks */
	     cval=&argv[n][i];
	     sscanf(cval,"%d", &SK->Sankoff_max_blocks);
             SK->backsampling_blocks=SK->Sankoff_max_blocks;
	 }

	 else if(strcmp(option,"-NON_PROTEIN")==0) {
	     /* protein sequence flag */
	     SK->flags.protein_sequence=FALSE;
	     M->IP->nAlphaLen=5;
	 }
         else if(strcmp(option,"-BACK_S_PROPORTIONAL")==0) {
	   /* back sample proportional */
	   /*           SK->flags.back_sampling_proportional=TRUE; */
         }
         else if(strcmp(option,"-EXACT_BLOCK")==0) {
	   /* don't back sample proportional */
             SK->flags.back_sampling_proportional=FALSE;
         }
         else if(strcmp(option,"-EXACT_POST")==0) {
	   /* don't back sample */
             SK->flags.exact_posterior=TRUE;
             SK->flags.back_sampling=FALSE;
             SK->flags.back_sampling_proportional=FALSE;
             SK->flags.min_memory=FALSE;
         }

	 else if(strcmp(option,"-MIN_MEMORY")==0) {
	     /*minimum memory hack */
	   /*	     SK->flags.min_memory=TRUE; */
	 }
	 else if(strcmp(option,"-MAX_SPEED")==0) {
	     /*speed  hack */
	     SK->flags.min_memory=FALSE;
	 }
	 else if(strcmp(option,"-IMPORT_SITE")==0) {
	     /* important sites in query sequence */
	     cval=&argv[n][i];
	     /* scan sites */
	     NInput_N=0;
	     j=0;
	     for(i=0;i<=strlen(cval);i++){
	         if(cval[i]==',' || cval[i]==(char ) 0){
		     /* add terminator the end of option string */
		     option[j]=(char) 0;
		     sscanf(option,"%d", &Input_N[NInput_N]);
		     NInput_N++;
		     j=0;
		 }
		 else{
		     option[j]=cval[i];
		     j++;
		 }
	     }
	     printf("NInput_N=%d\n",NInput_N);
	 }
	 else if(strcmp(option,"-SITE_WEIGHT")==0) {
	     /* Maximum Number of Blocks */
	     cval=&argv[n][i];
	     sscanf(cval,"%f", &Input_alpha);
	 }
	 else if(strcmp(option,"-MATRIX_NUM")==0) {
	     /* select matching matrix */
	     cval=&argv[n][i];
	     sscanf(cval,"%d", &SK->rel_matrix.matrix_num);
	     if(SK->rel_matrix.matrix_num >MAX_NUM_MATRIX_B+MAX_NUM_MATRIX_P+1
		 ||SK->rel_matrix.matrix_num <0)
	     {
	         printf("matrix_num should be [0-36].\n");
		 printf("matrix_num[0-7] BLOSUM30-100.\n");
		 printf("matrix_num[8] ALL BLOSUM.\n");
		 printf("matrix_num[9-35] PAM40-300, step 10.\n");
		 printf("matrix_num[36] ALL PAM.\n");
		 exit(-1);
	     }
	     if(SK->rel_matrix.matrix_num==MAX_NUM_MATRIX_B)
	         SK->rel_matrix.Nmatrix=MAX_NUM_MATRIX_B;
	     else if(SK->rel_matrix.matrix_num==MAX_NUM_MATRIX_B+
		     MAX_NUM_MATRIX_P+1)
	         SK->rel_matrix.Nmatrix=MAX_NUM_MATRIX_P;
	     else
	         SK->rel_matrix.Nmatrix=1;

	     if(SK->rel_matrix.matrix_num>MAX_NUM_MATRIX_B)
	         /* do not count all-blosum option */
	         SK->rel_matrix.matrix_num--;
	 }
	 else if(strcmp(option,"-DNA_MATRIX")==0) {
	     /* select DNA  matching matrix */
	     cval=&argv[n][i];
	     sscanf(cval,"%d-%d,%d", &matrix_start,&matrix_end,&matrix_step);
	     if(matrix_start >matrix_end || matrix_start <1
		 || matrix_end >500)
	     {
	         printf("dna_matrix_num should be [1-500].\n");
		 exit(-1);
	     }
	     if(matrix_step<0) {
	         printf("matrix_step should be greater or equal 1.\n");
		 exit(-1);
	     }

	     /* relation matrix */
	     SK->rel_matrix.matrix_num=matrix_start;
	     SK->rel_matrix.index_start=matrix_start;
	     SK->rel_matrix.index_step=matrix_step;
	     SK->rel_matrix.Nmatrix=(matrix_end-matrix_start)/matrix_step+1;
	 }
	 else if(strcmp(option,"-DEBUG")==0) { /* debug flag */
	     SK->flags.debug=TRUE;
	 }
	 else if(strcmp(option,"-OUT_SEQ")==0) {
	     /* output matching sequence flag */
	     SK->flags.output_sequence=TRUE;
	 }
	 else if(strcmp(option,"-BEST_ALIGN")==0) {
	     /* calculate best alignment */
	     SK->flags.best_alignment=TRUE;
	 }
	 else if(strcmp(option,"-BACK_SAMPLE")==0) {
	     /* back sampling flag */
	   /*	     SK->flags.back_sampling=TRUE; */
	 }
      	 else if(strcmp(option,"-BACK_S_SIZE")==0) {
	     /* back sampling number */
	     cval=&argv[n][i];
	     sscanf(cval,"%d", &SK->backsampling_number);
	 }
	 else if(strcmp(option,"-BACK_S_OUTSEQ")==0) {
	     /* output back sampling matching blocks flag */
	     SK->flags.back_sampling_sequence=TRUE;
	 }
	 else if(strcmp(option,"-BACK_CUTOFF")==0) {
	     /* find all alignment above a cutoff flag */
	     SK->flags.back_sampling_cutoff=TRUE;
	 }
	 else if(strcmp(option,"-BACK_C_VALUE")==0) {
	     /* cutoff value */
	     cval=&argv[n][i];
	     sscanf(cval,"%f",&SK->backsampling_cutoff_value);
	 }
         else if(strcmp(option,"-NO_PROFILE")==0) {
              SK->flags.profile=FALSE;
	 }
	 else if(strcmp(option,"-FULL_PROFILE")==0) {
	     /* protein sequence flag */
	     SK->flags.full_profile=TRUE;
         }
	 else {
	     printf("Illegal option %s \n",option);
	     /* print out the syntax */
	     print_usage_sankoff();
	     exit(-1);
	 }
      }

      /* check mandatory input */
      if(q_file==0) { /* no input file */
	  /* print out the syntax */
	  printf("NO query file!\n");
	  printf("Syntax for non-graphic users:\n");
	  printf("Gibbs -Pprogram_name [options]\n");
	  printf("for Sankoff method:\n");
	  print_usage_sankoff();
	  exit(-1);
      }

      if(d_file==0) { /* no input file */
	  /* print out the syntax */
	  printf("NO data file!\n");
	  printf("Syntax for non-graphic users:\n");
	  printf("Gibbs -Pprogram_name [options]\n");
	  printf("for Sankoff method:\n");
	  print_usage_sankoff();
	  exit(-1);
      }
      if (SK->flags.min_memory &&  SK->flags.exact_posterior) {
          printf("MIN_MEMORY can not be used with Exact option\n");
          printf("MIN_MEMORY option ignored\n");
          SK->flags.min_memory=FALSE;
      }

      if (SK->flags.min_memory &&  SK->flags.best_alignment) {
          printf("MIN_MEMORY can not be used with BEST_ALIGN option\n");
          printf("MIN_MEMORY option ignored\n");
          SK->flags.min_memory=FALSE;
      }
      /* the only intelligent option */
      if (SK->flags.back_sampling_proportional==FALSE) {
          SK->backsampling_blocks=SK->Sankoff_max_blocks;
      }

      Sankoff(M,SK);
	  return(0);

}



