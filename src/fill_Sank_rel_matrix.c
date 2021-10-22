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

/*************************************************************************/
/* $Id: fill_Sank_rel_matrix.c,v 1.1 1998/02/06 19:46:27 junz Exp junz $ */
/*                                                                       */
/* Author: Jun Zhu 10/1/97                                               */
/*                                                                       */
/* Fill in the relation matrix information.                              */
/*************************************************************************/

#include <string.h>
#include "common.h"
#include "matrix.h"
#include "matrix_dna.h"
#include "Sankoff.h"

void fill_Sank_rel_matrix(Model B, Sank_S SK)
{
   double ***pp_related;
   float  ***ff_related;
   char   **matrix_name;
   int    i,j,k,n;
   int    Nmatrix,matrix_num,AlphaLen;

   pp_related=SK->rel_matrix.pp_related;
   ff_related=SK->rel_matrix.ff_related;
   matrix_name=SK->rel_matrix.matrix_name;
   Nmatrix=SK->rel_matrix.Nmatrix;
   AlphaLen=B->IP->nAlphaLen;

   if(SK->flags.protein_sequence) { /* for protein */
       for(n=0; n<Nmatrix;  n++) {
           if(Nmatrix>1){/* multiple matrices */
	       if(SK->rel_matrix.matrix_num==MAX_NUM_MATRIX_B)
	           /* blosum matrix */
	           matrix_num=n;
	        else
	           matrix_num=n+MAX_NUM_MATRIX_B;
	   }
	   else {
	        matrix_num=SK->rel_matrix.matrix_num;
	   }

	   /* fill the relation matrices */
	   for(i=0; i<AlphaLen; i++) {
	       for(j=0; j<AlphaLen; j++){
		   pp_related[n][i][j]=p_related[matrix_num][i][j];
		   ff_related[n][i][j]=f_related[matrix_num][i][j];
	       }
	   }
	   strcpy(matrix_name[n],name_related[matrix_num]);
       }
   }
   else { /* for DNA */
       for(k=0; k<SK->rel_matrix.Nmatrix; k++) {
	   matrix_num=SK->rel_matrix.index_start+k*SK->rel_matrix.index_step-1;
	   for(i=0; i<AlphaLen; i++) {
	       for(j=0; j<AlphaLen; j++){
		   pp_related[k][i][j]=p_related_dna[matrix_num][i][j];
		   ff_related[k][i][j]=f_related_dna[matrix_num][i][j];
	       }
	   }
	   sprintf(matrix_name[k],"PAM_%d",matrix_num+1);
       }
   }
}
