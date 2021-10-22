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
/*                                                                        */
/* DESCRIPTION : This function reads in the sequences from the designated */
/*               input file.  The sequences are to be in fast-A format,   */
/*               meaning that there is a one line description for each    */
/*               sequence beginning with the greater than symbol (>).  It */
/*               is assumed that the elements can be stored in either     */
/*               upper case or lower case, and there can be white space   */
/*               between them.                                            */
/*               option=0, read query sequence and save it in Seq[0]      */
/*                      1, read data sequence and save it in Seq[1]       */
/**************************************************************************/

#include "common.h"
#include "Sankoff.h"


int Sankoff_get_string(Model B, Sank_S SK, int option)
{
   int loc;         /* postion in comment or in sequence   */
   char tmpch;      /* to hold the character read from file*/
   int  case_change;/* flag to indicate whether case change*/
   char tmpmsg[80]; /* temperory message.                  */
   int seqsize=50;  /* sequence size with default value 50 */
   int seq_ptr;     /* sequence indicator (0/1)            */
                    /* 0: query, 1: database               */
   FILE *fptr;      /* file pointer                        */
   Stype Seq=B->Seq;/* for quick access                    */


   if(option==0){ /* read query sequence */
       NEW(Seq->SeqLen, 2, int);
       NEWP(Seq->Orig,2,short);
       NEW(Seq->Orig[0], seqsize, short);
       Seq->Orig[1]=NULL;
       fptr=SK->queryfile_ptr;
       seq_ptr=0;
   }
   else{/* read data base sequences each at a time */
       if(Seq->Orig[1]!=NULL) free(Seq->Orig[1]);
       NEW(Seq->Orig[1], seqsize, short);
       fptr=B->IP->Datafiles->fpt;
       seq_ptr=1;
   }

   loc=0;
   tmpch=fgetc(fptr);
   while(tmpch != '>')
       tmpch=fgetc(fptr);
   if(tmpch == '>'){
       /** Once we see the '>', we know that we have reached  **/
       /** the beginning of the other sequence.               **/
       SK->comment[seq_ptr][loc]=tmpch;
       tmpch=fgetc(fptr);
       for(;tmpch!='\n';tmpch=fgetc(fptr)){
	   if(loc<48){
	     loc++;
	     SK->comment[seq_ptr][loc]=tmpch;
	   }
       }
       loc++;
       /* add an end to the string */
       SK->comment[seq_ptr][loc]=(char) 0;
   }

   loc=0;
   while(tmpch!= EOF && tmpch!='>') {
           case_change=0;
	   /** The next if statement allows us to ignore spaces, newlines, **/
	   /** and tabs.  If we want to ignore anything else, add it here  **/
	   if((tmpch != ' ') && (tmpch != '\n') && (tmpch != '\t')) {
	        if((int)tmpch < 97) {               /* Convert to lowercase */
		    tmpch = (char)((int)tmpch + 32);/* letters              */
		    case_change=1;
		}
                if(B->IP->nAlphaLen == 5) {/* option for DNA */
		    switch(tmpch) {
		    case 'a': break;
		    case 'c': break;
	            case 'n':               /* allow low complexity */
	            case 'x':                     /* regions              */
                      tmpch = 'e'; break;
	 	    case 't':
		      tmpch = 'b'; break;
		    case 'g':
		      tmpch = 'd'; break;
		    default:
		      fprintf(stderr,"\a\a\a%s\n", SK->comment[seq_ptr]);
		      fprintf(stderr, "%c is not a valid Nucleic Acid\n",
			     (char)((int)tmpch-case_change*32));
		      break;
		    }
		}
		else {   /* option for protein */
		    switch(tmpch) {      /* this will allow us to convert */
		    case 'v' :           /* to add into separate buckets  */
		      tmpch = 'b'; break;
		    case 'w' :
		      tmpch = 'j'; break;
		    case 'y' :
		      tmpch = 'o';  break;
		    case 'x':/* accept X as input */
		      tmpch = 'u';  break;
		    default :
		      break;
		    }
		}
		if((int)(tmpch) < 97 || (int)(tmpch) > 117){
		    fprintf(stderr,"%s\n", SK->comment[seq_ptr]);
		    sprintf(tmpmsg,
			  "%c is not a valid amino acid abbreviation",
			 (char)((int) tmpch-case_change*32));
		    fprintf(stderr,"%s\n",tmpmsg);
		    tmpch='x';

		}
		Seq->Orig[seq_ptr][loc] = (int)tmpch-97;
		loc++;

		/* if the string is going to be longer than the space      */
		/* allocated for it, go ahead and allocate some more space */
		if((loc+2)==seqsize){
		    seqsize+=50;
		    Seq->Orig[seq_ptr]=(short *)realloc(Seq->Orig[seq_ptr],
						   (seqsize)*sizeof(short));
		}
	   }
	   tmpch=fgetc(fptr);
   }
   Seq->SeqLen[seq_ptr]=loc;
   if(tmpch==EOF)
       return FALSE;
   else{
       fseek(fptr,-1L,SEEK_CUR); /* put ">" back to stream */
       return TRUE;
   }
}


