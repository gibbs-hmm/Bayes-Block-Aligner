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

/**************************************************************/
/* $Id: Gibbs.c,v 1.2 1996/11/12 16:54:43 junzhu Exp junzhu $ */
/*                                                            */
/* Author:     Jun Zhu, June 25,1996                          */
/*                                                            */
/**************************************************************/

#include <stdint.h>
#include "routines.h"

int main(int argc, char **argv) {
   /* check if we want GUI */
   if( argc >= 2) {
      Gibbs_NON_GUI(argc,argv);
   }
#ifdef GUI
   else { /* start GUI */
      Gibbs_GUI(argc,argv);               
   }
#else
   else {
     print_usage_sankoff();
   }
#endif
 return(0);
   
}
