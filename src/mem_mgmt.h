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

/* FILE NAME : mem_mgmt.h */

#ifndef MEM_MGMT_H
#define MEM_MGMT_H
#include "common.h"
#include "Sankoff.h"

void init_maxdata(MaxResults *tmp);
void free_maxdata(MaxResults *tmp, IPtype IP);
Model alloc_Model();
void alloc_Counts(Model B);
void alloc_Sank_BW(Model,Sank_S);
void alloc_Sank(Model,Sank_S);
void alloc_Sank_S(Model,Sank_S);
void alloc_Sank_N(Model,Sank_S);
void alloc_Sank_N(Model,Sank_S);
void alloc_Sank_rel_matrix(Model,Sank_S);
void alloc_Sank_W(Model,Sank_S);

void free_Sank_BW(Model,Sank_S);
void free_Sank(Model,Sank_S);
void free_Sank_S(Model,Sank_S);
void free_Sank_N(Model,Sank_S);
void free_sank_N(Model,Sank_S);
void free_Sank_rel_matrix(Model,Sank_S);
void free_Sank_W(Model,Sank_S);
#endif

