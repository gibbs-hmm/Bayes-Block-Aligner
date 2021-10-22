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

#ifndef  ROUTINES_H
#define  ROUTINES_H

#include "common.h"
#include "Sankoff.h"

int  Gibbs_NON_GUI(int,char**);


void fill_Sank_rel_matrix(Model,Sank_S);
void Sankoff(Model, Sank_S);
void Sankoff_Back_Sampling(Model,Sank_S);
void Sankoff_Back_Sampling_Min(Model,Sank_S);
void Sankoff_Back_Sampling_cutoff(Model,Sank_S);
void Sankoff_Bayesian_align(Model B, Sank_S SK);
void Sankoff_Make_S(Stype, Sank_S);
int  Sankoff_get_string(Model B, Sank_S SK, int option);
void Sankoff_probability_N(Model,Sank_S);
void Sankoff_probability_W(Model,Sank_S);
void Sankoff_probability_BW(Model,Sank_S);
void Sankoff_Trace_S(Model, Sank_S);
double SC(Model,Sank_S,int);
float max3(float,float,float);


#endif

#ifdef GUI

#include "GUI_common.h"
void Gibbs_GUI(int, char**);
void GUI_Bernoulli();
void GUI_Help(Widget, void*, void*);
void GUI_Sankoff(Widget, void*, void*);
void GUI_Sankoff_Search();
void GUI_win(Widget);
void Exit(Widget,void*,void*);

#else

int GUI_argc;   /* non gui version conditional */
char **GUI_argv; 

#endif

