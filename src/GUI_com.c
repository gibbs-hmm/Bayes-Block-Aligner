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
/* $Id: $                                                      */
/*                                                                 */
/* Author:  Jun Zhu                                                */
/*  Date:  July 1, 1996                                            */
/* Description: Create menu buttons for Gibbs                      /
/*                                                                 */
/*******************************************************************/

#include <stdlib.h>
#include "GUI_common.h"
#include "routines.h"
 
void GUI_com(Widget menubar)
{
 char         *buttons[]={"Commands", "Help"}; /* contents of menubar  */
 char         menu_attri;                      /* attribute for button */
 Widget       pulldown[2];                     /* pulldown buttons     */
 Widget       pbutton;                         /* buttons in each      */
                                               /* pulldown menupane    */

 /*set button for Command */
 menu_attri='C';
 pulldown[0]=wid_pulldown(menubar, buttons[0],menu_attri);;

 /* create menu for Command */
 /*pbutton=wid_pushg(pulldown[0],(Widget)0,"Bernoulli",
                       GUI_Bernoulli,"11",0,0);*/
 pbutton=wid_pushg(pulldown[0],(Widget)0,"Bayes_align",
		   GUI_Sankoff,"12",0,0);
 pbutton=wid_pushg(pulldown[0],(Widget)0,"Exit",Exit,"13",0,0);


 /*set button for Help */
 menu_attri='H';
 pulldown[1]=wid_pulldown(menubar, buttons[1],menu_attri);;

 /* set Help button at the end */

 /* create menu for Help */
 pbutton=wid_pushg(pulldown[1],(Widget)0,"For Bernoulli",
                       GUI_Help,"21",0,0);
 pbutton=wid_pushg(pulldown[1],(Widget)0,"For Sankoff",
                       GUI_Help,"22",0,0);

 /* display all the buttons */ 
 XtManageChildren(pulldown,XtNumber(buttons));

}


/* callback function for help menu */
void GUI_Help(Widget iw_temp, void* data, void* call_data)
{ 
  
   if(strcmp((char *)data,"21")==0)     /* help for bernoulli */
      system("netscape /home/junzhu1/public_html/manual/bernoulli\
/bernoulli.html &");
   
}

/* callback function for exit */
void Exit(Widget iw_temp, void* data, void* call_data)
{   XtDestroyWidget(toplevel);
    exit(0);   
}
