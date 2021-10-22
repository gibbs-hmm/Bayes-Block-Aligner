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

/***********************************************************************/
/* $Id: $                                                              */
/*                                                                     */
/* Author: Jun Zhu                                                     */
/* Date:  July 1, 1996                                                 */
/* Description: Create child widgets of Panewindow                     */
/*               All the widgets in the main window are created        */
/*               menubar(XmMenuBar):                                   */
/*	            at the top of main window to contain all the menu  */
/*                  and submenu                                        */
/*               iw_message(XmScrolledText):                           */
/*                  at the bottom of the main window to give dialog or */
/*                  warning messages                                   */
/*		 canvas(XmDrawingArea):                                */
/*                  at the center of the main window. It is the working*/
/* 		    area for drawing data flow graph                   */
/*                                                                     */
/***********************************************************************/

#include "GUI_common.h"
#include <Xm/RowColumn.h>
#include <Xm/Text.h>
#include <Xm/ScrollBar.h>
#include <Xm/DrawingA.h>

/* routines used */
void GUI_com(Widget);
void GUI_Sankoff(Widget iw_temp, void* data, void* call_data);

void GUI_win(Widget pw) 
{

  Widget       menubar;   /* menubar to contain all the buttons */
  Arg          wargs[20]; /* argument list */
  int          n;         /* argument register */
   int          data, call_data; 
 
  /* create a menubar to contain buttons */
  n = 0;
  XtSetArg(wargs[n], XmNtopAttachment,  XmATTACH_FORM);  n++;
  XtSetArg(wargs[n], XmNleftAttachment, XmATTACH_FORM);  n++;
  XtSetArg(wargs[n], XmNrightAttachment,XmATTACH_FORM);  n++;
  /* menubar =XmCreateMenuBar(pw, "menubar", wargs, n); */

  /* create menu buttons */
  /*  GUI_com(menubar); 

  XtManageChild(menubar);*/

  /* create a message Widget */
  /*  n=0;
  XtSetArg(wargs[n], XmNscrollingPolicy,    XmAUTOMATIC);       n++;
  XtSetArg(wargs[n], XmNscrollBarDisplayPolicy, XmSTATIC);      n++;
  XtSetArg(wargs[n], XmNscrollBarPlacement, XmBOTTOM_RIGHT);    n++;
  XtSetArg(wargs[n], XmNvisualPolicy,       XmCONSTANT);        n++;
  XtSetArg(wargs[n], XmNeditMode,           XmMULTI_LINE_EDIT); n++;
  XtSetArg(wargs[n], XmNeditable,           False);             n++;
  XtSetArg(wargs[n], XmNleftAttachment,     XmATTACH_FORM);     n++;
  XtSetArg(wargs[n], XmNrightAttachment,    XmATTACH_FORM);     n++; 
  XtSetArg(wargs[n], XmNbottomAttachment,   XmATTACH_FORM);     n++;
  XtSetArg(wargs[n], XmNheight,             90);                n++;
  XtSetArg(wargs[n], XmNtraversalOn,   False);        n++; 
  iw_message = XmCreateScrolledText(pw, "message", wargs, n);
  XtManageChild(iw_message); */

  /* print out first message */
  /* message("Gibbs>\n"); */

  /* Create a widget in which to display image. */
     n = 0;
       /* XtSetArg(wargs[n], XmNtopWidget,        menubar);            n++;*/
  XtSetArg(wargs[n], XmNtopAttachment,    XmATTACH_WIDGET);    n++;
  /* XtSetArg(wargs[n], XmNbottomWidget,     iw_message);         n++; */
  XtSetArg(wargs[n], XmNbottomAttachment, XmATTACH_WIDGET);    n++;
  XtSetArg(wargs[n], XmNleftAttachment,   XmATTACH_FORM);      n++;
  XtSetArg(wargs[n], XmNrightAttachment,  XmATTACH_FORM);      n++;

  /* canvas = XtCreateWidget("canvas", xmDrawingAreaWidgetClass, 
                         pw, wargs, n); */
  /* XtManageChild(canvas); */
  /* GUI_Sankoff(pw,data, call_data);*/

}
