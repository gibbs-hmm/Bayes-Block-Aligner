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

#ifndef GUI_COMMON_H
#define GUI_COMMON_H

/* include files for all GUI subroutines */
#include <X11/Xlib.h>
#include <X11/StringDefs.h>
#include <X11/Intrinsic.h> 
#include <Xm/Xm.h>

/* common parameters for GUI */
 Display        *idispl;           /* current display         */
 Window         iwtop;             /* top window identifier */
 Window         iroot;             /* root window identifier */
 Screen         *iscreen;          /* current screen          */ 
 Widget         toplevel;          /* top widget              */
 Widget         menubar;           /* menu bar widget         */
 Widget         canvas;            /* Draw Area widget    */
 Widget         dialog_box;        /* progress box            */ 
 Widget         iw_file;           /* file selection widget   */      
 Widget         iw_message;        /* message Widget */
 XtAppContext   app_context;       /* application context   */
 int            idepth;            /* depth of current screen */
 int            iposix,iposiy;     /* current position of pixmap      */
 GC             icontx;            /* default graphic context */
 GC canvas_gc;                     /* GC of the drawing area */

 typedef struct{
       XFontStruct *fontptr;  /* XFont pointer */
       int font_height, font_width; 
 } font_data, *font_data_ptr;
 font_data fontdata; 

static XtResource resources[] = {
  { XtNfont, XtCFont, XmRFontStruct, sizeof(XFontStruct *),
    XtOffset(font_data_ptr, fontptr),XmRString, "Fixed"}
};

/* command line data from GUI */
int GUI_argc;
char **GUI_argv; 

#define MAXCOLOR 256 
XColor       Colors[MAXCOLOR];   /* color cells */
Colormap     mapdef;             /* default color map */ 
static String  rgb_file="/usr/openwin/lib/X11/rgb.txt"; 

/* subroutine used by GUI */
void GUI_strip_space(char *);
void message (char *);
Widget wid_error(Widget, Widget, char*, int, int);
Widget wid_dialog(Widget, Widget, char*, int, int);
Widget wid_labelg(Widget, Widget, char*, int, int);
Widget wid_pulldown(Widget, char*, char);
Widget wid_pushg(Widget, Widget, char*, XtCallbackProc,char*,  int, int);
Widget wid_rowcol  (Widget, char, int, int);
Widget wid_scale   (Widget, Widget, char *, int, int, int,  int, int);
Widget wid_textboxb(Widget, Widget, char*, char*,  int);
Widget wid_toggleg (Widget, Widget, char*, int, XtCallbackProc, 
		    char*, int, int);


#endif



