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

/*********************************************************************/
/*                                                                   */
/* $Id: /home/junzhu/proj/RCS/Gibbs_GUI.c,v 1.3 1997                 */
/*                                                                   */
/* Author: Jun Zhu                                                   */
/* Date:   July 1, 1997                                              */
/* Description: graphic user interface for Gibbs                     */
/*                                                                   */
/*********************************************************************/

#include <Xm/Form.h>
#include "routines.h"

void Gibbs_GUI(int argc, char **argv) {

  int            i;
  Widget         pw;                           /* paned widget for contain  */
                                               /* all other widgets         */
  Arg            args[20];                    /* argument list             */
  int            n;                            /* argument register         */
  char           string[80];
  char           colorname[80];
  char           ch;
  FILE           *fpt;
  XTextProperty  windowName,iconName;          /* window name and icon name */
  char           * icon_name ="Bayes Aligner";
  char           * window_name="Bayes Aligner";
  toplevel = XtAppInitialize(&app_context,"Bayes Aligner",NULL,0,
       &argc,argv,NULL,NULL,0);

  n = 0;
  XtSetArg(args[n], XmNwidth,  750);  n++; 
  XtSetArg(args[n], XmNheight, 750);  n++;
  XtSetValues(toplevel,args,n);

  /*  XtToolkitInitialize(); */
  /* create an application context */
  /*   app_context=XtCreateApplicationContext();

   if ((idispl = XtOpenDisplay(app_context, NULL, argv[0], 
              "Bayes Aligner", NULL, 0, & argc, argv)) == NULL)
    { XtWarning("*** Bayes Aligner can not open display! "); exit (-1); } 
      toplevel = XtAppCreateShell(argv[0], "Bayes", 
                applicationShellWidgetClass, idispl, NULL, 0); */

  /* XtGetApplicationResources(toplevel, (char *)&fontdata, resources,
             XtNumber(resources), NULL, 0);  */
  /* allocate space for GUI_argv */
  NEWP(GUI_argv,30,char);
  for(n=0;n<30;n++) {
     NEW(GUI_argv[n],50,char);
  }


  /* pop up the windows */
    XtRealizeWidget(toplevel);

   idispl=XtDisplay(toplevel); 
  /* find default screen pointer */
   if((iscreen = DefaultScreenOfDisplay(idispl)) == 0)
    {XtWarning("*** Gibbs  can not determine screen!"); exit (-1);} 

  /* find depth of screen (number of bit planes) */
  /* if((idepth = DefaultDepthOfScreen(iscreen)) < 8)
    {XtWarning("*** Screen depth is too shallow!"); exit (-1); } */

  /* find root window pointer */
  /*  iroot=RootWindowOfScreen(iscreen);  */

  /* allocate color map */
  /* mapdef=DefaultColormap(idispl,DefaultScreen(idispl)); */
  /* allocate the cell for default colors defined in X11 directory */
  /* if((fpt=fopen(rgb_file,"r"))!=NULL) {
     i=0;
     while((ch=getc(fpt))!=EOF) {
        string[i]=ch;
	while((ch=getc(fpt))!='\n') {
	   i++;
	   string[i]=ch;
	}
	i=0;
	sscanf(string,"%hd %hd %hd %s", &Colors[0].red,
	       &Colors[0].green,&Colors[0].blue,colorname);	   
	XAllocNamedColor(idispl,mapdef,colorname,&Colors[0],&Colors[1]);
     }
     fclose(fpt);       
  }
  else {
     printf("Could not open \"%s\"!\n",rgb_file);
  }
  */
  /* initiate individual color cell */
  /*  for(i=0;i<=MAXCOLOR-1;i++) {
     Colors[i].pixel=i;
     Colors[i].flags=DoRed|DoGreen|DoBlue;
  }
  XQueryColors(idispl,mapdef,Colors,MAXCOLOR);*/
      
  /* set WM property for window's name and icon's name */
   iwtop=XtWindow(toplevel);
  XStringListToTextProperty(&icon_name,1,&iconName);
  XStringListToTextProperty(&window_name,1,&windowName);
  XSetWMName(idispl,iwtop,&windowName);
  XSetWMIconName(idispl,iwtop,&iconName); 
  GUI_Sankoff(toplevel,NULL, NULL);
  XtAppMainLoop(app_context);
}

void quit() {
   exit(0);
}




