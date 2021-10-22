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
/* $Id: $                                                          */
/* Author: Jun Zhu                                                 */
/* Date: 1996/08/09                                                */
/* Description: Subroutines used by GUI                            */
/*                                                                 */
/*******************************************************************/

#include "GUI_common.h"
#include <Xm/BulletinB.h>
#include <Xm/CascadeBG.h>
#include <Xm/MessageB.h>
#include <Xm/LabelG.h>
#include <Xm/PushB.h>
#include <Xm/RowColumn.h>
#include <Xm/Scale.h>
#include <Xm/Text.h>
#include <Xm/ToggleBG.h>


/********************   GUI_strip_space  ***************************/
/*                                                                 */
/*  strip off space at both ends                                   */
/*******************************************************************/

void GUI_strip_space (char *context)
{
      register int i,j,L;   

      L=strlen(context);
      j=0;
      for(i=0;i<=L-1;i++) {
	  if(context[i]!=' ') {
	      context[j]=context[i];
	      j++;
	  }
      }
      /* add end deliminator */
      context[j]='\0';      
}


/********************   message  ***********************************/
/*                                                                 */
/*   prints out a message in message Widget                        */
/*******************************************************************/

void message (char *context)
{
     static XmTextPosition igo  = 0; 
     static XmTextPosition iend = 0;
     /* insert the text after the last position */

     XmTextInsert(iw_message,iend, context);

     /* add a newline to the string */
     /*    XmTextInsert(iw_message, XmTextGetLastPosition(iw_message), 
                       "\n"); */

     /* update next line starting position */ 
     igo = XmTextGetLastPosition(iw_message);

     /* remember position of end of this line */
     iend = igo;

}

/************************  wid_error   ******************************/ 
/*                                                                  */
/*          create an error widget                                  */
/********************************************************************/

Widget wid_error(Widget iw_parent, Widget iw_error, char * message, 
                   int ix, int iy)
{ 
 Arg           args[10];         /* arg list */
 register      int  n;           /* arg count */
 Position      ixp,iyp;          /* parent widget's position */

 /* need to find location relative to present position of top window */
 XtSetArg(args[0], XmNx, &ixp);
 XtSetArg(args[1], XmNy, &iyp);
 XtGetValues(toplevel,args,2);

 if (ix < 0 || iy < 0) {
    /*  use default placement for dialog box */
    ix = 30;
    iy = 30;
 }

 n = 0;
 XtSetArg(args[n], XmNmessageString,
	  XmStringCreateLtoR(message,XmSTRING_DEFAULT_CHARSET)); n++;
 XtSetArg(args[n], XmNx, (Position) ix+ixp); n++;
 XtSetArg(args[n], XmNy, (Position) iy+iyp); n++;
 XtSetArg(args[n], XmNdefaultPosition, False); n++;
 XtSetArg(args[n], XmNtitle, "Input Error"); n++;

 /* ring the bell */
 printf("\a\a\a\n");

 if (iw_error == (Widget)0) {
    /* create widget first */
    iw_error = XmCreateErrorDialog(iw_parent,"iw_error",
                                            args,n);
    XtManageChild(iw_error);

    /* get rid of help and cnacel buttons in widget*/
    XtUnmanageChild(XmMessageBoxGetChild(
                  iw_error,XmDIALOG_HELP_BUTTON));
    XtUnmanageChild(XmMessageBoxGetChild(
                  iw_error,XmDIALOG_CANCEL_BUTTON));

 }
 else { 
    XtSetValues(iw_error,args,n);
 }
 return iw_error;
}
 
/************************  wid_dialog  ******************************/ 
/*                                                                  */
/*          create a dialog widget                                  */
/********************************************************************/

Widget wid_dialog(Widget iw_parent, Widget iw_dialog, char * title, 
                   int ix, int iy)
{ 
 Arg           args[10];         /* arg list */
 register      int  n;           /* arg count */
 Position      ixp,iyp;          /* parent widget's position */

 /* need to find location relative to present position of top window */
 XtSetArg(args[0], XmNx, &ixp);
 XtSetArg(args[1], XmNy, &iyp);
 XtGetValues(toplevel,args,2);

 if (ix < 0 || iy < 0) {
    /*  use default placement for dialog box */
    ix = 30;
    iy = 30;
 }

 n = 0;
 XtSetArg(args[n], XmNx, (Position) ix+ixp); n++;
 XtSetArg(args[n], XmNy, (Position) iy+iyp); n++;
 XtSetArg(args[n], XmNdefaultPosition, False); n++;
 XtSetArg(args[n], XmNresizePolicy, True); n++;
 XtSetArg(args[n], XmNtitle, title); n++;
  XtSetArg(args[n], XmNwidth,  750);  n++; 
  XtSetArg(args[n], XmNheight, 750);  n++;


 if (iw_dialog == (Widget)0) {
    /* create widget first */
   /*  iw_dialog = XmCreateBulletinBoardDialog(iw_parent,title,
                                            args,n);*/
    iw_dialog = XmCreateBulletinBoard(iw_parent,title,
                                            args,n);
    if (iw_dialog == (Widget)0)
    {
       message("*** WID_DIALOG Cannot open dialog widget.");
       return iw_dialog;
    }
 }
 else { 
    XtSetValues(iw_dialog,args,n);
 }
 return iw_dialog;
}
 

/*******************  wid_labelg  **************************************/
/*                                                                     */
/*     creal a label                                                   */
/*                                                                     */
/***********************************************************************/
Widget wid_labelg(Widget iw_parent, Widget iw_labelg, 
               char *label, int ix, int iy)
{ 
 Arg           args[3];          /* arg list */
 int           n;                /* arg count */
 XmString      str_label;        /* String label */

 n = 0;

 str_label = XmStringCreateLtoR(label,XmSTRING_DEFAULT_CHARSET);
 XtSetArg(args[n], XmNlabelString, str_label); n++;
 
 if (iw_labelg <= (Widget)0)
 {     /* create new label gadget */
       if (ix > 0 && iy > 0)
       {   /* set label position */
           XtSetArg(args[n], XmNx, (Position) ix); n++;
           XtSetArg(args[n], XmNy, (Position) iy); n++;
       }
       else
       {   /* use default alignment */
           XtSetArg(args[n], XmNalignment, XmALIGNMENT_BEGINNING);  n++;
       }

       iw_labelg = XmCreateLabelGadget(iw_parent,"iw_labelg",args,n);
 }
 else  {
       XtSetValues(iw_labelg,args,n);
 }

 /* free the String (not doing this may result in a memory leak) */
 XmStringFree(str_label);
 XtManageChild(iw_labelg);
 return iw_labelg;

 }
 

/*******************   wid_pulldown   *******************************/
/*                                                                  */
/*   creates pull-down menu                                         */
/*                                                                  */
/********************************************************************/

Widget wid_pulldown(Widget iw_parent, char *label, char quick)
{ 
 Arg           args[10];         /* arg list  */
 register      int  n;           /* arg count */
 Widget        iw_pull, iw_button;

 /*  create menu pane */

 iw_pull = XmCreatePulldownMenu(iw_parent,"iw_menup",NULL,0);

 n = 0;
 XtSetArg(args[n], XmNsubMenuId, iw_pull); n++;
 XtSetArg(args[n], XmNlabelString, XmStringCreate(label,
                                 XmSTRING_DEFAULT_CHARSET)); n++;
 XtSetArg(args[n], XmNmnemonic , quick); n++;
 /* create menu cascade */
 iw_button = XmCreateCascadeButtonGadget(iw_parent,"iw_cascade",args,n);
 XtManageChild(iw_button);

 return iw_pull;

}
 

/***********************  wid_pushg ********************************/
/*                                                                 */
/*                                                                 */
/*   creates a push button widget                                  */
/*******************************************************************/

Widget wid_pushg(Widget iw_parent, Widget iw_push, 
               char *label, void (*cb)(Widget, void*, void*),
               char *data,  int ix, int iy)
{ 

 Arg           args[10];   /* arg list  */
 register      int  n;     /* arg count */

 n = 0;

 XtSetArg(args[n], XmNx, ix); n++;
 XtSetArg(args[n], XmNy, iy); n++;
 XtSetArg(args[n], XmNlabelString, XmStringCreate(label,
                              XmSTRING_DEFAULT_CHARSET)); n++;

 if (iw_push <= (Widget)0) {
    /* create widget first */
    iw_push = XmCreatePushButton(iw_parent,"iw_pb",args,n);
    XtAddCallback(iw_push,XmNactivateCallback,cb,data);
 }
 else {
    XtRemoveAllCallbacks(iw_push,XmNactivateCallback); 
    XtSetValues(iw_push,args,n);
    XtAddCallback(iw_push,XmNactivateCallback,cb,data);
 }

 XtManageChild(iw_push);

 return iw_push;
}



/********************  wid_rowcol  ***************************************/
/*                                                                       */
/*            create a Row_Column widget                                 */
/*                                                                       */ 
/*************************************************************************/
 
Widget wid_rowcol(Widget iw_parent, char type, int ix, int iy) 
{ 
 Arg           args[10];         /* arg list */
 register      int  n;           /* arg count */
 Position      ixp,iyp;          /* parent widget's position */
 Widget        iw_rowcol;

 n = 0;
 if (ix >= 0 || iy >= 0)
 {
    /* need to find location relative to present position of top window */
    XtSetArg(args[0], XmNx, &ixp);
    XtSetArg(args[1], XmNy, &iyp);
    XtGetValues(toplevel,args,2);

    XtSetArg(args[n], XmNx, (Position) ix+ixp); n++;
    XtSetArg(args[n], XmNy, (Position) iy+iyp); n++;
    XtSetArg(args[n], XmNdefaultPosition, False); n++;
 }

 if (type == 'h')
 {
    XtSetArg(args[n], XmNorientation, XmHORIZONTAL);    n++;
 }
 else
 {
    XtSetArg(args[n], XmNspacing, 6);    n++;
 }

 /* create a RowColumn widget as a parent for text input box  */
 iw_rowcol = XmCreateRowColumn(iw_parent,"iw_rowcol", args, n);

 if (iw_rowcol == (Widget)0)
 {
    message("*** wid_rowcol cannot create widget.");
    return iw_rowcol;
 }

 XtManageChild(iw_rowcol);

 return iw_rowcol;
 
}

/********************  wid_scale  ***************************************/
/*                                                                      */
/*          create a scale bar                                          */
/*                                                                      */
/************************************************************************/

Widget wid_scale(Widget iw_parent, Widget iw_its, char * label,
              int imin, int imax, int inow,  int iwid, int ihi)

{   
 Arg          args[10];    /* arg list  */
 register int n;           /* argument count */
 Widget       iw_lab;
 Widget       iw_rowcolh;

 if (iw_its == (void *)0) {

    /* must create new widget  */
    if (strlen(label) > 0) {
       /* need rowcol widget to hold label and scale widget */
       iw_rowcolh  = wid_rowcol(iw_parent,'h',-1,-1);

       /* create label widget*/
       iw_lab      = wid_labelg(iw_rowcolh, 0, label, -1, -1);
    }
    else
       iw_rowcolh = iw_parent; 

    /* create scale widget */
    n=0;
    XtSetArg(args[n], XmNwidth,         iwid);           n++;        
    XtSetArg(args[n], XmNheight,        ihi);            n++; 
    XtSetArg(args[n], XmNscaleWidth,    iwid-10);        n++;         
    XtSetArg(args[n], XmNscaleHeight,   ihi/2);          n++;     
    XtSetArg(args[n], XmNminimum,       imin);           n++;            
    XtSetArg(args[n], XmNmaximum,       imax);           n++;         
    XtSetArg(args[n], XmNdecimalPoints, 0);              n++;       
    XtSetArg(args[n], XmNvalue,         inow);           n++;        
    XtSetArg(args[n], XmNorientation,   XmHORIZONTAL);   n++;       
    XtSetArg(args[n], XmNshowValue,     TRUE);           n++;
 
    iw_its = XmCreateScale(iw_rowcolh, "iw_scale", args, n);
 }
 else {
    /* alter scale size in existing widget */

    XtSetArg(args[0], XmNminimum, imin);       
    XtSetArg(args[1], XmNmaximum, imax);       
    XtSetArg(args[2], XmNvalue,   inow);
       
    XtSetValues(iw_its, args, 3);
 }

 XtManageChild(iw_its);
 return iw_its;

}


/*****************  wid_textboxb  ****************************************/ 
/*                                                                       */
/*      create a text input                                              */
/*                                                                       */
/*************************************************************************/

Widget wid_textboxb(Widget iw_parent, Widget iw_text,  
               char *prompt, char *string,  int icol)

{  
 Arg           args[10];    /* arg list  */
 register      int n;       /* arg count */
 Widget        iw_rowcol;

 n = 0;
 
 if (iw_text == (Widget)0)
    {
    if (strlen(prompt) > (size_t) 0)
       {    
       /* create a new RowColumn widget as a parent for text input box  */
       XtSetArg(args[n], XmNorientation, XmHORIZONTAL);    n++;
       iw_rowcol = XmCreateRowColumn(iw_parent, "iw_rowcol", args, n);
       XtManageChild(iw_rowcol);

       /* create label widget for the prompt */
       wid_labelg(iw_rowcol, 0, prompt, -1,-1);
       }
    else
       iw_rowcol = iw_parent;

    n = 0;
    XtSetArg(args[n], XmNvalue, string); n++;                  
    XtSetArg(args[n], XmNcolumns, (short) icol); n++; 
    XtSetArg(args[n], XmNeditMode, XmSINGLE_LINE_EDIT); n++; 
 
    /* XmTextField is not used because it results in a memory leak */ 
    iw_text = XmCreateText(iw_rowcol, "iw_text", args, n);
    }

 else 
    {    /* alter existing textbox string value */
    XtSetArg(args[n], XmNvalue, string); n++;                  
    XtSetArg(args[n], XmNcolumns, (short) icol); n++; 
    XtSetValues(iw_text, args, n);

    /* note that iw_text is child of iw_rowcol which may be unmanaged */
    XtManageChild(XtParent(iw_text));
    }

 XtManageChild(iw_text);
 return (iw_text);
 }


/*******************  wid_toggleg  *************************************/
/*                                                                     */
/*       create a toggle button                                        */
/*                                                                     */
/***********************************************************************/

Widget wid_toggleg(Widget iw_parent, Widget iw_toggleg, 
		   char *label, int value,
		   void (*cb)(Widget, void*, void*), 
		   char *data, int ix, int iy)
{ 
 Arg           args[10];          /* arg list */
 register      int  n;           /* arg count */
 XmString      str_temp;

 n = 0;

 if (ix > 0 || iy > 0)
    {
    XtSetArg(args[n], XmNx, (Position) ix); n++;
    XtSetArg(args[n], XmNy, (Position) iy); n++;
    }

 /* set the color (red) */
 XtSetArg(args[n], XmNselectColor, 30); n++;
        
 str_temp =  XmStringCreate(label, XmSTRING_DEFAULT_CHARSET); 
 XtSetArg(args[n], XmNlabelString,  str_temp); n++;

 if ( iw_toggleg == (Widget) 0)
    {   
    /* create widget first time */
    iw_toggleg = XmCreateToggleButtonGadget(iw_parent,"iw_tog",args,n);
    }

 else 
    {   /* update callbacks in case they are changed */
    XtRemoveAllCallbacks(iw_toggleg, XmNvalueChangedCallback);

    /* set label & position */ 
    XtSetValues(iw_toggleg,args,n);
    }

 /* can free the string now */
 XmStringFree(str_temp);

 /* add call back */
 XtAddCallback(iw_toggleg,XmNvalueChangedCallback,cb,data);

 /* set the toggle, do not notify other toggles */
 XmToggleButtonGadgetSetState(iw_toggleg,value,FALSE);

 XtManageChild(iw_toggleg);
 return iw_toggleg;
 }
 


























